#include "CellCalculator.h"
#include "CellContainerIterators.h"
#include "utils/ArrayUtils.h"
#include <iostream>
#include <spdlog/spdlog.h>
#include <list>
#include <algorithm>

int chunk_size = 16;

std::vector<std::vector<double>> sigma_mixed{{1.0}};
std::vector<std::vector<double>> epsilon_mixed{{1.0}};

CellCalculator::CellCalculator(CellContainer &cellContainer, double delta_t, double cutoff,
                               double r_l_, std::array<boundary_conditions,6> boundaries_cond,
                               force_type forceType,double gravity_factor, concurrency_strategy strategy)
    : parallelization(strategy), cellContainer(cellContainer), gravity_factor(gravity_factor), delta_t(delta_t), cutoff(cutoff), r_l(r_l_),
    domain_bounds(cellContainer.getDomainBounds()), domain_max_dim(CellContainer::domain_max_dim), boundaries(boundaries_cond),
    particles(CellContainer::particles)
    {
    if(forceType == force_type::smoothedLJ){
        force = forceSmoothedLennJonesPotentialFunction(sigma_mixed,epsilon_mixed,cutoff,r_l);
    }else if(forceType == force_type::LJ){
        force = forceLennJonesPotentialFunction(sigma_mixed,epsilon_mixed,cutoff);
    }else if(forceType == force_type::gravitational){
        force = forceSimpleGravitational(cutoff);
    }else if(forceType == force_type::Membrane){
        force = forceHarmonicForce(sigma_mixed,epsilon_mixed,cutoff);
    }else{
        throw std::invalid_argument("Force Type was not correctly specified");
    }
}

void CellCalculator::calculateX(){
    instructions cell_updates;

    #pragma omp parallel for default(none) shared(cell_updates, chunk_size) \
                            schedule(static,chunk_size) \
                            if(parallelization == concurrency_strategy::first_method)

    for (auto cell = begin_CellIterator(); cell != end_CellIterator(); ++cell) {

        //iterate trough particles of cell
        for (auto particle_ptr = (*cell).begin(); particle_ptr != (*cell).end();) {
            // dereferencing iterator yields a pointer to a particle, therefore dereference again
            Particle &particle = *(*particle_ptr);
            std::array<double, 3> old_f, v, old_x, new_x;

            old_f = particle.getOldF();
            v = particle.getV();
            old_x = particle.getX();

            double mass = particle.getM();

            new_x = old_x + delta_t * v + (delta_t * delta_t / (2 * mass)) * old_f;

            particle.setX(new_x);

            std::array<dim_t, 3> new_cell;
            cellContainer.allocateCellFromPosition(new_x, new_cell);

            if (new_cell[0] != cell.x || new_cell[1] != cell.y || new_cell[2] != cell.z) {

                applyBoundaries(*particle_ptr, new_cell, cell_updates);
                particle_ptr = (*cell).erase(particle_ptr);
            } else {
                particle_ptr++;
            }
        }
    }

    updateCells(cell_updates);
}



void CellCalculator::calculateV(){

    #pragma omp parallel for default(none) shared(chunk_size) \
                            schedule(static,chunk_size) \
                            if(parallelization == concurrency_strategy::first_method)

    for (auto cell = begin_CellIterator(); cell != end_CellIterator(); ++cell) {
        for (auto particle_ptr: *cell) {
            Particle &particle = *particle_ptr;
            std::array<double, 3> old_f, f, old_v;

            old_f = particle.getOldF();
            f = particle.getF();
            old_v = particle.getV();

            double mass = particle.getM();

            auto new_v = old_v + (delta_t / (2 * mass)) * (f + old_f);
            particle.setV(new_v);
        }
    }
}

void CellCalculator::calculateF(){
    calculateInterCellF();
    calculatePeriodicF();
    #pragma omp parallel for default(none) schedule(dynamic) \
                            if(parallelization == concurrency_strategy::first_method)

    for (auto iterator = begin_CellIterator(); iterator != end_CellIterator(); ++iterator) {
        calculateFWithin(&(*iterator));
    }
}


void CellCalculator::shiftF(){

    #pragma omp parallel for default(none) shared(chunk_size) \
                            schedule(static,chunk_size) \
                            if(parallelization == concurrency_strategy::first_method)

    for (auto cell = begin_CellIterator(); cell != end_CellIterator(); ++cell) {
        for (auto particle_ptr: *cell) {
            particle_ptr->shiftF();
        }
    }
}

void CellCalculator::calculateInterCellF() {
    for(std::array<dim_t,3> pattern : CellContainer::patterns) {

        #pragma omp parallel for default(none) shared(pattern) \
                                schedule(dynamic) \
                            if(parallelization == concurrency_strategy::first_method)

        for (StartPointIterator it = begin_StartIterator(pattern); it != end_StartIterator(); ++it) {
            std::array<double, 3> F_ij{};
            std::vector<Particle*>* cell_1;
            std::vector<Particle*>* cell_2;

            std::array<dim_t, 3> current_cell = *it;

            cell_1 = &particles[current_cell[0]][current_cell[1]][current_cell[2]];
            current_cell[0] += pattern[0];
            current_cell[1] += pattern[1];
            current_cell[2] += pattern[2];

            while (0 < current_cell[0] && 0 < current_cell[1] && 0 < current_cell[2] &&
                   current_cell[0] <= domain_max_dim[0] &&
                   current_cell[1] <= domain_max_dim[1] &&
                   current_cell[2] <= domain_max_dim[2]) {

                cell_2 = &particles[current_cell[0]][current_cell[1]][current_cell[2]];

                for (auto &p_i: *cell_1) {
                    for (auto &p_j: *cell_2) {


                            F_ij = force(*p_i, *p_j, {0, 0, 0});

                            for (int i = 0; i < 3; i++) {
                                p_i->addF(i, F_ij[i]);
                                p_j->addF(i, -F_ij[i]);
                            }

                    }
                }

                cell_1 = cell_2;
                current_cell[0] += pattern[0];
                current_cell[1] += pattern[1];
                current_cell[2] += pattern[2];
            }
        }
    }
}

void CellCalculator::calculatePeriodicF() {
    for(std::array<dim_t,3> pattern : CellContainer::patterns) {

        #pragma omp parallel for default(none) shared(pattern, chunk_size) \
                                schedule(static,chunk_size) \
                            if(parallelization == concurrency_strategy::first_method)

        for (StartPointIterator it = begin_StartIterator(pattern); it != end_StartIterator(); ++it) {

            std::array<dim_t, 3> current_cell = it.outside();
            std::vector<Particle *> *cell_1 = &particles[current_cell[0] - pattern[0]]
                                                        [current_cell[1] - pattern[1]]
                                                        [current_cell[2] - pattern[2]];

            std::array<double, 3> particle_offset{0, 0, 0};
            std::array<double, 3> F_ij{};

            if (mirror(current_cell, particle_offset)) {

                std::vector<Particle *> *cell_2 = &particles[current_cell[0]][current_cell[1]][current_cell[2]];

                for (auto &p_i: *cell_1) {
                    for (auto &p_j: *cell_2) {


                            F_ij = force(*p_i, *p_j, particle_offset);

                            for (int i = 0; i < 3; i++) {
                                p_i->addF(i, F_ij[i]);
                                p_j->addF(i, -F_ij[i]);
                            }

                    }
                }
            }
        }
    }
}



void CellCalculator::applyBoundaries(Particle* particle_ptr, std::array<dim_t, 3>& new_cell_position, instructions& cell_updates) {
    //second method for reflective boundaries
    const std::array<double,3> &x = particle_ptr->getX();
    const std::array<double,3> &v = particle_ptr->getV();
    //this map is necessary, because in the boundaries member, the order is
    //{positive_z,negative_z,positive_x,negative_x,positive_y,negative_y}
    //but here
    //{positive_x,negative_y,negative_z,positive_x,positive_y,positive_z}
    // is wanted
    static std::array<unsigned short,6> map_boundaries{3,5,1,2,4,0};//{neg_X, neg_Y, neg_Z, pos_X, pos_Y, pos_Z}

    for (int i = 0; i < 3; ++i) {
        //check if position is outside the domain
        if(x[i] < 0) {
            //check negative boundaries: 0 -> neg_X, 1 -> neg_Y, 2 -> neg_Z
            if(boundaries[map_boundaries[i]] == boundary_conditions::reflective) {
                //apply reflection in pos direction
                particle_ptr->setX(i, -x[i]);
                particle_ptr->setV(i, -v[i]);
            } else if(boundaries[map_boundaries[i]] == boundary_conditions::periodic){
                //apply periodic in pos direction
                particle_ptr->addX(i,domain_bounds[i]);
                particle_ptr->decBoundariesCrossedI(i); //crossed boundary in negative i direction
            }
        } //check if position is outside the domain
        else if(domain_bounds[i] < x[i]) {
            //check positive boundaries: 3 -> pos_X, 4 -> pos_Y, 5 -> pos_Z
            if(boundaries[map_boundaries[i+3]] == boundary_conditions::reflective) {
                //apply reflection in neg direction
                particle_ptr->setX(i, 2 * domain_bounds[i] - x[i]);
                particle_ptr->setV(i, -v[i]);
            
            } else if(boundaries[map_boundaries[i+3]] == boundary_conditions::periodic) {
                //apply periodic in neg direction
                particle_ptr->addX(i,-domain_bounds[i]);
                particle_ptr->incBoundariesCrossedI(i); //crossed boundary in positive i direction
            }
        }
    }

    cellContainer.allocateCellFromPosition(x, new_cell_position);

    if( 0 < new_cell_position[0] &&  new_cell_position[0] <= domain_max_dim[0] &&
        0 < new_cell_position[1] &&  new_cell_position[1] <= domain_max_dim[1] &&
        0 < new_cell_position[2] && new_cell_position[2] <= domain_max_dim[2]) {

        #pragma omp critical
        {
            cell_updates.emplace_back(particle_ptr, new_cell_position);
        }

    }else{
        //outflow
        SPDLOG_INFO("new halo particle: " + (*particle_ptr).toString());
        //just delete them because we don't have halo particles anyway
        #pragma omp critical
        {
            auto &instances = cellContainer.particle_instances;
            auto found = std::find(instances.begin(), instances.end(), *particle_ptr);
            if (found != instances.end())
                instances.erase(found);
        }
    }
}

void CellCalculator::updateCells(instructions& cell_updates) {
    for(auto instruction : cell_updates){

        std::array<dim_t, 3> &new_cell_position = std::get<1>(instruction);

        std::vector<Particle *> *new_cell = &particles[new_cell_position[0]][new_cell_position[1]][new_cell_position[2]];

        new_cell->push_back(std::get<0>(instruction));
    }
}

//mirror the last position back into the domain, return true if all boundaries successful
inline bool CellCalculator::mirror(std::array<dim_t,3> &position, std::array<double,3> &offset) {
    static std::array<unsigned short,6> map_boundaries{3,5,1,2,4,0};//{neg_X, neg_Y, neg_Z, pos_X, pos_Y, pos_Z}
    //if all exits were on periodic boundaries
    bool mirrored_fully = true;
    //this map is necessary, because in the boundaries member, the order is
    //{positive_z,negative_z,positive_x,negative_x,positive_y,negative_y}
    //but here 
    //{positive_x,negative_y,negative_z,positive_x,positive_y,positive_z}
    // is wanted
    for (int i = 0; i < 3; ++i) {
        if(position[i] < 1) {
            if(boundaries[map_boundaries[i]] == boundary_conditions::periodic) {
                position[i] = position[i] + domain_max_dim[i];
                offset[i] = domain_bounds[i];

            } else {
                mirrored_fully = false;
            }
        }
        else if(domain_max_dim[i] < position[i]) {
            if(boundaries[map_boundaries[i+3]] == boundary_conditions::periodic) {
                position[i] = position[i] - domain_max_dim[i];
                offset[i] = -domain_bounds[i];

            } else {
                mirrored_fully = false;
            }
        }
    }
    return mirrored_fully;
}


void CellCalculator::calculateSpecialForce(){
    for (auto cell = begin_CellIterator(); cell != end_CellIterator(); ++cell) {
        for (auto particle_ptr: *cell) {
            Particle &particle = *particle_ptr;
            for (int i=0;i<cellContainer.getSpecialPositions().size();i++){
                if (particle.getGrid()[0]==cellContainer.getSpecialPositions()[i][0]&&
                    particle.getGrid()[1]==cellContainer.getSpecialPositions()[i][1]&&
                    particle.getGrid()[2]==cellContainer.getSpecialPositions()[i][2]){
                    particle_ptr->addF(cellContainer.getSpecialStrength()[i]);
                }
            }
        }
    }
}

void CellCalculator::calculateFWithin(std::vector<Particle*> *current_cell) {
    Particle* p_i;
    Particle* p_j;
    std::array<double, 3> F_ij{};

    for (auto it1 = current_cell->begin(); it1 != current_cell->end(); ++it1) {
        for(auto it2 = std::next(it1); it2 != current_cell->end(); ++it2) {
            p_i = *it1;
            p_j = *it2;


                F_ij = force(*p_i, *p_j, {0, 0, 0});

                for (int i = 0; i < 3; i++) {
                    p_i->addF(i, F_ij[i]);
                    p_j->addF(i, -F_ij[i]);
                }

        }

        //add gravity
        (*it1)->addF(2, (*it1)->getM() * gravity_factor);
    }
}


