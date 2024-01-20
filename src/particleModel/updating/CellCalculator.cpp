#include "CellCalculator.h"
#include "particleModel/storage/CellContainerIterators.h"
#include "utils/ArrayUtils.h"
#include <iostream>
#include <spdlog/spdlog.h>
#include <limits>
#include <list>
#include <algorithm>
#include <cmath>

double min_distance = 0.7;

std::vector<std::vector<double>> sigma_mixed{{1.0}};

std::vector<std::vector<double>> epsilon_mixed{{5.0}};

CellCalculator::CellCalculator(CellContainer &cellContainer, double delta_t, double cutoff,double r_l_,
      std::array<boundary_conditions,6> boundaries_cond,std::string forceType,double gravity_factor)
    : cellContainer(cellContainer), gravity_factor(gravity_factor), delta_t(delta_t), cutoff(cutoff), r_l(r_l_),
    domain_bounds(cellContainer.getDomainBounds()), domain_max_dim(CellContainer::domain_max_dim), boundaries(boundaries_cond),
    particles(CellContainer::particles)
    {
    ghost_reflection_is_off = true;
    for(auto b : boundaries) {
        if(b == boundary_conditions::ghost_reflective) ghost_reflection_is_off = false;
    }

    if(forceType == "smoothedLJ"){
        force = forceSmoothedLennJonesPotentialFunction(sigma_mixed,epsilon_mixed,cutoff,r_l);
    }else if(forceType == "LJ"){
        force = forceLennJonesPotentialFunction(sigma_mixed,epsilon_mixed,cutoff);
    }else if(forceType == "Gravity"){
        force = forceSimpleGravitational(cutoff);
    }else{
        throw std::invalid_argument("Force Type was not correctly specified");
    }

}

void CellCalculator::calculateX(){
    instructions cell_updates;

    //todo omp here
    for (auto cell = begin_CI(); cell != end_CI(); ++cell) {


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
    //todo omp here
    for (auto cell = begin_CI(); cell != end_CI(); ++cell) {
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
    calculateLinkedCellF();
    calculatePeriodicF();
    //todo omp here
    for (auto iter = begin_CI(); iter != end_CI(); ++iter) {
        finishF(&(*iter));
    }
}

void CellCalculator::shiftF(){
    for (auto cell = begin_CI(); cell != end_CI(); ++cell) {
        for (auto particle_ptr: *cell) {
            particle_ptr->shiftF();
        }
    }
}


void CellCalculator::calculateLinkedCellF() {
    static std::array<dim_t, 3> current_cell;
    static std::array<dim_t, 3> pattern;
    std::vector<Particle*>* cell_1;
    std::vector<Particle*>* cell_2;
    std::array<double, 3> F_ij{};

    //write new path inside current/start and pattern
    cellContainer.setNextPath(current_cell, pattern);

    //todo make start & pattern iterator
    while(current_cell[0] != dim_t_res) {
        cell_1 = &particles[current_cell[0]][current_cell[1]][current_cell[2]];
        current_cell[0] += pattern[0];
        current_cell[1] += pattern[1];
        current_cell[2] += pattern[2];

        while(0 < current_cell[0] && 0 < current_cell[1] && 0 < current_cell[2] && current_cell[0] <= domain_max_dim[0]
              && current_cell[1] <= domain_max_dim[1] && current_cell[2] <= domain_max_dim[2]) {

            cell_2 = &particles[current_cell[0]][current_cell[1]][current_cell[2]];

            for(auto & p_i : *cell_1) {
                for(auto & p_j : *cell_2) {

                    if(inCutoffDistance(*p_i, *p_j, {0, 0, 0})) {
                        F_ij = force(*p_i, *p_j, {0, 0, 0});

                        for (int i = 0; i < 3; i++) {
                            p_i->addF(i, F_ij[i]);
                            p_j->addF(i, -F_ij[i]);
                        }
                    }
                }
            }

            cell_1 = &particles[current_cell[0]][current_cell[1]][current_cell[2]];
            current_cell[0] += pattern[0];
            current_cell[1] += pattern[1];
            current_cell[2] += pattern[2];
        }

        cellContainer.setNextPath(current_cell, pattern);
    }
}

void CellCalculator::calculatePeriodicF() {
    for(std::array<dim_t,3> pattern : CellContainer::patterns) {

        //todo omp here
        for (PeriodIterator it = begin_PI(pattern); it != end_PI(); ++it) {

            std::array<dim_t, 3> current_cell = *it;
            std::vector<Particle *> *cell_1 = &particles[current_cell[0] - pattern[0]]
            [current_cell[1] - pattern[1]]
            [current_cell[2] - pattern[2]];

            std::array<double, 3> particle_offset{0, 0, 0};
            std::array<double, 3> F_ij{};

            if (mirror(current_cell, particle_offset)) {

                std::vector<Particle *> *cell_2 = &particles[current_cell[0]][current_cell[1]][current_cell[2]];

                for (auto &p_i: *cell_1) {
                    for (auto &p_j: *cell_2) {

                        if (inCutoffDistance(*p_i, *p_j, particle_offset)) {
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
            }
        }
    }

    cellContainer.allocateCellFromPosition(x, new_cell_position);

    if( 0 < new_cell_position[0] &&  new_cell_position[0] <= domain_max_dim[0] &&
        0 < new_cell_position[1] &&  new_cell_position[1] <= domain_max_dim[1] &&
        0 < new_cell_position[2] && new_cell_position[2] <= domain_max_dim[2]) {

        //todo lock here
        cell_updates.emplace_back(particle_ptr, new_cell_position);

    }else{
        SPDLOG_INFO("new halo particle: " + (*particle_ptr).toString());
        //todo lock here
        //just delete them because we don't have halo particles anyway
        //cellContainer.getHaloParticles().push_back(particle_ptr);
        auto& instances = cellContainer.particle_instances;
        auto it = std::find(instances.begin(),instances.end(),*particle_ptr);
        if(it != instances.end())
            instances.erase(it);
    }
}

void CellCalculator::updateCells(instructions& cell_updates) {
    //todo omp here vs simply sequential?
    for(auto ins : cell_updates){

      Particle* particle_ptr = std::get<0>(ins);
      std::array<dim_t, 3>& new_cell_position = std::get<1>(ins);

      std::vector<Particle*> *new_cell = &particles[new_cell_position[0]][new_cell_position[1]][new_cell_position[2]];
      //todo lock here
      new_cell->push_back(particle_ptr);
    }
}

//mirror the last position back into the domain, return true if all boundaries successful
bool CellCalculator::mirror(std::array<dim_t,3> &position, std::array<double,3> &offset) {
    static std::array<unsigned short,6> map_boundaries{3,5,1,2,4,0};//{neg_X, neg_Y, neg_Z, pos_X, pos_Y, pos_Z}
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


void CellCalculator::finishF(std::vector<Particle*> *current_cell) {
    Particle* p_i;
    Particle* p_j;
    std::array<double, 3> F_ij{};

    for (auto it1 = current_cell->begin(); it1 != current_cell->end(); ++it1) {
        for(auto it2 = std::next(it1); it2 != current_cell->end(); ++it2) {
            p_i = *it1;
            p_j = *it2;

            if(inCutoffDistance(*p_i, *p_j, {0, 0, 0})) {
                F_ij = force(*p_i, *p_j, {0, 0, 0});

                for (int i = 0; i < 3; i++) {
                    p_i->addF(i, F_ij[i]);
                    p_j->addF(i, -F_ij[i]);
                }
            }
        }

        //add gravity
        (*it1)->addF(1, (*it1)->getM() * gravity_factor);
    }
}

bool CellCalculator::inCutoffDistance(Particle &p1, Particle &p2, const std::array<double,3> &offset) const {
    static double compare_distance = cutoff * cutoff;
    const auto& x1 = p1.getX(), x2 = p2.getX();

    double dx = x1[0] - x2[0] + offset[0];
    double dy = x1[1] - x2[1] + offset[1];
    double dz = x1[2] - x2[2] + offset[2];

    return dx * dx + dy * dy + dz * dz <= compare_distance;
}



void CellCalculator::initializeFX() {
    std::array<dim_t, 3> current_position;
    instructions cell_updates;

    //calculate F
    calculateLinkedCellF();

    //write new coordinates in current_position array
    cellContainer.setNextCell(current_position);

    while(current_position[0] != dim_t_res) {
        std::vector<Particle*> *current_cell = &particles[current_position[0]][current_position[1]][current_position[2]];

        //finish F calculation
        finishF(current_cell);

        for (auto iter = current_cell->begin(); iter != current_cell->end();) {
            Particle& particle = *(*iter);

            calculateVX(particle, false);
            particle.shiftF();

            std::array<dim_t, 3> position;
            cellContainer.allocateCellFromPosition(particle.getX(), position);
            
            if (position[0] != current_position[0] ||
                position[1] != current_position[1] ||
                position[2] != current_position[2])
            {
                cell_updates.emplace_back(*iter,position);
                iter = current_cell->erase(iter);

            } else {
                iter++;
            }
        }

        cellContainer.setNextCell(current_position);
    }

    //update particle distribution in the cells
    updateCells(cell_updates);
}



void CellCalculator::calculateWithinFVX() {
    std::array<dim_t, 3> current_position;
    instructions cell_updates;

    //write new coordinates in current_position array
    cellContainer.setNextCell(current_position);

    while(current_position[0] != dim_t_res) {
        std::vector<Particle*> *current_cell = &particles[current_position[0]][current_position[1]][current_position[2]];
        //finish F calculation
        finishF(current_cell);

        for (auto iter = current_cell->begin(); iter != current_cell->end();) {
            Particle* particle_ptr = *iter;
            Particle& particle = *particle_ptr;

            calculateVX(particle, true);
            particle.shiftF();

            std::array<dim_t, 3> position;
            cellContainer.allocateCellFromPosition(particle.getX(), position);

            if (position[0] != current_position[0] ||
                position[1] != current_position[1] ||
                position[2] != current_position[2])
            {
                cell_updates.emplace_back(*iter,position);
                iter = current_cell->erase(iter);
            } else{
                iter++;
            }
        }
        cellContainer.setNextCell(current_position);
    }
    updateCells(cell_updates);
}



void CellCalculator::calculateVX(Particle &particle, bool calculateV) {
    const std::array<double, 3> &f = particle.getF();
    const std::array<double, 3> &x = particle.getX();
    const std::array<double, 3> &v = particle.getV();
    const double &m = particle.getM();

    if(calculateV) {
        const std::array<double, 3> &f_old = particle.getOldF();

        //calculate V
        for (int i = 0; i < 3; i++) {
            particle.setV(i, v[i] + delta_t * (f[i] + f_old[i]) / (2 * m));
        }
    }

    //calculate X
    double x_0 = x[0] + delta_t * v[0] + delta_t * delta_t * f[0] / 2.0 / m;
    double x_1 = x[1] + delta_t * v[1] + delta_t * delta_t * f[1] / 2.0 / m;
    double x_2 = x[2] + delta_t * v[2] + delta_t * delta_t * f[2] / 2.0 / m;

    particle.setX(0, x_0);
    particle.setX(1, x_1);
    particle.setX(2, x_2);
}

