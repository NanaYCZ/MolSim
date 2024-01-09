#include "CellCalculator.h"
#include "particleModel/storage/CellContainerIterators.h"
#include "utils/ArrayUtils.h"
#include <iostream>
#include <spdlog/spdlog.h>
#include <limits>

double min_distance = 0.7;

std::vector<std::vector<double>> sigma_mixed{{1.0}};

std::vector<std::vector<double>> epsilon_mixed{{5.0}};

CellCalculator::CellCalculator(CellContainer &cellContainer, double delta_t, double cutoff,
      std::array<boundary_conditions,6> boundaries_cond, std::optional<double> target_temp_param,
       std::optional<double> max_temp_diff_param,double gravity_factor)
    : cellContainer(cellContainer), gravity_factor(gravity_factor), delta_t(delta_t), cutoff(cutoff),
    domain_bounds(cellContainer.getDomainBounds()), domain_max_dim(CellContainer::domain_max_dim),
    target_temp(target_temp_param), max_temp_diff(max_temp_diff_param), boundaries(boundaries_cond),
    particles(CellContainer::particles)
    {
    ghost_reflection_is_off = true;
    for(auto b : boundaries) {
        if(b == boundary_conditions::ghost_reflective) ghost_reflection_is_off = false;
    }
}

void CellCalculator::calculateX(){
  instructions cell_updates;

  for(auto cell = begin_CI(); cell != end_CI();++cell){
    //iterate trough particles of cell
    for(auto particle_ptr = (*cell).begin(); particle_ptr != (*cell).end();){
        // dereferencing iterator yields a pointer to a particle, therefore dereference again
      Particle& particle = *(*particle_ptr);
      std::array<double,3>  old_f, v, old_x, new_x;

      old_f = particle.getOldF(); v = particle.getV(); old_x = particle.getX();
      double mass = particle.getM();

      new_x = old_x + delta_t * v + (delta_t * delta_t / (2 * mass)) * old_f;
      particle.setX(new_x);

      std::array<dim_t, 3> new_cell;
      cellContainer.allocateCellFromPosition(new_x,new_cell);

      if(new_cell[0] != cell.x || new_cell[1] != cell.y || new_cell[2] != cell.z){
          //todo lock
          cell_updates.emplace_back(*particle_ptr,new_cell);
          particle_ptr = (*cell).erase(particle_ptr);
      }else{
          particle_ptr++;
      }
    }
  }
  updateCells(cell_updates);
}



void CellCalculator::calculateV(){
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

  for(auto iter = begin_CI(); iter != end_CI();++iter){
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

        //apply force between the last two cells of the path, the cell_1 being
        //the last one in the domain and cell_2 being the mirrored position of
        //the previously out of domain one
        std::array<double,3> particle_offset{0,0,0};
        if(mirror(current_cell, particle_offset)) {

            //todo lock/delegate mirrored
            cell_2 = &particles[current_cell[0]][current_cell[1]][current_cell[2]];

            for(auto & p_i : *cell_1) {
                for(auto & p_j : *cell_2) {

                    if(inCutoffDistance(*p_i, *p_j, particle_offset)) {
                        F_ij = force(*p_i, *p_j, particle_offset);

                        for (int i = 0; i < 3; i++) {
                            p_i->addF(i, F_ij[i]);
                            p_j->addF(i, -F_ij[i]);
                        }
                    }
                }
            }
        }

        cellContainer.setNextPath(current_cell, pattern);
    }
}



void CellCalculator::updateCells(instructions& cell_updates) {
    for(auto ins : cell_updates){

      Particle* particle_ptr = std::get<0>(ins);
      std::array<dim_t, 3> new_cell_position = std::get<1>(ins);

  
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

          std::vector<Particle*> *new_cell = &particles[new_cell_position[0]][new_cell_position[1]][new_cell_position[2]];
          //todo lock here
          new_cell->push_back(particle_ptr);
      }else{
        SPDLOG_INFO("new halo particle: " + (*particle_ptr).toString());
        //todo lock here
        cellContainer.getHaloParticles().push_back(particle_ptr);
      }
      
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

std::array<double,3> CellCalculator::force(const Particle &p_i, const Particle &p_j, const std::array<double,3> &offset) {
    const auto& x_i = p_i.getX(), x_j = p_j.getX();

    double sigma = sigma_mixed[p_i.getType()][p_j.getType()];
    double epsilon = epsilon_mixed[p_i.getType()][p_j.getType()];

    double norm = ArrayUtils::L2Norm(x_i - x_j + offset);
    norm = std::max(min_distance, norm);

    double prefactor = (-24 * epsilon) / (std::pow(norm, 2));

    prefactor *= (std::pow(sigma / norm, 6) - 2 * std::pow(sigma / norm, 12));

    return prefactor * (x_i - x_j + offset);
}


std::array<double,3> CellCalculator::ghostParticleLennardJonesForce(const Particle &particle,std::array<double,3> ghost_position){
  auto x = particle.getX();

  double sigma = sigma_mixed[particle.getType()][particle.getType()];
  double epsilon = epsilon_mixed[particle.getType()][particle.getType()];

  double norm = ArrayUtils::L2Norm(x - ghost_position);
  norm = std::max(min_distance, norm);

  double prefactor = (-24 * epsilon) / (std::pow(norm, 2));

  prefactor *= (std::pow(sigma / norm, 6) - 2 * std::pow(sigma / norm, 12));

  return prefactor * (x - ghost_position);
}



void CellCalculator::applyThermostats(){
  double kinetic_energy = 0;
  size_t amt = 0;

  for(auto iter = begin_CI(); iter != end_CI(); ++iter){
    for(Particle* particle_ptr : *iter){
      const std::array<double,3> &v = particle_ptr->getV();
      double v_squared = v[0] * v[0]  + v[1] * v[1] + v[2] * v[2];
      double m = particle_ptr->getM();
      kinetic_energy += v_squared * m;
      amt++;
    }
  }
  double k_boltzman = 1;
  //calculate the current temperatur from the current kinetic energy in the system
  //assuming we only have two kinds of dimensions namely 2 or 3
  double current_temp = kinetic_energy/((cellContainer.hasThreeDimensions() ? 3 : 2) * amt * k_boltzman);
  double next_temp;
  if(target_temp.has_value())
    next_temp = target_temp.value();  
  else  
    throw std::invalid_argument("applyThermostats was called, altough target temp was not provided ");
  

  //if the temperatur diffference would be too big cap it
  double temp_diff = next_temp - current_temp;
  if(max_temp_diff.has_value() && std::abs(temp_diff) > max_temp_diff){
    next_temp = (std::signbit(temp_diff) ? -1 : 1) * max_temp_diff.value() + current_temp;
  }


  // the scaling factor to reach the target temperature
  // assuming this is never negative because we only calculate with kelvin
  double temp_scaling = sqrt(next_temp/current_temp);

  //apply scaling
  for(auto iter = begin_CI(); iter != end_CI(); ++iter){
      for(Particle* particle_ptr : *iter){
        std::array<double,3> v = particle_ptr->getV();
        particle_ptr->setV(temp_scaling * v);
      }
  }

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


/**
 * @brief add force from Ghost Particles for every particle in the cell
 * 
 * @param i corresonds to direction i = 0 means in x dir, i = 1 means in y dir, i = 2 means in z dir.
 *          Using an i that is not part of {0,1,2} is undefined behaviour.
*/
void CellCalculator::addGhostParticleForcesInDir_i(int i,double boundary,
                                                std::vector<Particle*>& cell){
  for(auto particle_ptr : cell){
    Particle& particle = *particle_ptr;
    std::array<double,3> x = particle.getX();
    double distance = x[i] - boundary;
    double ref_size = std::pow(2,1.0/6.0) * sigma_mixed[particle.getType()][particle.getType()];
    if (std::abs(distance) < ref_size) {
      // calculate repulsing force with Halo particle
      double ghost_particle_i = x[i] - 2 * distance;
      std::array<double,3> ghost_particle_x = x;
      ghost_particle_x[i] = ghost_particle_i;
      std::array<double,3> F_particle_ghost = ghostParticleLennardJonesForce(particle,ghost_particle_x);
      particle.addF(F_particle_ghost);
    }
  }
}

/*
function images coordinate system like this:

  z ^ 
    |  ^ x
    | / 
    |/--------> y

*/
void CellCalculator::applyReflectiveBoundaries() {
    if(ghost_reflection_is_off) return;

    dim_t z_max = cellContainer.hasThreeDimensions() ? domain_max_dim[2] : 1;
    dim_t comparing_depth = cellContainer.getComparingdepth();

    if (cellContainer.hasThreeDimensions()) {
        //TOP SIDE
        //boundaries[0] corresponds to boundary_conditions in positiveZ direction
        if (boundaries[0] == boundary_conditions::ghost_reflective) {
            auto iter = cellContainer.begin_custom(
                    1, domain_max_dim[0], // iteration cuboid in x dim
                    1, domain_max_dim[1], // iteration cuboid in y dim
                    domain_max_dim[2] - comparing_depth, domain_max_dim[2]);  // iteration cuboid in z dim
            for (; iter != cellContainer.end_custom(); ++iter) {
                addGhostParticleForcesInDir_i(2, domain_bounds[2], *iter);
            }
        }

        //BOTTOM SIDE
        //boundaries[1] corresponds to boundary_conditions in negativeZ direction
        if (boundaries[1] == boundary_conditions::ghost_reflective) {
            auto iter = cellContainer.begin_custom(
                    1, domain_max_dim[0], // iteration cuboid in x dim
                    1, domain_max_dim[1], // iteration cuboid in y dim
                    1, comparing_depth);  // iteration cuboid in z dim
            for (; iter != cellContainer.end_custom(); ++iter) {
                addGhostParticleForcesInDir_i(2, 0, *iter);
            }
        }
    }

    //BACK SIDE
    /**
     * custom iterator, iterates over the cells on this side:

        /-------------/
       / |###########/|
      /  |##########/#|
     |   |#########|##|
     |   |#########|##|
     |  /          | /
     | /           |/
     |-------------|

    */
    //boundaries[2] corresponds to boundary_conditions in positiveX direction
    if (boundaries[2] == boundary_conditions::ghost_reflective) {
        auto iter = cellContainer.begin_custom(
                domain_max_dim[0] - comparing_depth, domain_max_dim[0], // iteration cuboid in x dim
                1, domain_max_dim[1], // iteration cuboid in y dim
                1, z_max);  // iteration cuboid in z dim
        for (; iter != cellContainer.end_custom(); ++iter) {
            addGhostParticleForcesInDir_i(0, domain_bounds[0], *iter);
        }
    }

    //FRONT SIDE
    /**
     * custom iterator, iterates over the cells on this side:

        /-------------/
       / |           /|
      /_____________/ |
     |#############|  |
     |#############|__|
     |#############| /
     |#############|/
     |-------------|

    */
    //boundaries[3] corresponds to boundary_conditions in negativeX direction

    if (boundaries[3] == boundary_conditions::ghost_reflective) {
        auto iter = cellContainer.begin_custom(
                1, comparing_depth, // iteration cuboid in x dim
                1, domain_max_dim[1], // iteration cuboid in y dim
                1, z_max);  // iteration cuboid in z dim

        for (; iter != cellContainer.end_custom(); ++iter) {
            addGhostParticleForcesInDir_i(0, 0, *iter);
        }
    }

    //RIGHT SIDE / POSITIVE Y DIRECTION
    /**
     * custom iterator, iterates over the cells on this side:

        /-------------/
       / |           /|
      /  |          /#|
     |   |         |##|
     |   |         |##|
     |  /          |#/
     | /           |/
     |-------------|

    */
    //boundaries[4] corresponds to boundary_conditions in positiveY direction
    if (boundaries[4] == boundary_conditions::ghost_reflective) {
        auto iter = cellContainer.begin_custom(
                1, domain_max_dim[0], // iteration cuboid in x dim
                domain_max_dim[1] - comparing_depth, domain_max_dim[1], // iteration cuboid in y dim
                1, z_max);  // iteration cuboid in z dim

        for (; iter != cellContainer.end_custom(); ++iter) {
            addGhostParticleForcesInDir_i(1, domain_bounds[1], *iter);
        }
    }

    //LEFT SIDE / NEGATIVE Y DIRECTION
    /**
     * custom iterator, iterates over the cells on this side:

        /-------------/
       /#|           /|
      /##|          / |
     |###|         |  |
     |###|         |  |
     |##/          | /
     |#/           |/
     |-------------|

    */
    //boundaries[5] corresponds to boundary_conditions in negativeY direction
    if (boundaries[5] == boundary_conditions::ghost_reflective) {
        auto iter = cellContainer.begin_custom(
                1, domain_max_dim[0], // iteration cuboid in x dim
                1, comparing_depth, // iteration cuboid in y dim
                1, z_max);  // iteration cuboid in z dim

        for (; iter != cellContainer.end_custom(); ++iter) {
            addGhostParticleForcesInDir_i(1, 0, *iter);
        }
    }
}
