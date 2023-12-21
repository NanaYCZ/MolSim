#include <gtest/gtest.h>
#include <cmath>
#include "particleModel/updating/CellCalculator.h"
#include "particleModel/storage/CellContainerIterators.h"

/**
 * @file contains the old methods for calculating Ghost particle boundary conditions 
 *       for validating the new one
*/

/**
 * @brief iterates of the Top or the Bottom side of the domain and calculates forces
 *        for boundary ghost particles
 * 
 * 
 * This function iterates over cells of cellContainer within the domain.
 * Depending on lower_z and upper_z iterates either over Top or Bottom side of domain.
 * iterates for a fixed z over the plane 
 * {
 *  1(minimum domain cell in y direction) until domain_max_dim[1](maximum domain cell in y direction),
 *  1(minimum domain cell in x direction) until domain_max_dim[0](maximum domain cell in x direction)
 * }
 * (this plane corresponds to one layer of cells -> iterates over all these cells)
 * iterates over a plane of this form for every z in [lower_z,upper_z]  (both are inclusive borders)
 * 
 * -> for every plane (horizontal plane of cells) all Particles in the plane of cells are iterated 
 * and it is checked if they are closer than sigma * 2^(1/6) to
 * 'z_border' (either lower domain border or upper domain border). 
 * If a Particle p1 is close enough, forces are calculated between this Particle p1 and a Particle, that 
 * has a position that corresponds to the position of p1, but mirrored at z_border.
 *  
 * @param lower_z lower bound for the x-y-planes of cells that are iterated
 * @param upper_z upper bound for the x-y-planes of cells that are iterated
 * @param z_border real valued domain boundary, forces are calculated relative to this border
 * 
 * 
*/
void calculateBoundariesTopOrBottom(dim_t lower_z,dim_t upper_z, double z_border,
                             CellCalculator& cellCalculator ){

    auto domain_max_dim = cellCalculator.getDomain_Max();
    auto particles = cellCalculator.getParticles();
    dim_t x, y;
    x = y =  1;

  // bottom boundary
  while(lower_z <= upper_z){
  while (y < domain_max_dim[1]) {
    std::vector<Particle*>& cell = particles[x][y][lower_z];

    for (auto particle_pointer : cell) {
      Particle& particle = *particle_pointer;
      double x_dim = particle.getX()[0];
      double y_dim = particle.getX()[1];
      double z_dim = particle.getX()[2];

      // a assume that we have an offset of 1 everywhere
      double distance =z_dim - z_border;
      double ref_size = std::pow(2,1.0/6.0) * sigma_mixed[particle.getType()][particle.getType()];
      if (std::abs(distance) < ref_size) {
        // calculate repulsing force with Halo particle
        double ghost_particle_z = z_dim - 2*distance;

        std::array<double,3> F_particle_ghost = cellCalculator.ghostParticleLennardJonesForce(particle,{x_dim,y_dim,ghost_particle_z});
        particle.addF(0, F_particle_ghost[0]);
        particle.addF(1, F_particle_ghost[1]);
        particle.addF(2, F_particle_ghost[2]);
      }
    }

    x++;
    if (x >= domain_max_dim[0]) {
      x = 1;
      y++;
    }
  }
  x = 1;
  y = 1;
  lower_z++;
  }

};

/**
 * @brief iterates of the Front or the Back side of the domain and calculates forces
 *        for boundary ghost particles
 * 
 * This function iterates over cells of cellContainer within the domain.
 * Depending on lower_x and upper_x iterates either over Front or Back side of domain.
 * iterates for a fixed x over the plane 
 * {
 *  1(minimum domain cell in y direction) until domain_max_dim[1](maximum domain cell in y direction),
 *  1(minimum domain cell in z direction) until domain_max_dim[2](maximum domain cell in z direction)
 * }
 * (this plane corresponds to one layer of cells -> iterates over all these cells)
 * iterates over a plane of this form for every x in [lower_x,upper_x]  (both are inclusive borders)
 * 
 * -> for every plane (vertical plane of cells) all Particles in the plane of cells are iterated 
 * and it is checked if they are closer than sigma * 2^(1/6) to
 * 'x_border' (either lower domain border or upper domain border). 
 * If a Particle p1 is close enough, forces are calculated between this Particle p1 and a Particle, that 
 * has a position that corresponds to the position of p1, but mirrored at x_border.
 *  
 * Additionaly for 2D boundaries have the parameter z_until, that ensures, that not the whole
 * plane 
 * {
 *  1(minimum domain cell in y direction) until domain_max_dim[1](maximum domain cell in y direction),
 *  1(minimum domain cell in z direction) until domain_max_dim[2](maximum domain cell in z direction)
 * }
 * is iterated, but only 
 * {
 *  1(minimum domain cell in y direction) until domain_max_dim[1](maximum domain cell in y direction),
 *  1 until z_until (z_until will be 1 for 2D Boundary Conditions)
 * }
 * 
 * @param lower_x lower bound for the y-z-planes of cells that are iterated
 * @param upper_x upper bound for the y-z-planes of cells that are iterated
 * @param x_border real valued domain boundary, forces are calculated relative to this border
 * @param z_until determines how much of the plane will be calculated
 * 
 * 
*/
void calculateBoundariesFrontOrBack(dim_t lower_x,dim_t upper_x ,double x_border, dim_t z_until,
                                                    CellCalculator& cellCalculator){

    auto domain_max_dim = cellCalculator.getDomain_Max();
    auto particles = cellCalculator.getParticles();
    

    dim_t y, z;
    z = y =  1;

  // bottom boundary
  while(lower_x <= upper_x){
  while (z <= z_until) {
    std::vector<Particle*>& cell = particles[lower_x][y][z];

    for (auto particle_pointer : cell) {
      Particle& particle = *particle_pointer;
      double x_dim = particle.getX()[0];
      double y_dim = particle.getX()[1];
      double z_dim = particle.getX()[2];

      // a assume that we have an offset of 1 everywhere
        double distance = x_dim - x_border;
      double ref_size = std::pow(2,1.0/6.0) * sigma_mixed[particle.getType()][particle.getType()];
      if (std::abs(distance) < ref_size) {
        // calculate repulsing force with Halo particle
        double ghost_particle_x = x_dim - 2 * distance;

        std::array<double,3> F_particle_ghost = cellCalculator.ghostParticleLennardJonesForce(particle,{ghost_particle_x,y_dim,z_dim});
        particle.addF(0, F_particle_ghost[0]);
        particle.addF(1, F_particle_ghost[1]);
        particle.addF(2, F_particle_ghost[2]);
      }
    }

    y++;
    if (y >= domain_max_dim[1]) {
      y = 1;
      z++;
    }
  }
  y = 1;
  z = 1;
  lower_x++;
  }
}; //Front and Back



/**
 * @brief iterates of the Left or the Right side of the domain and calculates forces
 *        for boundary ghost particles
 * 
* 
* This function iterates over cells of cellContainer within the domain.
* Depending on lower_y and upper_y iterates either over Left  or Right  side of domain.
* iterates for a fixed y over the plane 
* {
*  1(minimum domain cell in x direction) until domain_may_dim[0](maximum domain cell in x direction),
*  1(minimum domain cell in z direction) until domain_max_dim[2](maximum domain cell in z direction)
* }
* (this plane corresponds to one layer of cells -> iterates over all these cells)
* iterates over a plane of this form for every y in [lower_y,upper_y]  (both are inclusive borders)
* 
* -> for every plane (also vertical plane of cells) all Particles in the plane of cells are iterated 
* and it is checked if they are closer than sigma * 2^(1/6) to
* 'y_border' (will be either lower domain border or upper domain border). 
* If a Particle p1 is close enough, forces are calculated between this Particle p1 and a Particle, that 
* has a position that corresponds to the position of p1, but mirrored at y_border.
*  
* Additionaly for 2D boundaries have the parameter z_until, that ensures, that not the whole
* plane 
* {
*  1(minimum domain cell in x direction) until domain_may_dim[0](maximum domain cell in x direction),
*  1(minimum domain cell in z direction) until domain_max_dim[2](maximum domain cell in z direction)
* }
* is iterated, but only 
* {
*  1(minimum domain cell in x direction) until domain_max_dim[0](maximum domain cell in x direction),
*  1 until z_until (z_until will be 1 for 2D Boundary Conditions)
* }
* 
* @param lower_y lower bound for the x-z-planes of cells that are iterated
* @param upper_y upper bound for the x-z-planes of cells that are iterated
* @param y_border real valued domain boundary, forces are calculated relative to this border
* @param z_until determines how much of the plane will be calculated
* 
* 
*/
void calculateBoundariesLeftOrRight(dim_t lower_y,dim_t upper_y ,double y_border, dim_t z_until,
                                                    CellCalculator& cellCalculator){
    
    auto domain_max_dim = cellCalculator.getDomain_Max();
    auto particles = cellCalculator.getParticles();
    
    dim_t x, z;
    z = x =  1;

  // bottom boundary
  while(lower_y <= upper_y){
  while (z <= z_until) {
    std::vector<Particle*>& cell = particles[x][lower_y][z];

    for (auto particle_pointer : cell) {
      Particle& particle = *particle_pointer;
      double x_dim = particle.getX()[0];
      double y_dim = particle.getX()[1];
      double z_dim = particle.getX()[2];

      // a assume that we have an offset of 1 everywhere
        double distance = y_dim - y_border;
      double ref_size = std::pow(2,1.0/6.0) * sigma_mixed[particle.getType()][particle.getType()];
      if (std::abs(distance) < ref_size) {
        // calculate repulsing force with Halo particle
        double ghost_particle_y = y_dim - 2 * distance;

        std::array<double,3> F_particle_ghost = cellCalculator.ghostParticleLennardJonesForce(particle,{x_dim,ghost_particle_y,z_dim});
        particle.addF(0, F_particle_ghost[0]);
        particle.addF(1, F_particle_ghost[1]);
        particle.addF(2, F_particle_ghost[2]);
      }
    }

    x++;
    if (x >= domain_max_dim[0]) {
      x = 1;
      z++;
    }
  }
  x = 1;
  z = 1;
  lower_y++;
  }
}; //Left and Right





/**
 * @brief Uses two Containers and two Calculators and fills them with exactly the same particles.
 *        Those particlse are particles in every boundary region. Then the new and the old ghost
 *        particle boundary conditions are applied. The new ghost particle boundary conditions
 *        to the first Container (container) and the old ones to the second (other_container).
 *        Then it is compared, if that particles in the previously exact same containers
 *        are still exactly the same (meaning the both boundary conditions applied the same forces)
 * 
 * 
*/
TEST(test_new_Boundaries,test_basic){
   CellContainer cellContainer(4,4,1,1.0,1.0);

    CellCalculator cellCalculator(cellContainer,0.0014,1.0,
            {boundary_conditions::reflective,boundary_conditions::reflective,
        boundary_conditions::reflective,boundary_conditions::reflective,
        boundary_conditions::reflective,boundary_conditions::reflective});

    //adding Partciles in every relevant region
    cellContainer.addParticle({0.1,0.1,0},{0,0,0},1); //corner
    cellContainer.addParticle({2,0.1,0},{0,0,0},1); //x-border in negative direction in the middle 
    cellContainer.addParticle({2,3.99,0},{0,0,0},1); //x-border in positive direction in the middle 
    cellContainer.addParticle({3.99,2,0},{0,0,0},1); //y-border in positive direction in the middle 
    cellContainer.addParticle({0.1,2,0},{0,0,0},1); //y-border in negative direction in the middle

    cellContainer.addParticle({2,2,0},{0,0,0},1); //particle not at the boundary (further away then 1.1225 from boundary)
    //all of their forces will be initalized to 0

    cellContainer.createPointers();

    /*
    setup a second equivalent cellContainer
    
    */


    CellContainer other_cellContainer(4,4,1,1.0,1.0);


    CellCalculator other_cellCalculator(other_cellContainer,0.0014,1.0,
            {boundary_conditions::reflective,boundary_conditions::reflective,
        boundary_conditions::reflective,boundary_conditions::reflective,
        boundary_conditions::reflective,boundary_conditions::reflective});


     //adding Partciles in every relevant region
    other_cellContainer.addParticle({0.1,0.1,0},{0,0,0},1); //corner
    other_cellContainer.addParticle({2,0.1,0},{0,0,0},1); //x-border in negative direction in the middle 
    other_cellContainer.addParticle({2,3.99,0},{0,0,0},1); //x-border in positive direction in the middle 
    other_cellContainer.addParticle({3.99,2,0},{0,0,0},1); //y-border in positive direction in the middle 
    other_cellContainer.addParticle({0.1,2,0},{0,0,0},1); //y-border in negative direction in the middle

    other_cellContainer.addParticle({2,2,0},{0,0,0},1); //particle not at the boundary (further away then 1.1225 from boundary)
    //all of their forces will be initalized to 0

    other_cellContainer.createPointers();

    std::cout << cellContainer.to_string() << std::endl;

    std::cout << other_cellContainer.to_string() << std::endl;

    cellCalculator.applyReflectiveBoundaries(); //now Boundary Conditions should have been applied to all particles

    // now do the same manually with the old functions
    // we know that reflective boundaries are applied on every side in 2D
    auto domain_max_dim = other_cellCalculator.getDomain_Max();
    auto domain_bounds = other_cellCalculator.getDomainBounds();
    auto comparing_depth = other_cellContainer.getComparingdepth();
    dim_t z_max = 1;

    calculateBoundariesFrontOrBack(domain_max_dim[0]-comparing_depth,domain_max_dim[0],
                                   domain_bounds[0],z_max,other_cellCalculator);

    calculateBoundariesFrontOrBack(1,comparing_depth+1,0,z_max,other_cellCalculator);

    calculateBoundariesLeftOrRight(domain_max_dim[1]-comparing_depth,domain_max_dim[1],
                                    domain_bounds[1],z_max,other_cellCalculator); 

    calculateBoundariesLeftOrRight(1,comparing_depth+1,0,z_max,other_cellCalculator); 

    auto iter = cellContainer.begin();
    auto other_iter = cellContainer.begin();

    for(; iter != cellContainer.end() && other_iter != other_cellContainer.end(); ++iter, ++other_iter){
        auto particle_iter = (*iter).begin();
        auto other_particle_iter = (*other_iter).begin();
        for(;particle_iter != (*iter).end() && other_particle_iter != (*other_iter).end();
                                                 ++particle_iter , ++other_particle_iter){
            Particle& particle = **particle_iter;
            Particle& other_particle = **other_particle_iter;
            std::cout << "Comparing particles: " << std::endl;
            std::cout << particle.toString() << std::endl;
            std::cout << other_particle.toString() << std::endl;
            //check that all particles are equal
            ASSERT_EQ(particle,other_particle);
        }
    }





}





