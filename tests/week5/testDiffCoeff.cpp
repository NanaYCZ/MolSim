#include <gtest/gtest.h>
#include "particleModel/storage/CellContainer.h"
#include "particleModel/updating/CellContainerIterators.h"

#include "particleModel/updating/CellCalculator.h"
#include "particleModel/updating/ThermoStats.h"
#include "utils/ArrayUtils.h"




/**
 * @brief Does the same as the previous case, only that the particle is in the corner
 *        or our 2D rectangle domain and should reappear on the opposite corner
 * 
 * 
 * 
*/
TEST(test_CellCalculation, test_periodic_corner){
    CellContainer container_corner{10, 10, 0, 2, 2};
    CellCalculator calculator(container_corner,0.001,2.0,1.9,
        {boundary_conditions::periodic,boundary_conditions::periodic,
        boundary_conditions::periodic,boundary_conditions::periodic,
        boundary_conditions::periodic,boundary_conditions::periodic},force_type::smoothedLJ);
    ThermoStats thermoStats(container_corner,0.001);


    //all particles will automatically be intialized with force={0,0,0}

    //particle p1 in middle of the side that is the negative X boundary
    container_corner.addParticle({0,0,0},{-10,-10,0},1);
    //will be particles[0]


    container_corner.createPointers();

    thermoStats.initDiffusionCoefficient();


    // as the particle is exactly at the boundary -> should move to 
    // the respective other sides within one step already
    calculator.calculateX();

    std::cout << container_corner.to_string() << std::endl;

    std::vector<Particle> particles = container_corner.getInstances();

    std::cout << particles[0].toString() <<  "\n";

    std::array<int,3> crossed_boundaries = particles[0].getBoundariesCrossed();
    ASSERT_EQ(-1,crossed_boundaries[0]);
    ASSERT_EQ(-1,crossed_boundaries[1]);

    double diffusion_coefficient_actual = thermoStats.getDiffusionCoefficient();

    double diffusion_coefficient_expected = (std::pow(-10 * 0.001,2) + std::pow(-10 * 0.001,2));
    
    ASSERT_NEAR(diffusion_coefficient_actual,diffusion_coefficient_expected,1e-7);
    
}



/**
 * @brief First a CellContainer an CellCalculator are created, then two particles
 *        are added such, that one particle p1 is on the middle of the left side of a 2D rectangle (domain)
 *        and the other one is on the middle of the lower side. Both particles have velocities
 *        in the direction of the boundaries, that they are next to. Then the calculateX() method
 *        is called once and the new positions of the particles that are right at the boundary
 *        are calculated. We then expect that the reappear on the respective other sides
 *        
 *        |--------|                        |---p2---|
 *        |        |                        |        |
 *        |p1      |  -should become->      |       p1
 *        |        |                        |        |
 *        |___p2___|                        |________|
*/  
TEST(test_DiffCoeff1, test_periodic_diff_Coeff_reappearing){
    CellContainer container{10, 10, 0, 2, 2};
    CellCalculator calculator(container,0.05,2.0,1.9,
        {boundary_conditions::periodic,boundary_conditions::periodic,
        boundary_conditions::periodic,boundary_conditions::periodic,
        boundary_conditions::periodic,boundary_conditions::periodic},force_type::smoothedLJ);
    ThermoStats thermoStats(container,0.05);

    std::cout << container.to_string();

    //all particles will automatically be intialized with force={0,0,0}

    //particle p1 in middle of the side that is the negative X boundary
    container.addParticle({0,5,0},{-10,0,0},1);
    //will be particles[0]


    //particle p2 in middle of the side that is the negative Y boundary
    container.addParticle({5,0,0},{0,-10,0},1);
    //will be particles[1]


    container.createPointers();

    thermoStats.initDiffusionCoefficient();

    // as the particles are exactly at the boundary, they both should move to 
    // the respective other side within one step already
    calculator.calculateX();

    std::cout << container.to_string();




    std::vector<Particle> particles = container.getInstances();

    std::cout << particles[0].toString() << "\n" << particles[1].toString() << "\n";

    ASSERT_EQ(-1,particles[0].getBoundariesCrossed(0));
    ASSERT_EQ(-1,particles[1].getBoundariesCrossed(1));

    double diffusion_coefficient_actual = thermoStats.getDiffusionCoefficient();

    double diffusion_coefficient_expected = ( std::pow(-10 * 0.05,2) +  std::pow(-10 * 0.05,2) ) / 2.0;

    ASSERT_NEAR(diffusion_coefficient_expected,diffusion_coefficient_actual,1e-10);

}




TEST(test_DiffCoeff2, diff_Coeff_twice_crossed_boundary){
    CellContainer container{10, 10, 0, 2, 2};
    CellCalculator calculator(container,0.6,2.0,1.9,
        {boundary_conditions::periodic,boundary_conditions::periodic,
        boundary_conditions::periodic,boundary_conditions::periodic,
        boundary_conditions::periodic,boundary_conditions::periodic},force_type::smoothedLJ);
    ThermoStats thermoStats(container,0.06);


    //all particles will automatically be intialized with force={0,0,0}

    //particle p1 in middle of the side that is the negative X boundary
    container.addParticle({0,5,0},{-10,0,0},1);
    //will be particles[0]



    container.createPointers();

    thermoStats.initDiffusionCoefficient();

    // as the particles are exactly at the boundary, they both should move to 
    // the respective other side within one step already
    calculator.calculateX();
    calculator.calculateX();

    //boundary should get crossed twice

    std::cout << container.to_string();

    std::vector<Particle> particles = container.getInstances();

    std::cout << particles[0].toString() << "\n";

    ASSERT_EQ(-2,particles[0].getBoundariesCrossed(0));

    double diffusion_coefficient_actual = thermoStats.getDiffusionCoefficient();

    auto& new_particles = container.getInstances();

    for(auto particle = new_particles.begin(); particle != new_particles.end(); particle++){
        for(int i = 0; i < 3; i++){
            ASSERT_EQ(particle->getBoundariesCrossed(i),0);
        }
    }

    double diffusion_coefficient_expected = ( std::pow(-10 * 0.6 * 2,2)) / 1.0;

    ASSERT_NEAR(diffusion_coefficient_expected,diffusion_coefficient_actual,1e-10);

}

TEST(test_DiffCoeff3, diff_Coeff_very_big){
    CellContainer container{10, 10, 10, 2, 2};
    CellCalculator calculator(container,0.6,2.0,1.9,
                              {boundary_conditions::periodic,boundary_conditions::periodic,
                               boundary_conditions::periodic,boundary_conditions::periodic,
                               boundary_conditions::periodic,boundary_conditions::periodic},force_type::smoothedLJ);
    ThermoStats thermoStats(container,0.06);

    auto old_positions = thermoStats.getParticlePositionsPreviousIteration();

    //all particles will automatically be intialized with force={0,0,0}

    //particle p1 in middle of the side that is the negative X boundary / positive X boundary
    container.addParticle({0,5,5},{-10,0,0},1);
    container.addParticle({9.9,5,5},{10,0,0},1);
    //will be particles[0] / particles[1]

    //particle p1 in middle of the side that is the negative Y boundary / positive Y boundary
    container.addParticle({5,0,5},{0,-10,0},1);
    container.addParticle({5,9.9,5},{0,10,0},1);
    //will be particles[2] / particles[3]

    //particle in middle of the side that is the negative Z boundary / positive Z boundary
    container.addParticle({5,5,0},{0,0,-10},1);
    container.addParticle({5,5,9.9},{0,0,10},1);
    //will be particles[4] / particles[5]


    container.createPointers();

    thermoStats.initDiffusionCoefficient();

    calculator.calculateX();
    calculator.calculateX();

    //boundary should get crossed twice

    std::cout << container.to_string();



    std::vector<Particle> particles = container.getInstances();

    std::cout << particles[0].toString() << "\n";

    ASSERT_EQ(-2,particles[0].getBoundariesCrossed(0));
    ASSERT_EQ(2,particles[1].getBoundariesCrossed(0));

    ASSERT_EQ(-2,particles[2].getBoundariesCrossed(1));
    ASSERT_EQ(2,particles[3].getBoundariesCrossed(1));

    ASSERT_EQ(-2,particles[4].getBoundariesCrossed(2));
    ASSERT_EQ(2,particles[5].getBoundariesCrossed(2));

    double diffusion_coefficient_actual = thermoStats.getDiffusionCoefficient();

    auto& new_particles = container.getInstances();

    for(auto particle = new_particles.begin(); particle != new_particles.end(); particle++){
        for(int i = 0; i < 3; i++){
            ASSERT_EQ(particle->getBoundariesCrossed(i),0);
        }
    }

    ASSERT_NE(old_positions,thermoStats.getParticlePositionsPreviousIteration());

    double diffusion_coefficient_expected = std::pow(10 * 0.6 * 2,2);

    ASSERT_NEAR(diffusion_coefficient_expected,diffusion_coefficient_actual,1e-10);

}

