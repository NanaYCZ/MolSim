#include <gtest/gtest.h>
#include "particleModel/storage/CellContainer.h"
#include "particleModel/updating/CellContainerIterators.h"

#include "particleModel/updating/CellCalculator.h"
#include "particleModel/updating/ThermoStats.h"
#include "utils/ArrayUtils.h"





/**
 * @brief First a CellContainer and a CellCalculator with all boundaries periodic are created.
 *        Then a particle is added in the lower left corner of a 2D domain with velocity into lower
 *        left direction. Then one step of calculating new positions is done, which will lead to
 *        the particle crossing the negative x and the negative y boundary. This is checked in the
 *        crossed_boundary tracker of the particle. We know that after that the particle is on the opposite
 *        side of the domain, although it only made a step of "velocity * step_size" into (negative) y dir and
 *        into (negative) x dir. Therefore we would expect, that the diffusion coefficient is
 *        (10 * 0.001)^2  + (10 * 0.001)^2 / 1 (#particles).
 *        Note: for every particle is ||x_new - x_old||^2 calculated, but this is just the scalar product
 *        of x_new - x_old, ( <x_new - x_old,x_new - x_old> ) so we calculate the scalar product of the difference
 *        and the difference is just the step in two directions (x and y) so 10 * 0.001 for x and y
 *
 *        |--------|
 *        |        |
 *        |        |
 *        |p1______|
 *       /
 *     / - velocity
 *    v
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
 * @brief First a CellContainer and a CellCalculator with all boundaries periodic are created.
 *        Then a particle is added at the middle of the negative x boundary and a particle at the
 *        middle of the negative y boundary is added. Both have velocity into the direction of the boundary
 *        at which they are. Again one step of calculating new particle positions is performed. Then again
 *        it is checked it the boundaries, that were crossed were correctly tracked in boundaries_crossed
 *        member of the particles. Again the diffusion coefficient should only be "velocity * step_size"
 *        (for both particles) although the particles moved to the other side of the domain.
 *        Therefore the diffusion coefficient is ( (10 * 0.001)^2  + (10 * 0.001)^2  )/ 2 (#particles)
 *
 *
 *        |--------|                        |---p2---|
 *        |        |                        |        |
 *   <----|p1      |  --------------->      |       p1
 *    |   |        |                        |        |
 *    |   |___p2___|                        |________|
 *    |        |
 * velocity -- |
 *             v
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





    std::vector<Particle> particles = container.getInstances();

    std::cout << particles[0].toString() << "\n" << particles[1].toString() << "\n";

    ASSERT_EQ(-1,particles[0].getBoundariesCrossed(0));
    ASSERT_EQ(-1,particles[1].getBoundariesCrossed(1));

    double diffusion_coefficient_actual = thermoStats.getDiffusionCoefficient();

    double diffusion_coefficient_expected = ( std::pow(-10 * 0.05,2) +  std::pow(-10 * 0.05,2) ) / 2.0;

    ASSERT_NEAR(diffusion_coefficient_expected,diffusion_coefficient_actual,1e-10);

}



/**
 * @brief First a CellContainer and a CellCalculator with all boundaries periodic are created.
 *        Then a particle is added at the middle of the negative x boundary.
 *        The particle again has a velocity into negative x direction, but his time the step size is quite
 *        big, namely 0.6. Then we perform  two times a step of calculating new velocities and because
 *        of the big step size, the boundary in negative X direction will be crossed twice, as
 *        (10 * 0.6 +  10 * 0.6 = 12). Therefore we check crossed_boundaries of the particle and
 *        the expected diffusions coefficient of (10 * 0.6)^2 +  (10 * 0.6)^2
 *
 *        at which they are. Again one step of calculating new particle positions is performed. Then again
 *        it is checked it the boundaries, that were crossed were correctly tracked in boundaries_crossed
 *        member of the particles. Again the diffusion coefficient should only be "velocity * step_size"
 *        (for both particles) although the particles moved to the other side of the domain.
 *        Therefore the diffusion coefficient is (10 * 0.001  + 10 * 0.001  )/ 2 (#particles)
 *
 *
 *        |--------|                        |--------|              |--------|
 *        |        |       calcX()          |        |   calcX()    |        |
 *   <----|p1      |  --------------->      |       p1  --------->  p1       |
 *    |   |        |                        |        |              |        |
 *    |   |________|                        |________|              |________|
 *    |
 * velocity
*/
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

    calculator.calculateX();
    calculator.calculateX();

    //boundary should get crossed twice


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



/**
 * @brief Lastly a big test is performed. First a CellContainer and a CellCalculator(all boundaries periodic)
 *        are created. Then for every direction x, y and z we place a particle at the respective negative boundary
 *        and one is placed at the respective positive boundary. So 6 particles in total. Each particle has
 *        velocity into the direction of the boundary at which it was placed. The step size is again big and
 *        two steps are calculated. Therefore every particle will cross its respective boundary twice. We
 *        perform the same checks as before(boundaries_crossed and diffusion coefficient) and additionally, that
 *        after calling the diffusion coefficient, the boundaries_crossed tracker is reset (should be 0) and
 *        that the particle_positions_previous_iteration member of Thermostats is updated. The member is used
 *        by the diffusion coefficient, to see what the previous positions of all the particles were and
 *        should get update with every iteration
 *
 *        no ASCII art because very difficult for 3D with arrows :(
 *
 *
*/
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

