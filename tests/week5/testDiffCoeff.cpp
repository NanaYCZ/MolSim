#include <gtest/gtest.h>
#include "particleModel/storage/CellContainer.h"
#include "particleModel/storage/CellContainerIterators.h"

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
        boundary_conditions::periodic,boundary_conditions::periodic},"smoothedLJ");
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

    std::list<Particle>& instances = container_corner.getInstances();

    std::vector<Particle> particles(instances.begin(),instances.end());

    std::cout << particles[0].toString() <<  "\n";

    std::array<int,3> crossed_boundaries = particles[0].getBoundariesCrossed();
    ASSERT_EQ(1,crossed_boundaries[0]);
    ASSERT_EQ(1,crossed_boundaries[1]);

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
        boundary_conditions::periodic,boundary_conditions::periodic},"smoothedLJ");
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


    std::list<Particle>& instances = container.getInstances();

    std::vector<Particle> particles(instances.begin(),instances.end());

    std::cout << particles[0].toString() << "\n" << particles[1].toString() << "\n";

    ASSERT_EQ(1,particles[0].getBoundariesCrossed(0));
    ASSERT_EQ(1,particles[1].getBoundariesCrossed(1));

    double diffusion_coefficient_actual = thermoStats.getDiffusionCoefficient();

    double diffusion_coefficient_expected = ( std::pow(-10 * 0.05,2) +  std::pow(-10 * 0.05,2) ) / 2.0;

    ASSERT_NEAR(diffusion_coefficient_expected,diffusion_coefficient_actual,1e-10);

}




TEST(test_DiffCoeff2, diff_Coeff_twice_crossed_boundary){
    CellContainer container{10, 10, 0, 2, 2};
    CellCalculator calculator(container,0.6,2.0,1.9,
        {boundary_conditions::periodic,boundary_conditions::periodic,
        boundary_conditions::periodic,boundary_conditions::periodic,
        boundary_conditions::periodic,boundary_conditions::periodic},"smoothedLJ");
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


    std::list<Particle>& instances = container.getInstances();

    std::vector<Particle> particles(instances.begin(),instances.end());

    std::cout << particles[0].toString() << "\n";

    ASSERT_EQ(2,particles[0].getBoundariesCrossed(0));

    double diffusion_coefficient_actual = thermoStats.getDiffusionCoefficient();

    double diffusion_coefficient_expected = ( std::pow(-10 * 0.6 * 2,2)) / 1.0;

    ASSERT_NEAR(diffusion_coefficient_expected,diffusion_coefficient_actual,1e-10);

}

