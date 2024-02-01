#include <gtest/gtest.h>
#include "particleModel/storage/CellContainer.h"
#include "particleModel/updating/CellContainerIterators.h"

#include "particleModel/updating/CellCalculator.h"
#include "utils/ArrayUtils.h"



/**
 * @brief First a CellContainer an CellCalculator are created, then three particles
 *        are added such, that one particle p1 is on the lower left corner of a 2D rectangle
 *        and the other two are at the left upper (p2) and right lower corner (p3). According to
 *        the logic of the periodic boundary conditions, there should be forces between p1 and p2 
 *        and forces between p1 and p3. The forces are calculated and it is checked, weather
 *        the particles have the expected forces
 * 
 *        |p3------|
 *        |        |
 *        |        |
 *        |p1____p2|
*/
TEST(test_CellCalculation, test_periodic_forces){
    CellContainer container{10, 10, 0, 2, 2};
    CellCalculator calculator(container,0.0014,2.0,1.9,
        {boundary_conditions::periodic,boundary_conditions::periodic,
        boundary_conditions::periodic,boundary_conditions::periodic,
        boundary_conditions::periodic,boundary_conditions::periodic},force_type::LJ);


    //particle p1 in corner 
    //0.1 so the particles don't overlap
    container.addParticle({0.5,0.5,0},{0,0,0},1);
    //will be particles[0]


    //particle p2 in lower right corner
    container.addParticle({9.0,0.5,0},{0,0,0},1);
    //will be particles[1]

    //particle p3 in upper left corner
    container.addParticle({0.5,9.0,0},{0,0,0},1);
    //will be particles[2]

    //-> there should be forces between p1 and p2 as well as between p1 and p3

    container.createPointers();

    std::cout << container.to_string() << "\n";

    //should calculate forces for the appropriate particles
    calculator.calculateF();

    std::vector<Particle> particles = container.getInstances();

    for(auto& particle : particles){
        //forces should not be zero for all particles
        ASSERT_TRUE(ArrayUtils::L2Norm(particle.getF()) > 0);
    }    
    std::array<double,3> p1_force = particles[0].getF();
    std::array<double,3> p2_force = particles[1].getF();
    std::array<double,3> p3_force = particles[2].getF();

    ASSERT_NEAR(p1_force[0],-p2_force[0],0.0000001);
    // operator- does operation on double -> use standard method for comparing two doubles
    ASSERT_NEAR(p1_force[1],-p3_force[1],0.0000001);
    // operator+ does operation on double -> use standard method for comparing two doubles 

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
TEST(test_CellCalculation, test_periodic_reappearing){
    CellContainer container{10, 10, 0, 2, 2};
    CellCalculator calculator(container,0.0014,2.0,1.9,
        {boundary_conditions::periodic,boundary_conditions::periodic,
        boundary_conditions::periodic,boundary_conditions::periodic,
        boundary_conditions::periodic,boundary_conditions::periodic},force_type::LJ);


    //all particles will automatically be intialized with force={0,0,0}

    //particle p1 in middle of the side that is the negative X boundary
    container.addParticle({0,5,0},{-10,0,0},1);
    //will be particles[0]


    //particle p2 in middle of the side that is the negative Y boundary
    container.addParticle({5,0,0},{0,-10,0},1);
    //will be particles[1]


    container.createPointers();


    // as the particles are exactly at the boundary, they both should move to 
    // the respective other side within one step already
    calculator.calculateX();

    std::vector<Particle> particles = container.getInstances();


    // test if the particle p1 was moved to the opposite side of 
    // the domain (x ~ 10, wheres in the beginning x[0] = 0)
    std::array<double,3> x_1 = particles[0].getX();
    ASSERT_GE(x_1[0],9);

    // test if the particle p1 was moved to the opposite side of 
    // the domain (y ~ 10 , whereas in the beginning y = x[1] = 0)
    std::array<double,3> x_2 = particles[1].getX();
    ASSERT_GE(x_2[1],9);
    
}



/**
 * @brief Does the same as the previous case, only that the particle is in the corner
 *        or our 2D rectangle domain and should reappear on the opposite corner
 * 
 * 
 * 
*/
TEST(test_CellCalculation, test_periodic_corner){
    CellContainer container{10, 10, 0, 2, 2};
    CellCalculator calculator(container,0.0014,2.0,1.9,
        {boundary_conditions::periodic,boundary_conditions::periodic,
        boundary_conditions::periodic,boundary_conditions::periodic,
        boundary_conditions::periodic,boundary_conditions::periodic},force_type::LJ);


    //all particles will automatically be intialized with force={0,0,0}

    //particle p1 in middle of the side that is the negative X boundary
    container.addParticle({0,0,0},{-10,-10,0},1);
    //will be particles[0]


    container.createPointers();


    // as the particle is exactly at the boundary -> should move to 
    // the respective other sides within one step already
    calculator.calculateX();

    std::vector<Particle> particles = container.getInstances();

    // test if the particle p1 was moved to the opposite side of 
    // the domain (x ~ 10, y ~ 10 wheres in the beginning x[0] = 0  and y = x[1] = 0)
    std::array<double,3> x = particles[0].getX();
    ASSERT_GE(x[0],9);
    ASSERT_GE(x[1],9);
    
}
