#include <gtest/gtest.h>
#include "particleModel/storage/CellContainer.h"
#include "particleModel/storage/CellContainerIterators.h"

#include "particleModel/updating/CellCalculator.h"
#include "utils/ArrayUtils.h"

/**
 * @brief
 */
TEST(test_CellCalculation,test_periodic){
    CellContainer periodic{10, 10, 10, 2, 2};
    CellContainer not_periodic{10, 100, 10, 2, 2};

    std::array<double,3> x1{5,5,5};
    std::array<double,3> x2{5,6.1225,5};
    std::array<double,3> x3{5,7.245,5};
    std::array<double,3> v{0,5,0};
    double m = 1;
    //std::array<

    periodic.addParticle(x1,v,m);
    periodic.addParticle(x2,v,m);
    periodic.addParticle(x3,v,m);
    not_periodic.addParticle(x1,v,m);
    not_periodic.addParticle(x2,v,m);
    not_periodic.addParticle(x3,v,m);
}


TEST(test_CellCalculation, test_periodic_simple){
    CellContainer container{10, 10, 0, 2, 2};
    CellCalculator calculator(container,0.0014,2.0,
        {boundary_conditions::periodic,boundary_conditions::periodic,
        boundary_conditions::periodic,boundary_conditions::periodic,
        boundary_conditions::periodic,boundary_conditions::periodic});


    //particle p1 in corner 
    //0.1 so the particles don't overlap
    container.addParticle({0.1,0.1,0},{0,0,0},1);
    //will be particles[0]


    //particle p2 in lower right corner
    container.addParticle({9.9,0,0},{0,0,0},1);
    //will be particles[1]

    //particle p3 in upper left corner
    container.addParticle({0,9.9,0},{0,0,0},1);
    //will be particles[2]

    //-> there should be forces between p1 and p2 as well as between p1 and p3

    container.createPointers();

    //should calculate forces for the appropriate particles
    calculator.calculateF();

    std::vector<Particle>& particles = container.getInstances();

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
