#include <gtest/gtest.h>
#include "particleModel/updating/CellCalculator.h"

/**
 * @brief helper method to compare two particles through their distance
*/
bool approx_equal(const Particle& p1, const Particle& p2) {
    for(int i = 0; i < 3 ; i++){
        //allow for error
        if (std::fabs(p1.getX()[i] - p2.getX()[i]) > 1e-2) {
            return false;
        }
    }
    return true;
}

/**
 * @brief both tests run the same simulation setup of 6 particles in serial
 * and parallel mode and compare the results.
 * Here both simulations are the same if each particle in one simulation has
 * a particle on it's position in the other simulation.
*/
TEST(test_first_parallel,compare_2D) {
    std::array<boundary_conditions,6> bc{boundary_conditions::periodic,
                                         boundary_conditions::periodic,
                                         boundary_conditions::periodic,
                                         boundary_conditions::periodic,
                                         boundary_conditions::periodic,
                                         boundary_conditions::periodic};

    //run serial
    CellContainer cellContainer_serial{30, 30, 0, 3.0, 3.0};
    cellContainer_serial.addParticle({1,1,0},{-1,0,0},1);
    cellContainer_serial.addParticle({1,2,0},{-1,0,0},1);
    cellContainer_serial.addParticle({3,1,0},{0,-1,0},1);
    cellContainer_serial.addParticle({29,1,0},{1,0,0},1);
    cellContainer_serial.addParticle({29,2,0},{1,0,0},1);
    cellContainer_serial.addParticle({3,29,0},{0,1,0},1);
    cellContainer_serial.createPointers();

    CellCalculator cellCalculator_serial(cellContainer_serial,0.0005, 3.0, 1.9, bc,"LJ",0,concurrency_strategy::serial);
    ASSERT_EQ(cellCalculator_serial.parallelization, concurrency_strategy::serial);

    cellCalculator_serial.calculateF();
    cellCalculator_serial.shiftF();

    for (int i = 0; i < 100; ++i) {
        cellCalculator_serial.calculateX();
        cellCalculator_serial.calculateF();
        cellCalculator_serial.calculateV();
        cellCalculator_serial.shiftF();
    }

    std::vector<Particle> result_serial = cellContainer_serial.getInstances();

    //run parallel
    CellContainer cellContainer_parallel{30, 30, 0, 3.0, 3.0};
    cellContainer_parallel.addParticle({1,1,0},{-1,0,0},1);
    cellContainer_parallel.addParticle({1,2,0},{-1,0,0},1);
    cellContainer_parallel.addParticle({3,1,0},{0,-1,0},1);
    cellContainer_parallel.addParticle({29,1,0},{1,0,0},1);
    cellContainer_parallel.addParticle({29,2,0},{1,0,0},1);
    cellContainer_parallel.addParticle({3,29,0},{0,1,0},1);
    cellContainer_parallel.createPointers();

    CellCalculator cellCalculator_parallel(cellContainer_parallel,0.0005, 3.0, 1.9, bc,"LJ",0,concurrency_strategy::first_method);
    ASSERT_EQ(cellCalculator_parallel.parallelization, concurrency_strategy::first_method);

    cellCalculator_parallel.calculateF();
    cellCalculator_parallel.shiftF();

    for (int i = 0; i < 100; ++i) {
        cellCalculator_parallel.calculateX();
        cellCalculator_parallel.calculateF();
        cellCalculator_parallel.calculateV();
        cellCalculator_parallel.shiftF();
    }

    std::vector<Particle> result_parallel = cellContainer_parallel.getInstances();

    //compare results
    for(Particle& p1 : result_serial) {
        auto found = std::find_if(result_parallel.begin(), result_parallel.end(), [p1](Particle& p2) {return approx_equal(p1,p2);});
        ASSERT_NE(found, result_parallel.end());
    }
    ASSERT_EQ(result_parallel.size(), result_serial.size());
}

TEST(test_first_parallel,compare_3D){
    std::array<boundary_conditions,6> bc{boundary_conditions::periodic,
                                         boundary_conditions::periodic,
                                         boundary_conditions::periodic,
                                         boundary_conditions::periodic,
                                         boundary_conditions::periodic,
                                         boundary_conditions::periodic};

    //run serial
    CellContainer cellContainer_serial{30, 30, 30, 3.0, 3.0};
    cellContainer_serial.addParticle({1,1,0},{-1,0,0},1);
    cellContainer_serial.addParticle({1,2,0},{-1,0,0},1);
    cellContainer_serial.addParticle({3,1,0},{0,-1,0},1);
    cellContainer_serial.addParticle({29,1,0},{1,0,0},1);
    cellContainer_serial.addParticle({29,2,0},{1,0,0},1);
    cellContainer_serial.addParticle({3,29,0},{0,1,0},1);
    cellContainer_serial.createPointers();

    CellCalculator cellCalculator_serial(cellContainer_serial,0.0005, 3.0, 1.9, bc,"LJ",0,concurrency_strategy::serial);
    ASSERT_EQ(cellCalculator_serial.parallelization, concurrency_strategy::serial);

    cellCalculator_serial.calculateF();
    cellCalculator_serial.shiftF();

    for (int i = 0; i < 100; ++i) {
        cellCalculator_serial.calculateX();
        cellCalculator_serial.calculateF();
        cellCalculator_serial.calculateV();
        cellCalculator_serial.shiftF();
    }

    std::vector<Particle> result_serial = cellContainer_serial.getInstances();

    //run parallel
    CellContainer cellContainer_parallel{30, 30, 30, 3.0, 3.0};
    cellContainer_parallel.addParticle({1,1,0},{-1,0,0},1);
    cellContainer_parallel.addParticle({1,2,0},{-1,0,0},1);
    cellContainer_parallel.addParticle({3,1,0},{0,-1,0},1);
    cellContainer_parallel.addParticle({29,1,0},{1,0,0},1);
    cellContainer_parallel.addParticle({29,2,0},{1,0,0},1);
    cellContainer_parallel.addParticle({3,29,0},{0,1,0},1);
    cellContainer_parallel.createPointers();

    CellCalculator cellCalculator_parallel(cellContainer_parallel,0.0005, 3.0, 1.9, bc,"LJ",0,concurrency_strategy::first_method);
    ASSERT_EQ(cellCalculator_parallel.parallelization, concurrency_strategy::first_method);

    cellCalculator_parallel.calculateF();
    cellCalculator_parallel.shiftF();

    for (int i = 0; i < 100; ++i) {
        cellCalculator_parallel.calculateX();
        cellCalculator_parallel.calculateF();
        cellCalculator_parallel.calculateV();
        cellCalculator_parallel.shiftF();
    }

    std::vector<Particle> result_parallel = cellContainer_parallel.getInstances();

    //compare results
    for(Particle& p1 : result_serial) {
        auto found = std::find_if(result_parallel.begin(), result_parallel.end(), [p1](Particle& p2) {return approx_equal(p1,p2);});
        ASSERT_NE(found, result_parallel.end());
    }
    ASSERT_EQ(result_parallel.size(), result_serial.size());
}