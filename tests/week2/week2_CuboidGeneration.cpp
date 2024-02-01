#include <gtest/gtest.h>
#include "inputHandling/FileReader.h"
#include "inputHandling/generators/CuboidGeneration.h"
#include "particleModel/CellContainer.h"
#include "inputHandling/FileReaderProgramArgs.h"

/**
 * @brief Tests the CuboidGenerator. Four cuboids are specified by CuboidData structs,
 *        then they are added/ created with the CuboidGenerator. Then in the beginning
 *        it is checked if there was no random velocity added in the 3D dimension, because
 *        the CuboidGenerator should automatically recognize only two dimensional Cuboids
 *        and then create random velocities with only two dimensions. After adding three-dimensional
 *        cuboids, the CuboidGenerator should then also initalize random velocities in 
 *        three dimensions
*/
TEST(cuboidgeneration,test_cuboid_dimension){
    FileReader::CuboidData c1{{0, 0, 0}, {0, 0, 0}, 300, 5, 1, 1, 1.1225, 1, 5};
    FileReader::CuboidData c2{{0, 0, 0}, {0, 0, 0}, 1, 1, 1, 1, 1.1225, 1, 5, 0};
    FileReader::CuboidData c3{{0, 0, 0}, {0, 0, 0}, 1, 1, 2, 1, 1.1225, 1, 5, 0};
    FileReader::CuboidData c4{{0, 0, 1}, {0, 0, 0}, 1, 1, 1, 1, 1.1225, 1, 5, 0};
    std::list<FileReader::CuboidData> cuboids{c1, c2};

    //check if velocity dimension is set to 2
    CellContainer particleContainer(1000,1000,0,2.0,2.0);
    addCuboids(particleContainer, cuboids);

    for(int i = 0; i < particleContainer.size(); i++) {
        ASSERT_EQ(particleContainer.getInstances()[i].getV().at(2), 0.0);
    }

    //any further dimension tests don't make any sense for the CellContainer
}


/**
 * @brief Tests the Cuboidgenerator by first creating two cuboids   
 *        specified by c1 and c2, then tests if first and last particles position
 *        from the second cuboid got created correctly
 * 
 * 
*/
TEST(cuboidgeneration,test_cuboidgeneration){
    double h{1.1225};
    FileReader::CuboidData c1{{10, 10, 10}, {0, 0, 0}, 10, 5, 10, 1, h, 1, 5};
    FileReader::CuboidData c2{{0, 0, 0}, {0, 0, 0}, 10, 15, 20, 1, h, 1, 5};
    std::list<FileReader::CuboidData> cuboids{c1, c2};

    //check if first and last particles position from second cuboid got created correctly
    CellContainer particleContainer(100,100,0,1.0,1.0);
    addCuboids(particleContainer, cuboids);
    std::array<double, 3> x_f{0.0, 0.0, 0.0};
    std::array<double, 3> x_s{9 * h, 14 * h, 19 * h};

    bool found_first{false};
    bool found_second{false};
    auto& particles = particleContainer.getInstances();
    for(int i = 0; i < particles.size(); i++) {
        if(particles[i].getX() == x_f) {
            found_first = true;
        } else if(particles[i].getX() == x_s) {
            found_second = true;
        }
    }
    ASSERT_TRUE(found_first);
    ASSERT_TRUE(found_second);
}
