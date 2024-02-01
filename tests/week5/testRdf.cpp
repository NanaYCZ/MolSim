#include <gtest/gtest.h>
#include "particleModel/updating/ThermoStats.h"
#include "particleModel/storage/CellContainer.h"

#include "inputHandling/FileReaderProgramArgs.h"



/**
 * @brief first a equal-sided rectangle of particles is added to a CellContainer, where the particles
 *        are the corner points of the equal-sided rectangle. Then we know, that when looking at all
 *        distances between two particle pairs, we have 4 times the side length of 4 and two times the
 *        distance of the diagonal of the equal-sided rectangle, which have length sqrt(32). Knowing
 *        this information, the rdf is calculated for the intervals, in which we expect a count higher
 *        than zero. Then all intervals, that the rdf produces are checked for correctness.
 *
 */
TEST(test_Rdf,test_rectangle){
    CellContainer container(15,15,0,3.0,3.0);
    ThermoStats thermoStats(container,0.1,0,0);

    // basic rectangle -> 4 * (side of length 4) + 2 * (sqrt(32) (diagonals))
    container.addParticle({0,0,0},{0,0,0},3);
    container.addParticle({4,0,0},{0,0,0},4);
    container.addParticle({0,4,0},{0,0,0},4);
    container.addParticle({4,4,0},{0,0,0},7);


    std::vector<double> rdf = thermoStats.getRadialDistributionFunction(0.1);


    double expected_in_40 = 4 / ( (4.0 * M_PI / 3.0) *
                                   ( std::pow(4.1,3) - std::pow(4.0,3)) );

    double expected_in_56 = 2 / ( (4.0 * M_PI / 3.0) *
                              ( std::pow(5.7,3) - std::pow(5.6,3)) );



    std::cout << ArrayUtils::to_string(rdf) << "\n";

    for(size_t i = 0; i < rdf.size();i++){
    if(i == 40){ // 20 * 0.1(step size) = 2.0 -> bucket 2.0 - 2.1
        ASSERT_NEAR(expected_in_40,rdf[i],1e-7);
    }else if(i == 56){ // 28 * 0.1(step size) = 2.8 -> bucket 2.8 - 2.9
        ASSERT_NEAR(expected_in_56,rdf[i],1e-7);
    }else{
        ASSERT_EQ(0,rdf[i]);
    }
    }
}


/**
 * @brief first a equal-sided pyramid of particles is added to a CellContainer, where the particles
 *        are the corner points of the equal-sided pyramid. Then we know, that when looking at all
 *        distances between two particle pairs, there is only one distance that should exist, namely
 *        the side length of the pyramid, which is 2 * sqrt(2) or sqrt(8). Knowing
 *        this information, the rdf is calculated for the intervals, in which we expect a count higher
 *        than zero. Then all intervals, that the rdf produces are checked for correctness.
 */
TEST(test_Rdf,test_pyramid){
    CellContainer container(15,15,0,3.0,3.0);
    ThermoStats thermoStats(container,0.1,0,0);

    // basic equal-sided pyramid (every side has length 2 * sqrt(2) ~ 2.828 )
    container.addParticle({2,2,2},{0,0,0},3);
    container.addParticle({2,0,0},{0,0,0},4);
    container.addParticle({0,2,0},{0,0,0},4);
    container.addParticle({0,0,2},{0,0,0},7);

    std::vector<double> rdf = thermoStats.getRadialDistributionFunction(0.1);
    //-> there should be 6 entries be counted in the 28th "backet"
    //because there are 6 pairs with the distance 2.828 -> they belong in 2.8 - 2.9 bucket

    double expected = 6 / ( (4.0 * M_PI / 3.0) * 
                      ( std::pow(2.9,3) - std::pow(2.8,3)) );
    
    std::cout << ArrayUtils::to_string(rdf) << "\n";

    for(size_t i = 0; i < rdf.size();i++){
        if(i == 28){ // 28 * 0.1(step size) = 2.8 -> bucket 2.8 - 2.9
            ASSERT_NEAR(expected,rdf[i],1e-7);
        }else{
            ASSERT_EQ(0,rdf[i]);
        }
    }

}



/**
 * @brief first a cuboid of particles is added to a CellContainer, where the particles
 *        are the corner points of the cuboid. Then we know, that when looking at all
 *        distances between two particle pairs, There are exactly three different distances,
 *        the side length of the cuboid, which is 2 and this distance exists 12 times,
 *        the diagonals of the side areas of the cuboid, which is sqrt(8) ans this exists 12 times as well
 *        and the space diagonals, that go 'through' the volume of the cuboid and that have
 *        length sqrt(8 + 4). There are 4 of these space diagonals. Knowing this information,
 *        the rdf is calculated for the intervals, in which we expect a count higher
 *        than zero. Then all intervals, that the rdf produces are checked for correctness.
 */
TEST(test_Rdf,test_cuboid){
    CellContainer container(15,15,0,3.0,3.0);
    ThermoStats thermoStats(container,0.1,0,0);

    // basic cuboid: 12 * (side of length 2) + 12 * (sqrt(8) (2 diagonals for every side) ) +
    // 4 * (sqrt(8 + 4)  (4 space diagonals in total)) )
    container.addParticle({0,0,0},{0,0,0},3);
    container.addParticle({2,0,0},{0,0,0},4);
    container.addParticle({0,2,0},{0,0,0},4);
    container.addParticle({2,2,0},{0,0,0},7);

    container.addParticle({0,0,2},{0,0,0},3);
    container.addParticle({2,0,2},{0,0,0},4);
    container.addParticle({0,2,2},{0,0,0},4);
    container.addParticle({2,2,2},{0,0,0},7);

    std::vector<double> rdf = thermoStats.getRadialDistributionFunction(0.1);


    double expected_in_20 = 12 / ( (4.0 * M_PI / 3.0) *
                            ( std::pow(2.1,3) - std::pow(2.0,3)) );

    double expected_in_28 = 12 / ( (4.0 * M_PI / 3.0) *
                               ( std::pow(2.9,3) - std::pow(2.8,3)) );

    double expected_in_34 = 4 / ( (4.0 * M_PI / 3.0) *
                               ( std::pow(3.5,3) - std::pow(3.4,3)) );

    std::cout << ArrayUtils::to_string(rdf) << "\n";

    for(size_t i = 0; i < rdf.size();i++){
        if(i == 20){ // 20 * 0.1(step size) = 2.0 -> bucket 2.0 - 2.1
            ASSERT_NEAR(expected_in_20,rdf[i],1e-7);
        }else if(i == 28){ // 28 * 0.1(step size) = 2.8 -> bucket 2.8 - 2.9
            ASSERT_NEAR(expected_in_28,rdf[i],1e-7);
        }else if(i == 34){ // 34 * 0.1(step size) = 3.4 -> bucket 3.4 - 3.5
            ASSERT_NEAR(expected_in_34,rdf[i],1e-7);
        }else{
            ASSERT_EQ(0,rdf[i]);
        }
    }

}





