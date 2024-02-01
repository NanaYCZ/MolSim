#include <gtest/gtest.h>
#include "particleModel/updating/CellCalculator.h"
#include "utils/ForceCalculations.h"



/**
 * @brief For various different distances between particles it is checked whether
 *        the smoothed Lennard-Jones potential is calculated correctly. The Force
 *        Calculation methods require an array of sigma and epsilon values, that
 *        will be chosen depending on the type of the particle. By default every
 *        particle has type 0. The expected values were calculated with a python
 *        script 'smoothedLJ_comparison' that is within this folder. Whenever
 *        changes are made to the force calculation the new version can be verified
 *        against these tests.
 *
*/
TEST(force_smoothed_LJ_potential,basic_big){

    std::vector<std::vector<double>> sigma_mixed = {{1}};
    std::vector<std::vector<double>> epsilon_mixed = {{1}};



    ForceCalculation force = forceSmoothedLennJonesPotentialFunction(sigma_mixed,epsilon_mixed,3.0,1.5);

    Particle p1({1.0,1,1},{0,0,0},1);
    Particle p2({2.0,1.1,1.2},{0,0,0},1);
    auto result1 = force(p1,p2,{0,0,0});

    //[5.7,4.5,1],[5.9,4.0,0.5]
    p1 = Particle({5.7,4.5,1},{0,0,0},1);
    p2 = Particle({5.9,4.0,0.5},{0,0,0},1);
    auto result2 = force(p1,p2,{0,0,0});

    //[0.0,0.1,0.01],[0.1,0.7,0.3]
    p1 = Particle({0.0,0.1,0.01},{0,0,0},1);
    p2 = Particle({0.1,0.7,0.3},{0,0,0},1);
    auto result3 = force(p1,p2,{0,0,0});

    p1 = Particle({0.34,0.231,0.13},{0,0,0},1);
    p2 = Particle({1.4,0.0,0.3},{0,0,0},1);
    auto result4 = force(p1,p2,{0,0,0});

    p1 = Particle({20.3,34.1,45.14},{0,0,0},1);
    p2 = Particle({19.7,33.9,45.6},{0,0,0},1);
    auto result5 = force(p1,p2,{0,0,0});


    //values that were calculated with python comparison
    std::array<double,3> expected_force1,expected_force2, expected_force3, expected_force4, expected_force5;
    expected_force1 = {-14.36784445 , -1.43678445 , -2.87356889};
    expected_force2 = {-660.542705  , 1651.35676251 , 1651.35676251};
    expected_force3 = {-1149.09513034, -6894.57078205, -3332.37587799};
    expected_force4 = {-1.69181557 ,0.36868811, -0.27132891};
    expected_force5 = {796.82817265 , 265.60939088, -610.90159903};

    std::cout << "res1: " << ArrayUtils::to_string(result1) << "\n"
              << "res2: " <<  ArrayUtils::to_string(result2) << "\n"
              << "res3: " <<  ArrayUtils::to_string(result3) << "\n"
              << "res4: " <<  ArrayUtils::to_string(result4) << "\n"
              << "res5: " <<  ArrayUtils::to_string(result5) << "\n";


    for(int i = 0; i < 3; i++ ){
        ASSERT_NEAR(result1[i],expected_force1[i],1e-5);
    }

    for(int i = 0; i < 3; i++ ){
        ASSERT_NEAR(result2[i],expected_force2[i],1e-5);
    }

    for(int i = 0; i < 3; i++ ){
        ASSERT_NEAR(result3[i],expected_force3[i],1e-5);
    }

    for(int i = 0; i < 3; i++ ){
        ASSERT_NEAR(result4[i],expected_force4[i],1e-5);
    }

    for(int i = 0; i < 3; i++ ){
        ASSERT_NEAR(result5[i],expected_force5[i],1e-5);
    }


}


/**
 * @brief This test case is for a edge case(that encountered to be problematic)
 *        Calculating the smoothed Lennard Jones Potential between two particles,
 *        whose distance is exactly the cutoff radius should yields zero.
 *
 *
 */
TEST(force_smoothed_LJ_potential,edge_case){

    std::vector<std::vector<double>> sigma_mixed = {{1}};
    std::vector<std::vector<double>> epsilon_mixed = {{1}};

    ForceCalculation force = forceSmoothedLennJonesPotentialFunction(sigma_mixed,epsilon_mixed,3.6,1.5);

    //these two particles exactly have 3.6 units distance between them
    //meaning they should have no force or almost no force between them
    Particle a{{4.2,1.8,1.8},{0,0,0},1};
    Particle b{{6.6,3,4.2},{0,0,0},1};
    auto result = force(a,b,{0,0,0});

    ASSERT_NEAR(result[0],0,1e-5);
    ASSERT_NEAR(result[1],0,1e-5);
    ASSERT_NEAR(result[2],0,1e-5);

}
/**
 * @brief By definition of the smoothed Lennard-Jones Potential, it should be equal to
 *        the regular Lennard-Jones-Potential if two particles are closer than r_l(
 *        which we chose to be 2.0 in this case). Therefore the force between several
 *        particle pairs (that are closer than 2.0 to each other) is calculated and
 *        it is checked whether for these particle pairs the smoothed Lennard-Jones
 *        potential is equal to the regular Lennard-Jones Potential
 */
TEST(force_smoothed_LJ_potential,compare_to_normal_LJ){

    std::vector<std::vector<double>> sigma_mixed = {{1}};
    std::vector<std::vector<double>> epsilon_mixed = {{1}};

    ForceCalculation force_smoothed = forceSmoothedLennJonesPotentialFunction(sigma_mixed,epsilon_mixed,3.0,2.0);
    ForceCalculation force_not_smoothed = forceLennJonesPotentialFunction(sigma_mixed,epsilon_mixed,3.0);


    //particles are closer than r_l=2.0 to each other
    Particle p1({1.2,3,4},{0,0,0},1);
    Particle p2({2.2,3.1,4.2},{0,0,0},1);
    auto result_smoothed = force_smoothed(p1,p2,{0,0,0});
    auto result_not_smoothed = force_not_smoothed(p1,p2,{0,0,0});

    for(int i = 0;i < 3; i++){
        ASSERT_NEAR(result_smoothed[i],result_not_smoothed[i],1e-7);
    }


    //particles are closer than r_l=2.0 to each other
    p1 = Particle({10.5,7,8},{0,0,0},1);
    p2 = Particle({10.5,6.7,8.9},{0,0,0},1);
    result_smoothed = force_smoothed(p1,p2,{0,0,0});
    result_not_smoothed = force_not_smoothed(p1,p2,{0,0,0});

    for(int i = 0;i < 3; i++){
    ASSERT_NEAR(result_smoothed[i],result_not_smoothed[i],1e-7);
    }

    //particles are closer than r_l=2.0 to each other
    p1 = Particle({0.1,0.001,0.00001},{0,0,0},1);
    p2 = Particle({1,1,0.00001},{0,0,0},1);
    result_smoothed = force_smoothed(p1,p2,{0,0,0});
    result_not_smoothed = force_not_smoothed(p1,p2,{0,0,0});

    for(int i = 0;i < 3; i++){
        ASSERT_NEAR(result_smoothed[i],result_not_smoothed[i],1e-7);
    }

}


