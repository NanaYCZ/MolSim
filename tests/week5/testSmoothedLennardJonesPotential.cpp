#include <gtest/gtest.h>
#include "particleModel/updating/CellCalculator.h"
#include "utils/ForceCalculations.h"



/**
 * @brief Test the forceCalculation with hand-calculated values
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

    //values that were calculated with python comparison
    std::array<double,3> expected_force1,expected_force2, expected_force3;
    expected_force1 = {-14.36784445 , -1.43678445 , -2.87356889};
    expected_force2 = {-660.542705  , 1651.35676251 , 1651.35676251};
    expected_force3 = {-1149.09513034, -6894.57078205, -3332.37587799};

    std::cout << "res1: " << ArrayUtils::to_string(result1) << "\n"
              << "res2: " <<  ArrayUtils::to_string(result2) << "\n"
              << "res3: " <<  ArrayUtils::to_string(result3) << "\n";


    for(int i = 0; i < 3; i++ ){
        ASSERT_NEAR(result1[i],expected_force1[i],1e-5);
    }

    for(int i = 0; i < 3; i++ ){
        ASSERT_NEAR(result2[i],expected_force2[i],1e-5);
    }

    for(int i = 0; i < 3; i++ ){
        ASSERT_NEAR(result3[i],expected_force3[i],1e-5);
    }

}

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

}


