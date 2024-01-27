#include <gtest/gtest.h>
#include "particleModel/storage/CellContainer.h"

//each test checks if the CellContainer, when initialized, calculates the patterns and
//stores them correctly, by comparing them with an example of an expected vector of arrays.

TEST(test_precalculating_patterns,test_depth_1_2D){
    double width = 30;
    double height = 30;
    double depth = 0;
    double r_cutoff = 3.0;
    double cell_size = 3.0;
    CellContainer cellContainer{width, height, depth, r_cutoff, cell_size};
    std::vector<std::array<dim_t,3>> patterns = {{1,0,0},{1,1,0},{0,1,0},{-1,1,0}};

    for(std::array<dim_t,3> pattern : cellContainer.getPatterns()) {
            auto found = std::find(patterns.begin(), patterns.end(), pattern);
            ASSERT_NE(found, patterns.end());
    }
}

TEST(test_precalculating_patterns,test_depth_2_2D){
    double width = 30;
    double height = 30;
    double depth = 0;
    double r_cutoff = 3.0;
    double cell_size = 1.5;
    CellContainer cellContainer{width, height, depth, r_cutoff, cell_size};
    std::vector<std::array<dim_t,3>> patterns = {{1,0,0},{1,1,0},{0,1,0},{-1,1,0},
                                                 {2,0,0},{2,1,0},{2,2,0},{1,2,0},
                                                 {0,2,0},{-1,2,0},{-2,-2,0},{-2,-1,0}};

    for(std::array<dim_t,3> pattern : cellContainer.getPatterns()) {
        auto found = std::find(patterns.begin(), patterns.end(), pattern);
        ASSERT_NE(found, patterns.end());
    }
}

TEST(test_precalculating_patterns,test_depth_1_3D){
    double width = 30;
    double height = 30;
    double depth = 30;
    double r_cutoff = 3.0;
    double cell_size = 3.0;
    CellContainer cellContainer{width, height, depth, r_cutoff, cell_size};
    std::vector<std::array<dim_t,3>> patterns = {{1,0,0},{1,1,0},{0,1,0},{-1,1,0},
                                                 {0,0,1},{1,0,1},{1,1,1},{0,1,1},
                                                 {-1,1,1},{1,0,-1},{1,1,-1},{0,1,-1},
                                                 {-1,1,-1}};

    for(std::array<dim_t,3> pattern : cellContainer.getPatterns()) {
        auto found = std::find(patterns.begin(), patterns.end(), pattern);
        ASSERT_NE(found, patterns.end());
    }
}

TEST(test_precalculating_patterns,test_depth_2_3D){
    double width = 30;
    double height = 30;
    double depth = 30;
    double r_cutoff = 3.0;
    double cell_size = 1.5;
    CellContainer cellContainer{width, height, depth, r_cutoff, cell_size};
    std::vector<std::array<dim_t,3>> patterns = {{1,0,0},{1,1,0},{0,1,0},{-1,1,0},
                                                 {2,0,0},{2,1,0},{2,2,0},{1,2,0},
                                                 {0,2,0},{-1,2,0},{-2,-2,0},{-2,-1,0},
                                                 {0,0,1},{1,0,1},{1,1,1},{0,1,1},
                                                 {-1,1,1},{1,0,-1},{1,1,-1},{0,1,-1},
                                                 {-1,1,-1}, {0,0,2},{-1,1,2},{0,1,2},
                                                 {2,0,1},{2,1,1},{2,2,1},{1,2,1},
                                                 {0,2,1},{-1,2,1},{-2,-2,1},{-2,-1,1},
                                                 {2,0,-1},{2,1,-1},{2,2,-1},{1,2,-1},
                                                 {0,2,-1},{-1,2,-1},{-2,-2,-1},{-2,-1,-1},
                                                 {2,0,2},{2,1,2},{2,2,2},{1,2,2},
                                                 {0,2,2},{-1,2,2},{-2,-2,2},{-2,-1,2},
                                                 {2,0,-2},{2,1,-2},{2,2,-2},{1,2,-2},
                                                 {0,2,-2},{-1,2,-2},{-2,-2,-2},{-2,-1,-2},
                                                 {1,0,-2},{1,1,-2},{0,1,-2},{-1,1,-2},
                                                 {1,0,2},{1,1,2}};

    for(std::array<dim_t,3> pattern : cellContainer.getPatterns()) {
        auto found = std::find(patterns.begin(), patterns.end(), pattern);
        ASSERT_NE(found, patterns.end());
    }
}