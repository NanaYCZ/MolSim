#include <gtest/gtest.h>
#include "particleModel/updating/CellContainerIterators.h"

//each test checks if the iterator provides the correct points,
//by comparing the result with an example of an expected collection of points.


/**
 * @brief test iteration over all cells of a 2D domain
*/
TEST(test_iterators,cell_iterator_2D){
    double width = 9;
    double height = 9;
    double depth = 0;
    double r_cutoff = 3.0;
    double cell_size = 3.0;
    CellContainer cellContainer{width, height, depth, r_cutoff, cell_size};
    std::vector<std::array<dim_t,3>> cells = {{1,1,1},{2,1,1},{3,1,1},
                                              {1,2,1},{2,2,1},{3,2,1},
                                              {1,3,1},{2,3,1},{3,3,1}};

    int counter = 0;
    for(auto iterator = begin_CI(); iterator != end_CI(); ++iterator) {
        ++counter;
        auto found = std::find(cells.begin(), cells.end(), iterator.position());
        ASSERT_NE(found, cells.end());
    }

    ASSERT_EQ(counter, cells.size());
}

/**
 * @brief test iteration over all cells of a 3D domain
*/
TEST(test_iterators,cell_iterator_3D){
    double width = 9;
    double height = 9;
    double depth = 9;
    double r_cutoff = 3.0;
    double cell_size = 3.0;
    CellContainer cellContainer{width, height, depth, r_cutoff, cell_size};
    std::vector<std::array<dim_t,3>> cells = {{1,1,1},{2,1,1},{3,1,1},
                                              {1,2,1},{2,2,1},{3,2,1},
                                              {1,3,1},{2,3,1},{3,3,1},
                                              {1,1,2},{2,1,2},{3,1,2},
                                              {1,2,2},{2,2,2},{3,2,2},
                                              {1,3,2},{2,3,2},{3,3,2},
                                              {1,1,3},{2,1,3},{3,1,3},
                                              {1,2,3},{2,2,3},{3,2,3},
                                              {1,3,3},{2,3,3},{3,3,3}};
    int counter = 0;
    for(auto iterator = begin_CI(); iterator != end_CI(); ++iterator) {
        ++counter;
        auto found = std::find(cells.begin(), cells.end(), iterator.position());
        ASSERT_NE(found, cells.end());
    }

    ASSERT_EQ(counter, cells.size());
}

/**
 * @brief test iteration over all starting points in a 2D domain,
 * for a pattern in one dimension(direction)
*/
TEST(test_iterators,start_point_iterator_2D_1DM){
    double width = 12;
    double height = 12;
    double depth = 0;
    double r_cutoff = 3.0;
    double cell_size = 3.0;
    CellContainer cellContainer{width, height, depth, r_cutoff, cell_size};
    std::array<dim_t,3> pattern = {2,0,0};
    std::vector<std::array<dim_t,3>> points = {{1,1,1},{1,2,1},{1,3,1},{1,4,1},
                                               {2,1,1},{2,2,1},{2,3,1},{2,4,1}};

    int counter = 0;
    for(auto iterator = begin_SI(pattern); iterator != end_SI(); ++iterator) {
        ++counter;
        auto found = std::find(points.begin(), points.end(), *iterator);
        ASSERT_NE(found, points.end());
    }

    ASSERT_EQ(counter, points.size());
}

/**
 * @brief test iteration over all starting points in a 2D domain,
 * for a pattern in two dimensions(direction)
*/
TEST(test_iterators,start_point_iterator_2D_2DM){
    double width = 12;
    double height = 12;
    double depth = 0;
    double r_cutoff = 3.0;
    double cell_size = 3.0;
    CellContainer cellContainer{width, height, depth, r_cutoff, cell_size};
    std::array<dim_t,3> pattern = {1,-2,0};
    std::vector<std::array<dim_t,3>> points = {{1,1,1},{1,2,1},{1,3,1},{1,4,1},
                                               {2,4,1},{3,4,1},{4,4,1},{2,3,1},
                                               {3,3,1},{4,3,1}};

    int counter = 0;
    for(auto iterator = begin_SI(pattern); iterator != end_SI(); ++iterator) {
        ++counter;
        auto found = std::find(points.begin(), points.end(), *iterator);
        ASSERT_NE(found, points.end());
    }

    ASSERT_EQ(counter, points.size());
}

/**
 * @brief test iteration over all starting points in a 3D domain,
 * for a pattern in one dimension(direction)
*/
TEST(test_iterators,start_point_iterator_3D_1DM){
    double width = 9;
    double height = 9;
    double depth = 9;
    double r_cutoff = 3.0;
    double cell_size = 3.0;
    CellContainer cellContainer{width, height, depth, r_cutoff, cell_size};
    std::array<dim_t,3> pattern = {1,0,0};

    std::vector<std::array<dim_t,3>> points = {{1,1,1},{1,2,1},{1,3,1},
                                               {1,1,2},{1,2,2},{1,3,2},
                                               {1,1,3},{1,2,3},{1,3,3}};

    int counter = 0;
    for(auto iterator = begin_SI(pattern); iterator != end_SI(); ++iterator) {
        ++counter;
        auto found = std::find(points.begin(), points.end(), *iterator);
        ASSERT_NE(found, points.end());
    }

    ASSERT_EQ(counter, points.size());
}

/**
 * @brief test iteration over all starting points in a 3D domain,
 * for a pattern in two dimensions(direction)
*/
TEST(test_iterators,start_point_iterator_3D_2DM){
    double width = 9;
    double height = 9;
    double depth = 9;
    double r_cutoff = 3.0;
    double cell_size = 3.0;
    CellContainer cellContainer{width, height, depth, r_cutoff, cell_size};
    std::array<dim_t,3> pattern = {1,1,0};

    std::vector<std::array<dim_t,3>> points = {{1,1,1},{1,2,1},{1,3,1},
                                               {2,1,1},{3,1,1},{1,1,2},
                                               {1,2,2},{1,3,2},{2,1,2},
                                               {3,1,2},{1,1,3},{1,2,3},
                                               {1,3,3},{2,1,3},{3,1,3}};

    int counter = 0;
    for(auto iterator = begin_SI(pattern); iterator != end_SI(); ++iterator) {
        ++counter;
        auto found = std::find(points.begin(), points.end(), *iterator);
        ASSERT_NE(found, points.end())<<"("<<(*iterator)[0]<<","<<(*iterator)[1]<<","<<(*iterator)[2]<<") not found";
    }

    ASSERT_EQ(counter, points.size());
}

/**
 * @brief test iteration over all starting points in a 3D domain,
 * for a pattern in three dimensions(direction)
*/
TEST(test_iterators,start_point_iterator_3D_3DM){
    double width = 9;
    double height = 9;
    double depth = 9;
    double r_cutoff = 3.0;
    double cell_size = 3.0;
    CellContainer cellContainer{width, height, depth, r_cutoff, cell_size};
    std::array<dim_t,3> pattern = {1,1,-2};

    std::vector<std::array<dim_t,3>> points = {{1,1,1},{1,2,1},{1,3,1},
                                               {2,1,1},{3,1,1},{1,1,2},
                                               {1,2,2},{1,3,2},{2,1,2},
                                               {3,1,2},{1,1,3},{1,2,3},
                                               {1,3,3},{2,1,3},{3,1,3},
                                               {2,2,3},{3,2,3},{2,3,3},
                                               {3,3,3},{2,2,2},{3,2,2},
                                               {2,3,2},{3,3,2}};

    int counter = 0;
    for(auto iterator = begin_SI(pattern); iterator != end_SI(); ++iterator) {
        ++counter;
        auto found = std::find(points.begin(), points.end(), *iterator);
        ASSERT_NE(found, points.end())<<"("<<(*iterator)[0]<<","<<(*iterator)[1]<<","<<(*iterator)[2]<<") not found";
    }

    ASSERT_EQ(counter, points.size());
}

/**
 * @brief test iteration over all points outside the 2D domain,
 * which are the last step of a path for a given pattern in one
 * dimension(direction)
*/
TEST(test_iterators,periodic_iterator_2D_1DM){
    double width = 12;
    double height = 12;
    double depth = 0;
    double r_cutoff = 3.0;
    double cell_size = 3.0;
    CellContainer cellContainer{width, height, depth, r_cutoff, cell_size};
    std::array<dim_t,3> pattern = {-1,0,0};

    std::vector<std::array<dim_t,3>> points = {{0,1,1},{0,2,1},{0,3,1},{0,4,1}};

    int counter = 0;
    for(auto iterator = begin_SI(pattern); iterator != end_SI(); ++iterator) {
        ++counter;
        auto found = std::find(points.begin(), points.end(), iterator.outside());
        ASSERT_NE(found, points.end());
    }

    ASSERT_EQ(counter, points.size());
}

/**
 * @brief test iteration over all points outside the 2D domain,
 * which are the last step of a path for a given pattern in two
 * dimensions(direction)
*/
TEST(test_iterators,periodic_iterator_2D_2DM){
    double width = 12;
    double height = 12;
    double depth = 0;
    double r_cutoff = 3.0;
    double cell_size = 3.0;
    CellContainer cellContainer{width, height, depth, r_cutoff, cell_size};
    std::array<dim_t,3> pattern = {-1,1,0};

    std::vector<std::array<dim_t,3>> points = {{0,2,1},{0,3,1},{0,4,1},{0,5,1},
                                               {1,5,1},{2,5,1},{3,5,1}};

    int counter = 0;
    for(auto iterator = begin_SI(pattern); iterator != end_SI(); ++iterator) {
        ++counter;
        auto found = std::find(points.begin(), points.end(), iterator.outside());
        ASSERT_NE(found, points.end());
    }

    ASSERT_EQ(counter, points.size());
}

/**
 * @brief test iteration over all points outside the 3D domain,
 * which are the last step of a path for a given pattern in one
 * dimension(direction)
*/
TEST(test_iterators,periodic_iterator_3D_1DM){
    double width = 12;
    double height = 12;
    double depth = 12;
    double r_cutoff = 3.0;
    double cell_size = 3.0;
    CellContainer cellContainer{width, height, depth, r_cutoff, cell_size};
    std::array<dim_t,3> pattern = {-1,0,0};

    std::vector<std::array<dim_t,3>> points = {{0,1,1},{0,2,1},{0,3,1},{0,4,1},
                                               {0,1,2},{0,2,2},{0,3,2},{0,4,2},
                                               {0,1,3},{0,2,3},{0,3,3},{0,4,3},
                                               {0,1,4},{0,2,4},{0,3,4},{0,4,4}};

    int counter = 0;
    for(auto iterator = begin_SI(pattern); iterator != end_SI(); ++iterator) {
        ++counter;
        auto found = std::find(points.begin(), points.end(), iterator.outside());
        ASSERT_NE(found, points.end());
    }

    ASSERT_EQ(counter, points.size());
}

/**
 * @brief test iteration over all points outside the 3D domain,
 * which are the last step of a path for a given pattern in two
 * dimensions(direction)
*/
TEST(test_iterators,periodic_iterator_3D_2DM){
    double width = 12;
    double height = 12;
    double depth = 12;
    double r_cutoff = 3.0;
    double cell_size = 3.0;
    CellContainer cellContainer{width, height, depth, r_cutoff, cell_size};
    std::array<dim_t,3> pattern = {-1,1,0};

    std::vector<std::array<dim_t,3>> points = {{0,2,1},{0,3,1},{0,4,1},{0,5,1},
                                               {1,5,1},{2,5,1},{3,5,1},{3,5,4},
                                               {0,2,2},{0,3,2},{0,4,2},{0,5,2},
                                               {1,5,2},{2,5,2},{3,5,2},{2,5,4},
                                               {0,2,3},{0,3,3},{0,4,3},{0,5,3},
                                               {1,5,3},{2,5,3},{3,5,3},{1,5,4},
                                               {0,2,4},{0,3,4},{0,4,4},{0,5,4}};

    int counter = 0;
    for(auto iterator = begin_SI(pattern); iterator != end_SI(); ++iterator) {
        ++counter;
        auto found = std::find(points.begin(), points.end(), iterator.outside());
        ASSERT_NE(found, points.end());
    }

    ASSERT_EQ(counter, points.size());
}

/**
 * @brief test iteration over all points outside the 3D domain,
 * which are the last step of a path for a given pattern in three
 * dimensions(direction)
*/
TEST(test_iterators,periodic_iterator_3D_3DM){
    double width = 12;
    double height = 12;
    double depth = 12;
    double r_cutoff = 3.0;
    double cell_size = 3.0;
    CellContainer cellContainer{width, height, depth, r_cutoff, cell_size};
    std::array<dim_t,3> pattern = {-1,1,1};

    std::vector<std::array<dim_t,3>> points = {{0,2,2},{0,3,2},{0,4,2},{3,4,5},
                                               {0,5,2},{1,5,2},{2,5,2},{3,5,2},
                                               {0,2,3},{0,3,3},{0,4,3},{2,4,5},
                                               {0,5,3},{1,5,3},{2,5,3},{3,5,3},
                                               {0,2,4},{0,3,4},{0,4,4},{3,3,5},
                                               {0,5,4},{1,5,4},{2,5,4},{3,5,4},
                                               {0,2,5},{0,3,5},{0,4,5},{2,3,5},
                                               {0,5,5},{1,5,5},{2,5,5},{3,5,5},
                                               {1,2,5},{2,2,5},{3,2,5},{1,3,5},
                                               {1,4,5}};

    int counter = 0;
    for(auto iterator = begin_SI(pattern); iterator != end_SI(); ++iterator) {
        ++counter;
        auto found = std::find(points.begin(), points.end(), iterator.outside());
        ASSERT_NE(found, points.end())<<"("<<iterator.outside()[0]<<","<<iterator.outside()[1]<<","<<iterator.outside()[2]<<") not found";
    }

    ASSERT_EQ(counter, points.size());
}