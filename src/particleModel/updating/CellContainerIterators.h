#pragma once

#include "particleModel/storage/CellContainer.h"

typedef int dim_t;


/**
 * @class CellIterator
 *
 * @brief OpenMP "parallel for" compatible iterator for all domain cells
 */
class CellIterator {
public:
    /**
     * @brief constructor for any starting position ( begin: {1,1,1}, end: {-1,-1,-1})
     *
     * calculates the total amount of iterations for parallel distribution
     *
     * @param x index of the domain
     * @param y index of the domain
     * @param z index of the domain
     */
    CellIterator(dim_t x, dim_t y, dim_t z);

    /**
     * @brief required overload for OpenMP
     *
     * @param amount of iterations
     *
     * @return CellIterator with the new state
     */
    CellIterator &operator+=(int amount);

    /**
     * @brief sets next cell index
     *
     * @return CellIterator with the new state
     */
    CellIterator &operator++();

    /**
     * @brief extra method to get the indexes of a cell instead of a pointer to it
     *
     * @return array of indexes of the current cell
     */
    std::array<dim_t,3> position();

    /**
     * @brief dereferences the iterator to a pointer of the current cell
     *
     * @return pointer to the current cell
     */
    std::vector<Particle*>& operator*() const;

    /**
     * @brief compares iterators based on current indexes
     *
     * overload used for checking if the iterator is done,
     * primarily by sequential for loops, ignored by OpenMP scheduling
     *
     * @param other CellIterator to compare with, usually a CellIterator with end state
     *
     * @return bool based on the inequality
     */
    bool operator!=(CellIterator other) const;

    dim_t x = 1;
    dim_t y = 1;
    dim_t z = 1;

    int remaining;

    /**
     * @brief increments the cell index to the next one,
     * sets to end state {-1,-1,-1} when done
     */
    void next_index();
};

/**
 * @brief overload required by OpenMP scheduling
 *
 * gives the amount of iterations left for CellIterator a,
 * used by OpenMP "parallel for" to distribute iterations instead
 * of using the != overload
 *
 * @param a CellIterator to get the amount of iterations left from
 * @param b CellIterator, required for overload but not used
 *
 * @return amount of remaining iterations of CellIterator a
 */
int operator-(CellIterator a, CellIterator b);

/**
 * @brief provides a CellIterator in the begin state {1,1,1}
 *
 * @return instance of CellIterator
 */
CellIterator begin_CellIterator();

/**
 * @brief provides a CellIterator in the end state {-1,-1,-1}
 *
 * @return instance of CellIterator
 */
CellIterator end_CellIterator();

enum state_StartIterator {
    x_axis = 0,
    y_axis = 1,
    z_axis = 2,
    reset = -1,
    finished = -2
};


/**
 * @class StartPointIterator
 *
 * @brief OpenMP "parallel for" compatible iterator for all starting points of a "pattern"
 *
 * rebuild of legacy method SetNextPath(...)
 *
 * keeps track of a mapping, used to determine a point outside the domain, which is the result
 * of a "path" reaching the end (necessary for periodic boundaries)
 *
 * reminder: a "pattern" is the direction which will be traversed through the domain in the
 * linked cell algorithm. In order to cover all the domain cells with these "paths", certain
 * starting points have to be selected.
 * Each "pattern" has a individual set of starting points which are provided by this iterator.
 *
 * And all "paths" combined result in all the necessary combinations of cells for the linked
 * cell algorithm.
 */
class StartPointIterator {
public:
    /**
     * @brief constructor of the first iteration for a given "pattern"
     *
     * calculates the total amount of iterations for parallel distribution
     *
     * @param pattern determines the direction for which the starting points
     * will be provided
     */
    StartPointIterator(std::array<dim_t, 3> pattern);

    /**
     * @brief constructor for the end state for any "pattern"
     */
    StartPointIterator();

    /**
     * @brief required overload for OpenMP
     *
     * @param amount of iterations
     *
     * @return StartPointIterator with the new state
     */
    StartPointIterator &operator+=(int p);

    /**
      * @brief sets next cell index
      *
      * @return StartPointIterator with the new state
      */
    StartPointIterator operator++();

    /**
     * @brief compares iterators based on current indexes
     *
     * overload used for checking if the iterator is done,
     * primarily by sequential for loops, ignored by OpenMP scheduling
     *
     * @param other StartPointIterator to compare with, usually a
     * StartPointIterator with end state
     *
     * @return bool based on the inequality
     */
    bool operator!=(StartPointIterator other);

    /**
     * @brief allows iteration over end points for periodic boundaries
     *
     * for every starting point exists a point outside the domain, which is
     * the tail of a path.
     *
     * this method applies the mapping on the current starting point, which
     * provides a point outside the domain, which is the result of a "path"
     * exiting the domain by one step - no fancy math proof, just empirical :)
     *
     * @return current starting point mapped into an outside point
     */
    std::array<dim_t, 3> outside();

    /**
     * @brief dereference overload
     *
     * @return array of indexes for the current starting point
     */
    std::array<dim_t,3> operator*();

    std::array<dim_t, 3> min{1,1,1};
    std::array<dim_t, 3> max{CellContainer::domain_max_dim};
    //initial pattern which gets decremented to keep track of the remaining layers/planes
    std::array<dim_t, 3> progress;
    std::array<dim_t, 3> current{};
    std::array<dim_t, 3> mapping{};

    state_StartIterator plane_axis = reset;

    int remaining;

    /**
     * @brief increments indexes to the next starting point
     *
     * imagine the starting points are distributed as walls in a cube(domain),
     * each wall may consist multiple layers/planes and there can be in total
     * 3 walls, one for each axis.
     *
     * this method increments to the next point or the next plane + next point,
     * using the following helper methods: next_point_on_plane() and next_plane_corner()
     *
     * this way "current" is set to the next starting point within a previously mentioned plane
     */
    void next_index();

    /**
     * @brief helper method for next_index()
     *
     * iterates over every point of a plane
     *
     * increments "current" to the next starting point within the current plane a * b
     *
     * @param a value 0-2, representing one axis of the current plane
     * @param b value 0-2, representing the other axis of the current plane
     * (example: (0,1) is equal to the plane created by the x and y axis)
     */
    void next_point_on_plane(short a, short b);

    /**
     * @brief helper method for next_index()
     *
     * iterates over every plane corner, without overlaps
     *
     * sets "current" to the first position of plane for
     * next_point_on_plane(...) to iterate over and uses
     * progress attribute to keep track of the planes.
     *
     * updates mapping with the current offsets to turn
     * a starting point into a point outside of the domain
     * for periodic boundaries
     */
    void next_plane_corner();

    /**
     * @brief helper method for next_plane_corner()
     *
     * updates progress attribute, sets necessary index of "current" according
     * to plane_axis to make it a corner of the current plane and updates the
     * border for remaining planes to avoid overlap
     */
    void set_axis();
};

/**
 * @brief overload required by OpenMP scheduling
 *
 * gives the amount of iterations left for StartPointIterator b,
 * used by OpenMP "parallel for" to distribute iterations instead
 * of using the != overload
 *
 * @param a StartPointIterator to get the amount of iterations left from
 * @param b StartPointIterator, required for overload but not used
 *
 * @return amount of remaining iterations of StartPointIterator b
 */
int operator-(StartPointIterator a, StartPointIterator b);

/**
 * @brief provides a StartPointIterator with the first iteration of a starting point
 * for a given pattern
 *
 * @param pattern to determine the starting points for
 *
 * @return instance of StartPointIterator
 */
StartPointIterator begin_StartIterator(std::array<dim_t,3> pattern);

/**
 * @brief provides a StartPointIterator in the end state (progress {0,0,0}, plane_axis(finished))
 *
 * @return instance of StartPointIterator
 */
StartPointIterator end_StartIterator();
