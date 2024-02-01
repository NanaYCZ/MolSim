#pragma once

#include "utils/ForceCalculations.h"
#include "particleModel/storage/CellContainer.h"
#include <optional>

enum class boundary_conditions{
    outflow,
    reflective,
    periodic
};

enum class concurrency_strategy{
    serial,
    first_method,
    second_method
};

enum class force_type{
    gravitational,
    LJ,
    smoothedLJ
};

extern double min_distance;

extern int chunk_size;

extern std::vector<std::vector<double>> sigma_mixed;

extern std::vector<std::vector<double>> epsilon_mixed;

/**
 * @brief a list of tuples that contain information to change a particles location
 *
 * a tuple contains the pointer of the particle, the cell the pointer should be moved to.
 */
typedef std::vector<std::tuple<Particle*, std::array<dim_t,3>>> instructions;

/**
 * @class CellCalculator
 * @brief offers important functions for particle manipulation.
 *
 * This class simulates the interaction between particles inside the cell storage structure.
 * It offers the functionality to calculate the forces, velocities and position of particles.
 */
class CellCalculator {

public:
    CellCalculator(CellContainer &cellContainer, double delta_t, double cutoff, double r_l_,
                   std::array<boundary_conditions,6> boundaries_cond, force_type forceType = force_type::LJ,
                   double gravity_factor = 0, concurrency_strategy strategy = concurrency_strategy::serial);

    /**
     * @brief calculate the forces acting between the particles in two cells each, along given paths
     *
     * in order to consider all the particles in cutoff distance, we have to look at a certain amount
     * of cells around a cell depending on the cell size. That amount is determent with the
     * "comparing_depth" attribute, which is the amount of layers we are comparing around each cell.
     * Now instead of iterating through all of the cells within the cutoff radius for each cell directly,
     * we are covering them by iterating over so called "paths" created with a starting position and a
     * shifting pattern that shows which cells should be calculated next.
     * This way we are calculating all the forces within the cutoff radius for each cell, through divide
     * and conquer with "paths" that are non overlapping for each pattern, which allows parallelization
     * without using locks, with the hope of improving performance.
     *
     * OpenMP: since all the "paths" are discrete for each "pattern", we can safely apply parallelization
     * for all starting points, which is dynamic because of varying "path" lengths
     *
     * we delegated the periodic force calculation, because the mirroring resulted in overlapping of the "paths"
     */
    void calculateInterCellF();

    /**
     * @brief calculate the forces acting between particles through the domain border
     *
     * delegates the periodic force calculation, which was previously in "calculateLinkedCellF()"
     * to avoid race conditions, when using OpenMP "parallel for"
     *
     * it iterates over all last indexes, a "path" would reach, mirrors them based on the periodic boundaries
     * and calculates the forces. Just like the inter-cell calculation did, when reaching the end of a path
     *
     * the last indexes are determined through the starting point iterator, since every starting point has
     * a corresponding point outside the domain based on the current pattern/direction. More details are
     * provided in the StartingPointIterator documentation (f.e. outside())
     *
     * OpenMP: since all the points outside the domain are mirrored on the starting points of the "pattern",
     * which are discrete and each of them corresponds to a discrete "path", the cell combinations through the
     * periodic boundaries are discrete for each pattern, allowing safe parallelization
     */
    void calculatePeriodicF();




    /**
     * @brief calculates the new positions for all particles
     * and checks if their location in the cell structure needs
     * to be updated, calling for the boundary conditions in that
     * case
     *
     * OpenMP: since the cells don't share any particles, we can safely apply "parallel for", while taking
     * care of race conditions through shared access to "cell_updates"
    */
    void calculateX();


    /**
     * @brief calculates the velocities for all particles
     *
     * OpenMP: since the cells don't share any particles, we can safely apply "parallel for"
    */
    void calculateV();


    /**
     * @brief calculates the forces for all particles by using the LinkedCellAlgorithm
     * which is divided into 3 parts: 1.inter-cell calculation 2.periodic force calculation
     * 3.force calculation with each cell
     *
     * OpenMP: since the cells don't share any particles, we can safely apply "parallel for"
     * for the force calculation within each cell
    */
    void calculateF();


    /**
     * @brief shifts the forces of all particles for correct overwriting
     *
     * OpenMP: since the cells don't share any particles, we can safely apply "parallel for"
    */
    void shiftF();

    auto& getParticles(){
        return particles;
    }

    auto getDomain_Max(){
        return domain_max_dim;
    }

    auto getDomainBounds(){
        return domain_bounds;
    }

    ForceCalculation force;
    concurrency_strategy parallelization;


private:
    CellContainer &cellContainer;
    const double gravity_factor;
    const double delta_t;
    const double cutoff;
    const double r_l;
    std::array<double,3> domain_bounds;
    std::array<dim_t, 3> domain_max_dim;

    //{positive_z,negative_z,positive_x,negative_x,positive_y,negative_y}
    std::array<boundary_conditions,6> boundaries;

    std::vector<std::vector<std::vector<std::vector<Particle*>>>> &particles;

    /**
     * @brief helper method to change the location of particles within the cell structure
     *
     * when the positions of the particles get changed a list of instructions is being created
     * that summarizes all the changes of particles between cells that have to be made. That
     * list is being processed in this method.
     *
     * @param cell_updates list of instructions to change the location of particles
     */
    void updateCells(instructions& cell_updates);

    /**
     * @brief helper method to apply boundary conditions on a particle which was detected to
     * leave it's current cell
     *
     * This method takes care of:  boundary_conditions::reflective
     *                             boundary_conditions::periodic
     *                             boundary_conditions::outflow
     * by taking according actions, if the new position of a particle would be outside
     * the domain bounds
     *
     * @param particle_ptr pointer to detected particle
     * @param new_cell_position the index of the cell the particle "tries" to move into
     * @param cell_updates for adding the resulting instruction where to move the particle into
     */
    void applyBoundaries(Particle* particle_ptr, std::array<dim_t, 3>& new_cell_position, instructions& cell_updates);


    /**
     * @brief helper method to calculate the forces not covered in "calculateInterCellF()"
     *
     * since the "calculateInterCellF()" methods only calculates the forces between cells, we are
     * now calculating the forces between the particles within a cell.
     *
     * @param current_cell the cell to calculate the force within
     */
    void calculateFWithin(std::vector<Particle*> *current_cell);

    /**
     * @brief helper method for "calculatePeriodicF()" to mirror the cell position on the other side
     *
     * used to implement the periodic boundaries, writes the changed distance of the
     * particles through the mirroring in the offset for further calculations
     *
     * @param position index outside the domain to mirror
     * @param offset from shifting the positions
     *
     * @return bool indicating if the mirroring was successful with the current boundary conditions
     */
    inline bool mirror(std::array<dim_t,3> &position, std::array<double,3> &offset);
};