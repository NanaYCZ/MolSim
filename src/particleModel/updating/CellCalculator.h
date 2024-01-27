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

extern double min_distance;

extern int schedule_size;

extern std::vector<std::vector<double>> sigma_mixed;

extern std::vector<std::vector<double>> epsilon_mixed;

/**
 * @brief a list of tuples that contain information to change a particles location
 *
 * a tuple contains the index of a particle within it's cell, the current cell it's
 * located and the new cell to move it into.
 */
typedef std::vector<std::tuple<Particle*, std::array<dim_t,3>>> instructions;

/**
 * @class CellCalculator
 * @brief offers important functions for particle interactions.
 *
 * This class simulates the interaction between particles inside the cell storage structure.
 * It offers the functionality to calculate the forces, velocities and position of particles, by
 * requesting the next cells to consider, organised by the provided CellContainer, which allows
 * multiple CellCalculators to run in parallel based on the consumer-producer pattern.
 */
class CellCalculator {

public:
    CellCalculator(CellContainer &cellContainer, double delta_t, double cutoff, double r_l_,
                   std::array<boundary_conditions,6> boundaries_cond, std::string forceType,
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
     * and conquer with "paths" that are non overlapping for each pattern, which allows parallelisation
     * without using locks, with the hope of improving performance.
     */
    void calculateLinkedCellF();

    void calculatePeriodicF();

    /**
     * @brief applies a Thermostat iteration to the CellContainer of this CellCalculator
     * 
     *  First the current temperature @f$ T_{current} @f$ of the system (all Particles within the boundaries) is calculated. 
     *  Then the scaling factor @f$ \beta = \sqrt{ \frac{ T_{target} }{ T_{current} } } @f$ is calculated,
     *  which when applied to all particle velocities would change the temperature of the system 
     *  to `target_temp` (CellCalculator member). Then if a `max_temp_diff` is given, the absolute 
     *  value of @f$ \beta @f$ is capped by`max_temp_diff`. If capped the current temperature might 
     *  not be @f$ T_{target} @f$. Then the velocities of all particles are scaled by @f$ \beta @f$.
    */
    void applyThermostats();


    /**
     * @brief calculates the new positions for all particles
     *        and handles boundary_conditions::reflective
     *                    boundary_conditions::periodic
     *                    boundary_conditions::outflow
     *        by calling updateCells
    */
    void calculateX();


    /**
     * @brief calculates the velocities for all particles
    */
    void calculateV();


    /**
     * @brief calculates the forces for all particles by using the LinkedCell 
     *        Algorithm
    */
    void calculateF();


    /**
     * @brief shifts the forces of all particles
     * 
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

    concurrency_strategy parallelization;

    /**
     * @brief helper method to change the location of particles within the cell structure
     *
     * when the positions of the particles get changed a list of instructions is being created
     * that summarizes all the changes of particles between cells that have to be made. That
     * list is being processed in this method. 
     * This method then also takes care of:  boundary_conditions::reflective
     *                                       boundary_conditions::periodic
     *                                       boundary_conditions::outflow
     * by taking according actions, if the new position of a particle would be outside
     * the domain bounds
     *                  
     * 
     *
     * @param cell_updates list of instructions to change the location of particles
     */
    void updateCells(instructions& cell_updates);

    void applyBoundaries(Particle* particle_ptr, std::array<dim_t, 3>& new_cell_position, instructions& cell_updates);

    /**
     * @brief helper method to calculate the forces not covered in "calculateLinkedCellF()"
     *
     * since the "calculateLinkedCellF()" methods only calculates the forces between cells, we are
     * now calculating the forces between the particles within a cell.
     *
     * @param current_cell the cell to calculate the force within
     */
    void finishF(std::vector<Particle*> *current_cell);

    /**
     * @brief helper method to determine if the particles are in cutoff distance
     *
     * calculating the euclidean distance, but comparing the results squared to avoid calculating the square root
     */
    bool inCutoffDistance(Particle &p1, Particle &p2, const std::array<double,3> &offset) const;

    /**
     * @brief helper method to mirror the cell position on the other side
     *
     * used to implement the periodic boundaries, writes the changed distance of the
     * particles through the mirroring in the offset for further calculations
     */
    bool mirror(std::array<dim_t,3> &position, std::array<double,3> &offset);
};