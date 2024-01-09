#pragma once

#include "particleModel/storage/CellContainer.h"
#include <optional>

enum class boundary_conditions{
    outflow,
    reflective,
    ghost_reflective,
    periodic
};

extern double min_distance;

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
    CellCalculator(CellContainer &cellContainer, double delta_t, double cutoff,
                   std::array<boundary_conditions,6> boundaries_cond, std::optional<double> target_temp = std::nullopt, 
                   std::optional<double> max_temp_diff_param = std::nullopt,
                   double gravity_factor = 0);

    /**
     * @brief initializes the force and updates the position to remain the calculation order
     *
     * we set the calculation order inside the while loop to F - V - X, that's why we have to
     * calculate X before the loop in order to keep the pattern F - X - F - V - X - F - V - X...
     */
    void initializeFX();

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

    /**
     * @brief calculate the forces not covered in "calculateLinkedCellF()", V, X and update the cells
     *
     * We applied a new order of calculating because V and X can be calculated independently from
     * other particles as long as the force of them is all up-to-date. This means that after we
     * finished the force calculation within a cell, we can keep the cell and calculate V and X
     * on top. This is more elegant than having to iterate over all cells for X, then F + cell pairs
     * and then iterate over all cells for V again.
     */
    void calculateWithinFVX();

    /**
     * 
     * @brief calculates reflective boundary conditions with Ghost Particles
     * 
     * ! only acts on Particles if any boundary is boundary_condition::ghost_reflective !
     * 
     * This method assumes, that at the moment it is called all particles are within the 
     * domain boundaries. If a particle has coordinates outside the domain it is undefined behaviour.
     * applyBoundaries iterates over every side of the cuboid (that represents the valid domain) and 
     * applies boundary conditions the boundary conditions specified in the 'boundaries' member
     * of the CellCalculator class.
     * If we are not simulating with three dimensions, the claculation of boundary conditions
     * for the top and bottom of the domain do not happen (positive/ negative z direction).
     * 
     * If the boundary conditions are 'outflow', then just no 
     * boundary conditions are applied in that direction.
     * The outflow is automatically realized by updateCells, which deletes
     * Particles from the cellContainer, if they have moved outside the boundary.
     * (the particles are moved into the halo_particles vector)
     * 
     * If the boundary conditions are reflective, applyReflectiveBoundaries uses Ghost Partilces to create
     * opposing forces for particles in the direction of that boundary side. 
     * 
     * (for imagination(how i think):
     * 
     *  z ^   ^
     *    |  / x
     *    | /
     *    |/_________ > y
     * )
     * 
     * 
     * 
    */
    void applyReflectiveBoundaries();

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

    std::array<double,3> force(const Particle &p_i, const Particle &p_j, const std::array<double,3> &offset);

    std::array<double,3> ghostParticleLennardJonesForce(const Particle &particle,std::array<double,3> ghost_position);


private:
    CellContainer &cellContainer;
    const double gravity_factor;
    const double delta_t;
    const double cutoff;
    std::array<double,3> domain_bounds;
    std::array<dim_t, 3> domain_max_dim;

    std::optional<double> target_temp;
    std::optional<double> max_temp_diff;

    //{positive_z,negative_z,positive_x,negative_x,positive_y,negative_y}
    std::array<boundary_conditions,6> boundaries;
    bool ghost_reflection_is_off;

    std::vector<std::vector<std::vector<std::vector<Particle*>>>> &particles;

    void addGhostParticleForcesInDir_i(int i,double boundary,
                                    std::vector<Particle*>& particles);


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

    /**
     * @brief helper method calculate the Velocity and Position
     *
     * in order to keep "calculateWithinFVX()" and "initializeFX()" clean, we moved
     * the functionality to calculate the two attributes V and X into this method.
     *
     * @param particle the particle to update the attributes of
     * @param calculateV to make it applicable for both calculations "FVX" and "FX"
     */
    void calculateVX(Particle& particle, bool calculateV);

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