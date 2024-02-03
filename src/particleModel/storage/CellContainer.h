#pragma once

#include "outputWriter/VTKWriter.h"
#include <vector>

/**
 * @brief type for dimensions, used to apply changes easily in the code
 */
typedef int dim_t;

/**
 * @class CellContainer
 * @brief container structure to store the particles within cells
 */
class CellContainer {
public:
    /**
     * @brief constructor initializing the cell storage structure and the maximum domain cell index
     *
     * the cell storage structure contains the halo cells. Also the dimension will be set for the
     * linked cell algorithm and the comparing depth is being determined based on the cell size and
     * cutoff radius.
     *
     * precalculates patterns based on comparing depth, for CellCalculator to use
     *
     * @param domain_width width of the domain
     * @param domain_height height of the domain
     * @param domain_depth depth of the domain, set smaller than cell size to create 2d structure
     * @param r_cutoff cutoff to indicate when particles too far away can be ignored
     * @param cell_size should be either bigger than the cutoff or the cutoff should be a multiple of the cell_size
     *
     * @throws std::invalid_argument in case the comparing depth is too big for the domain size,
     * in that case the cell size needs to be adjusted for the linked cell algorithm
     */
    CellContainer(double domain_width, double domain_height, double domain_depth, double r_cutoff, double cell_size);

    virtual ~CellContainer();
    /**
     * @brief creates a particle instance with the given parameters
     *
     * stores the particle in particle_instances, and increments particle_amount.
     *
     * @param x_arg particle position
     * @param v_arg particle velocity
     * @param m_arg particle mass
     *
     * @throws std::invalid_argument if the particle position is outside the domain bounds
     */
    void addParticle(std::array<double, 3> x_arg, std::array<double, 3> v_arg, double m_arg);

    void addParticle(std::array<double, 3> x_arg, std::array<double, 3> v_arg, std::array<int,3> grid_index, double m_arg);

    void addParticle(std::array<double, 3> x_arg, std::array<double, 3> v_arg, double m_arg, double sigma, double epsilon);

    void addParticle(std::array<double, 3> x_arg, std::array<double, 3> v_arg, double m_arg, std::array<int,3> grid_index,  double sigma, double epsilon);

    void addParticle(const Particle& particle,double sigma, double epsilon);

    /**
     * @brief after all Particles were created and are stored, this function creates pointers to them
     * 
     * In the actual cells no Particles, but particle pointer are stored, because the particles
     * need to be frequently moved to other cells and moving pointers is cheaper.
    */
    void createPointers();

    void plotParticles(outputWriter::VTKWriter &writer);

    /**
     * @brief returns the string representation of the CellContainer
     * @returns string representation
    */
    std::string to_string();

    /**
     * @brief returns size of the CellContainer
     * @returns size of the CellContainer
    */
    size_t size();

    std::list<Particle> to_list();

    /**
     * Getter(we did not anotate each of those :3):
    */
    std::vector<Particle>& getInstances() {
        return particle_instances;
    }

    bool hasThreeDimensions(){
        return three_dimensions;
    }

    std::array<double,3> getDomainBounds(){
        return domain_bounds;
    }

    std::vector<std::array<dim_t,3>> getPatterns() {
        return patterns;
    }

    /**
    * @brief allocates a cell position in the domain for the given position
    *
    * @param x particle position to map a cell position to
    * @param cell_position array to write the results into
    */
    void allocateCellFromPosition(const std::array<double, 3> &x, std::array<dim_t , 3> &cell_position) {
        cell_position[0] = std::floor(x[0] / cell_size + 1);
        cell_position[1] = std::floor(x[1] / cell_size + 1);
        cell_position[2] = std::floor(x[2] / cell_size + 1);

        //cover edge case with last cell being less than cell_size
        if(domain_bounds[0] < x[0]) {
            cell_position[0] = std::ceil((x[0] - domain_bounds[0]) / cell_size + domain_max_dim[0]);
        }
        if(domain_bounds[1] < x[1]) {
            cell_position[1] = std::ceil((x[1] - domain_bounds[1]) / cell_size + domain_max_dim[1]);
        }
        if(domain_bounds[2] < x[2]) {
            cell_position[2] = std::ceil((x[2] - domain_bounds[2]) / cell_size + domain_max_dim[2]);
        }
    }

private:
    bool three_dimensions;
    const double cell_size;
    std::array<double, 3> domain_bounds;
    dim_t comparing_depth = 1;
    size_t particle_amount = 0;
    std::array<int,3> default_grid_index ={0, 0, 0};

    std::vector<Particle> particle_instances;

    static std::array<dim_t, 3> domain_max_dim;
    static std::vector<std::array<dim_t,3>> patterns;
    static std::vector<std::vector<std::vector<std::vector<Particle*>>>> particles;

    friend class ThermoStats;
    friend class CellCalculator;
    friend class CellIterator;
    friend class StartPointIterator;

    /**
      * @brief helper method to set the next 3d pattern in the pre-calculation
      *
      * reuse of legacy code
      *
      * iterates in total over 3 sets of directions, which together make up the
      * half of all possible directions from a starting point, since we are applying
      * newtons 3rd law and have to traverse in only one direction.
      *
      * @param pattern to store the next pattern iteration
      * @return bool indicating that last pattern is not reached
      */
    bool setNext3dPattern(std::array<dim_t, 3> &pattern);

    /**
      * @brief helper method to set the next 2d pattern in the pre-calculation
      *
      * reuse of legacy code
      *
      * iterates in total over 2 sets of directions, which together make up the
      * half of all possible directions from a starting point, since we are applying
      * newtons 3rd law and have to traverse in only one direction.
      *
      * example of the two sets:
      *
      *         - - - - -
      *         - - - - -
      *         - - x 1 1
      *         2 2 2 2 2
      *         2 2 2 2 2
      *
      * @param pattern to store the next pattern iteration
      * @return bool indicating that last pattern is not reached
      */
    bool setNext2dPattern(std::array<dim_t, 3> &pattern);


    /**
      * @brief helper method to check result of double modulo
      *
      * @param a first value
      * @param b second value
      * @param epsilon threshold
      */
    bool isApproximatelyEqual(double a, double b, double epsilon = 1e-8);
};
