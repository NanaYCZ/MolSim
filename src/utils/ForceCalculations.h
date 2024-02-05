#pragma once

#include <array>
#include <functional>

#include "particleModel/storage/Particle.h"
#include "utils/ArrayUtils.h"

/**
 * @file
 * 
 * contains the different functions that can be used for calculation of
 * forces between two concrete Particles, all functions share the same
 * function interface, such that they can be passed as arguments
 * to the calulateF() function, that calculates Forces between not only two, but
 * all Particles
 */



/**
 * using a lambda function interface for the Calculation of the force between two
 * Particles, takes in two Particles, an offset(which is relevant for periodic boundaries) and returns the vector3 of forces
 * acting between the two given Particles
 * simplified:
 * forceCalculation refers to such functions:
 * "std::array<double,3> func(const Particle&,const Particle&,const std::array<double,3> offset)"
 * uses constant references because forceCalculation mustn't change the Particles
 */
using ForceCalculation = std::function<std::array<double, 3>(const Particle &, const Particle &, const std::array<double,3> &)>;


/**
 * @brief Calculate force between neighboring particles
 *
 * Uses particle grids (the index of the particles in three dimensions) to determine if they are neighbors
 * If true, uses the harmonic potential to calculate, which further sets different factor of multiplication
 * depending on whether they are direct neighbors or diagonal neighbors.
 * If false, uses LJ force to calculate, which will be set to zero if the length is too large to save calculation
 * or if length is too small to avoid self penetration of the membrane. In the end only repulsive part will remain.
 *
 * @param sigma_mixed the matrix of mixed interaction parameters for the respective Particles (supplied by CellCalculator)
 * @param epsilon_mixed the matrix of mixed interaction parameters for the respective Particles (supplied by CellCalculator)
 * @param cutoff upper cutoff radius for the force
 *
 * @return Three-dimensional vector that corresponds to \f$ f_{ij} \f$
**/

ForceCalculation inline forceHarmonicForce (std::vector<std::vector<double>>& sigma_mixed,
                                            std::vector<std::vector<double>>& epsilon_mixed, double cutoff) {
    return [&sigma_mixed,&epsilon_mixed,cutoff](const Particle &p_i, const Particle &p_j, const std::array<double, 3> &offset) -> std::array<double, 3> {
//direct neighbor factor=1, diagonal neighbor factor=2
        double factor = sqrt(2);
//r_zero average bond length
        const auto r_zero = p_i.getRZ();
//position of two particles
        const auto &x_i = p_i.getX(), x_j = p_j.getX();
//position index of two particles
        const auto &grid_i = p_i.getGrid(), grid_j = p_j.getGrid();
//distance of the two particles
        std::array<double, 3> delta_x = x_j - x_i + offset;
        double scalar_product = ArrayUtils::scalarProduct(delta_x, delta_x);
        double norm = std::sqrt(scalar_product);
//direct neighbor
        if (abs(grid_i[0] - grid_j[0]) + abs(grid_i[1] - grid_j[1])
            + abs(grid_i[2] - grid_j[2]) == 1) {
            factor = 1;
            return p_i.getFP() * (norm - factor*r_zero) / norm * delta_x;
        }
//diagonal neighbor
        else if (abs(grid_i[0] - grid_j[0]) <= 1 && abs(grid_i[1] - grid_j[1]) <= 1
        && abs(grid_i[2] - grid_j[2]) <= 1){
            return p_i.getFP() * (norm - factor*r_zero) / norm* delta_x;
        }
//not neighbor, use lj force
        else{
            double sigma = sigma_mixed[p_i.getType()][p_j.getType()];
            double epsilon = epsilon_mixed[p_i.getType()][p_j.getType()];
            double r_c_squared = cutoff * cutoff;
 //norm>=cutoff, force = 0
            if(r_c_squared <= scalar_product) return {0,0,0};
 //norm>=2^1/6*sigma, force = 0, to avoid self penetration
            if(norm>sigma*1.12246204831) return {0,0,0};
            double prefactor = (-24 * epsilon) / (std::pow(norm, 2));
            prefactor *= (std::pow(sigma / norm, 6) - 2 * std::pow(sigma / norm, 12));
            return prefactor * (x_i - x_j + offset);
        }

    };
}

/**
 * @brief Calculate force between \f$ p_i \f$ and \f$ p_j \f$
 *
 * Uses particle \f$ p_i \f$ and particle \f$ p_j \f$ to calculate the
 * force \f$ f_{ij} \f$ between them.
 * Uses the formula from the first worksheet:
 * \f$ F_{ij} = \frac{m_i m_j}{{\| \mathbf{x}_i - \mathbf{x}_j \|}^3}
 * (\mathbf{x}_j - \mathbf{x}_i) \f$
 *
 * (offset is always applied, when the difference between the positions of p_i and p_j is calculated)
 *
 * @param sigma_mixed the matrix of mixed interaction parameters for the respective Particles (supplied by CellCalculator)
 * @param epsilon_mixed the matrix of mixed interaction parameters for the respective Particles (supplied by CellCalculator)
 * @param cutoff upper cutoff radius for the force
 *
 * @return Three-dimensional vector that corresponds to \f$ f_{ij} \f$
 */

ForceCalculation inline forceSimpleGravitational(double cutoff){
        return [cutoff](const Particle &p_i, const Particle &p_j, const std::array<double,3> &offset) -> std::array<double,3> {
        double m_i = p_i.getM(), m_j = p_j.getM();
        std::array<double, 3> x_i = p_i.getX(), x_j = p_j.getX();

        double r_c = cutoff;
        double r_c_squared = r_c * r_c;

        std::array<double, 3> delta_x = x_i - x_j + offset;
        double scalar_product = ArrayUtils::scalarProduct(delta_x,delta_x);

        /*instantly return 0 if r_c <= norm */
        if(r_c_squared <= scalar_product) {
            return {0, 0, 0};
        }else{
            double prefactor = (m_i * m_j) / std::pow(std::sqrt(scalar_product), 3);
            return prefactor * -1 * (x_i - x_j + offset);
        }
  };
}


/**
 * @brief Calculate force between \f$ p_i \f$ and \f$ p_j \f$
 *
 * Uses particle \f$ p_i \f$ and particle \f$ p_j \f$ to calculate the
 * force \f$ f_{ij} \f$ between them.
 * Uses the formula from the second worksheet:
 * \f$ F_{ij} = -\frac{24 \varepsilon}{{\| \mathbf{x}_i - \mathbf{x}_j \|}^2} \left(
 * \left( \frac{\sigma}{{\| \mathbf{x}_i - \mathbf{x}_j \|}} \right)^{6} - 2
 * \left( \frac{\sigma}{{\| \mathbf{x}_i - \mathbf{x}_j \|}} \right)^{12}
 * \right) (\mathbf{x}_i - \mathbf{x}_j) \f$
 *
 * (offset is always applied, when the difference between the positions of p_i and p_j is calculated)
 *
 *
 * @param sigma_mixed the matrix of mixed interaction parameters for the respective Particles (supplied by CellCalculator)
 * @param epsilon_mixed the matrix of mixed interaction parameters for the respective Particles (supplied by CellCalculator)
 * @param cutoff upper cutoff radius for the force
 *
 * @return Three-dimensional vector that corresponds to \f$ f_{ij} \f$
 */
ForceCalculation inline forceLennJonesPotentialFunction(std::vector<std::vector<double>>& sigma_mixed,
                                                        std::vector<std::vector<double>>& epsilon_mixed, double cutoff){
    return [&sigma_mixed,&epsilon_mixed,cutoff]
            (const Particle &p_i, const Particle &p_j, const std::array<double,3> &offset) -> std::array<double,3> {
        const auto& x_i = p_i.getX(), x_j = p_j.getX();
        double sigma = sigma_mixed[p_i.getType()][p_j.getType()];
        double epsilon = epsilon_mixed[p_i.getType()][p_j.getType()];
        //make formula more readable, compiler will optimize away
        double r_c = cutoff;
        double r_c_squared = r_c * r_c;

        std::array<double, 3> delta_x = x_i - x_j + offset;
        double scalar_product = ArrayUtils::scalarProduct(delta_x,delta_x);

        /*instantly return 0 if r_c <= norm */
        if(r_c_squared <= scalar_product)
            return {0,0,0};

        double norm = std::sqrt(scalar_product);

        double prefactor = (-24 * epsilon) / (std::pow(norm, 2));

        prefactor *= (std::pow(sigma / norm, 6) - 2 * std::pow(sigma / norm, 12));

        return prefactor * (x_i - x_j + offset);

    };
}


/**
 * @brief returns a function, that calculates the smoothed Lennard-Jones-Potential according to this formula:
 *
 *  if \f$ \|\mathbf{x}_i - \mathbf{x}_j\| \geq r_c \f$ then \f$ F(\mathbf{x}_i, \mathbf{x}_j) = 0\f$
 *
 *  if \f$ r_l < \|\mathbf{x}_i - \mathbf{x}_j\| < r_c \f$ then \f$ F(\mathbf{x}_i, \mathbf{x}_j) =
 *     -24 \times \frac{\sigma^6\varepsilon(r_c - \|\mathbf{x}_i - \mathbf{x}_j\|)}{\|\mathbf{x}_i - \mathbf{x}_j\|^14 (r_c - r_l)^3} \times \\
 *     \quad \left( r_c^2(2\sigma^6 - \|\mathbf{x}_i - \mathbf{x}_j\|^6) + r_c(3r_l - \|\mathbf{x}_i - \mathbf{x}_j\|)(\|\mathbf{x}_i - \mathbf{x}_j\|^6 - 2\sigma^6) + \right. \\
 *     \quad \left. \|\mathbf{x}_i - \mathbf{x}_j\|(5r_l\sigma^6 - 2r_l\|\mathbf{x}_i - \mathbf{x}_j\|^6 - 3\sigma^6\|\mathbf{x}_i - \mathbf{x}_j\| + \|\mathbf{x}_i - \mathbf{x}_j\|^7)\right) \times \\
 *     \quad (\mathbf{x}_i - \mathbf{x}_j)
 * \f$
 *
 * and if \f$ \|\mathbf{x}_i - \mathbf{x}_j\| \leq r_l \f$ then  \f$ F(\mathbf{x}_i, \mathbf{x}_j) =
 *  -24 \times \frac{ \varepsilon}{{\| \mathbf{x}_i - \mathbf{x}_j \|}^2} \left(
 * \left( \frac{\sigma}{{\| \mathbf{x}_i - \mathbf{x}_j \|}} \right)^{6} - 2
 * \left( \frac{\sigma}{{\| \mathbf{x}_i - \mathbf{x}_j \|}} \right)^{12}
 * \right) (\mathbf{x}_i - \mathbf{x}_j) \f$
 *
 * with the following definitions:
 * \f$ \( \mathbf{x}_i \) \f$ and \f$ \( \mathbf{x}_j \) \f$ are the position vectors of particles \f$ \( i \) \f$ and \f$ \( j \) \f$
 * \f$ \( \sigma \) \f$ and \f$ \( \varepsilon \) \f$ are mixed sigma / epsilon depending on the particles p_i and p_j
 * \f$ \( r_l \) \f$ is the lower cutoff radius
 * \f$ \( r_c \) \f$ is the upper cutoff radius
 *
 * (offset is always applied, when the difference between the positions of p_i and p_j is calculated)
 *
 * @param sigma_mixed the matrix of mixed interaction parameters for the respective Particles (supplied by CellCalculator)
 * @param epsilon_mixed the matrix of mixed interaction parameters for the respective Particles (supplied by CellCalculator)
 * @param cutoff upper cutoff radius for the force
 * @param r_l lower cutoff radius for the force
 *
 * @return Three-dimensional vector corresponding to the force \( f_{ij} \)
 */
ForceCalculation inline forceSmoothedLennJonesPotentialFunction(std::vector<std::vector<double>>& sigma_mixed,
                                                 std::vector<std::vector<double>>& epsilon_mixed, double cutoff, double r_l){
    return [&sigma_mixed,&epsilon_mixed,cutoff,r_l]
        (const Particle &p_i, const Particle &p_j, const std::array<double,3> &offset) -> std::array<double,3> {
              const auto& x_i = p_i.getX(), x_j = p_j.getX();
    double sigma = sigma_mixed[p_i.getType()][p_j.getType()];
    double epsilon = epsilon_mixed[p_i.getType()][p_j.getType()];

    //make formula more readable, compiler will optimize
    double r_c = cutoff;
    double r_c_squared = r_c * r_c;

    std::array<double, 3> delta_x = x_i - x_j + offset;
    double scalar_product = ArrayUtils::scalarProduct(delta_x,delta_x);

    /*instantly return 0 if r_c <= norm (sqrt not calculated !)*/
    if(r_c_squared <= scalar_product)
        return {0,0,0};

    double norm = std::sqrt(scalar_product);


    /* r_l <= norm  < r_c */
    if(r_l < norm){
        double norm_pow_6 = std::pow(norm,6), sigma_pow_6  = std::pow(sigma,6);


        //corresponds to -(24 * sigma^6 * epsilon * (r_c - norm)) / (norm^14 * (r_c - r_l)^3)
        //but makes some optimization to avoid unnecessary expensive powers
        double prefactor = -(24 * sigma_pow_6 * epsilon * (r_c - norm)) / 
                            (norm_pow_6 * norm_pow_6 * norm * norm  * std::pow((r_c - r_l),3) );

        //corresponds to r_c^2 * (2 * sigma^6 - norm^6)
        //             + r_c * (3 * r_l - norm) (norm^6 - 2 * sigma^6)
        //             + norm * (5 * r_l * sigma^6 - 2 * r_l *  norm^6 - 3 * sigma^6 * norm + norm^7)
        prefactor *=    r_c_squared * (2 * sigma_pow_6 - norm_pow_6)
                         + r_c * (3 * r_l - norm) * (norm_pow_6 - 2 * sigma_pow_6)
                         + norm * (5 * r_l * sigma_pow_6 - 2 * r_l * norm_pow_6 
                                  - 3 * sigma_pow_6 * norm + norm_pow_6 * norm);
        
        return prefactor * -1 * (x_i - x_j + offset); //make clear that here the x_j - x_i is returned
                                                      //but still the offset is applied in the correct direction,
                                                      //by multiplying by -1
    }else /* norm <= r_l */ {
        double prefactor = (-24 * epsilon) / (std::pow(norm, 2));

        prefactor *= (std::pow(sigma / norm, 6) - 2 * std::pow(sigma / norm, 12));

        return prefactor * (x_i - x_j + offset);
    }

    };
}





