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

ForceCalculation inline forceBetweenTwoNeighbors(){
    return [] (const Particle &p_i, const Particle &p_j, const std::array<double,3> &offset) -> std::array<double, 3> {


        auto delta_x = p_i.getX() - p_j.getX();
        double l2Norm = std::sqrt(delta_x[0] * delta_x[0] + delta_x[1] * delta_x[1] + delta_x[2] * delta_x[2]);

        double k = p_i.getFP();
        double r0 = p_i.getA();

        std::array<double, 3> force;
        double magnitude = k * (l2Norm - r0) / l2Norm;

        for (int dim = 0; dim < 3; ++dim) {
            force[dim] = magnitude * (p_j.getX()[dim] - p_i.getX()[dim]);
        }

        return force;
    };
}


ForceCalculation inline forceBetweenTwoDiagonalNeighbors(){
    return [] (const Particle &p_i, const Particle &p_j, const std::array<double,3> &offset) -> std::array<double, 3> {


        auto delta_x = p_i.getX() - p_j.getX();
        double l2Norm = std::sqrt(delta_x[0] * delta_x[0] + delta_x[1] * delta_x[1] + delta_x[2] * delta_x[2]);

        double k = p_i.getFP();
        double r0 = p_i.getA();

        std::array<double, 3> force;
        double magnitude = k * (l2Norm - sqrt(2)*r0) / l2Norm;

        for (int dim = 0; dim < 3; ++dim) {
            force[dim] = magnitude * (p_j.getX()[dim] - p_i.getX()[dim]);
        }

        return force;
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
 * @param sigma_mixed Matrix of mixed interaction parameters for the respective Particles (supplied by CellCalculator)
 * @param epsilon_mixed Matrix of mixed interaction parameters for the respective Particles (supplied by CellCalculator)
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
 * @param sigma_mixed Matrix of mixed interaction parameters for the respective Particles (supplied by CellCalculator)
 * @param epsilon_mixed Matrix of mixed interaction parameters for the respective Particles (supplied by CellCalculator)
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
        std::cout << "  A A  A";
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
 * @param sigma_mixed Matrix of mixed interaction parameters for the respective Particles (supplied by CellCalculator)
 * @param epsilon_mixed Matrix of mixed interaction parameters for the respective Particles (supplied by CellCalculator)
 * @param cutoff upper cutoff radius for the force
 * @param r_l Lower cutoff radius for the force
 *
 * @return Three-dimensional vector corresponding to the force \( F_{ij} \)
 */
ForceCalculation inline forceSmoothedLennJonesPotentialFunction(std::vector<std::vector<double>>& sigma_mixed,
                                                 std::vector<std::vector<double>>& epsilon_mixed, double cutoff, double r_l){
    return [&sigma_mixed,&epsilon_mixed,cutoff,r_l]
        (const Particle &p_i, const Particle &p_j, const std::array<double,3> &offset) -> std::array<double,3> {
              const auto& x_i = p_i.getX(), x_j = p_j.getX();
    double sigma = sigma_mixed[p_i.getType()][p_j.getType()];
    double epsilon = epsilon_mixed[p_i.getType()][p_j.getType()];
    //make formula more readable, compiler will optimize away
    std::cout << "  A A  A";
    double r_c = cutoff;
    double r_c_squared = r_c * r_c;

    std::array<double, 3> delta_x = x_i - x_j + offset;
    double scalar_product = ArrayUtils::scalarProduct(delta_x,delta_x);

    /*instantly return 0 if r_c <= norm */
    if(r_c_squared <= scalar_product)
        return {0,0,0};

    double norm = std::sqrt(scalar_product);
    //norm = std::max(min_distance, norm);

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





