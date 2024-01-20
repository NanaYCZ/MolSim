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
 * using a lambda function interface for the calcualtion of the force between two
 * Particles, takes in two Particles and returns the vector3 of forces
 * acting between the two given Particles
 * simplified:
 * forceCalculation refers to such functions "std::array<double,3> func(const Particle&,const Particle&)"
 * uses constant references because forceCalculation mustn't change the Particles
 */
using ForceCalculation = std::function<std::array<double, 3>(const Particle &, const Particle &, const std::array<double,3> &)>;


/**
 * @brief Calculate force between \f$ p_i \f$ and \f$ p_j \f$
 *
 * Uses particle \f$ p_i \f$ and particle \f$ p_j \f$ to calculate the
 * force \f$ f_{ij} \f$ between them.
 * Uses the formula from the first worksheet:
 * \f$ F_{ij} = \frac{m_i m_j}{{\| \mathbf{x}_i - \mathbf{x}_j \|}^3}
 * (\mathbf{x}_j - \mathbf{x}_i) \f$
 *
 *
 * @param p_i Particle \f$ i \f$ for force calculation
 * @param p_j Particle \f$ j \f$ for force calculation
 *
 * @return Three-dimensional vector that corresponds to \f$ f_{ij} \f$
 */
ForceCalculation inline forceSimpleGravitational(double cutoff){
        return [cutoff](const Particle &p_i, const Particle &p_j, const std::array<double,3> &offset) -> std::array<double,3> {
        double m_i = p_i.getM(), m_j = p_j.getM();

        std::array<double, 3> x_i = p_i.getX(), x_j = p_j.getX();

        double r_c = cutoff;

        double dx = x_i[0] - x_j[0] + offset[0];
        double dy = x_i[1] - x_j[1] + offset[1];
        double dz = x_i[2] - x_j[2] + offset[2];

        double r_c_squared = r_c * r_c;
        double scalar_product = dx * dx + dy * dy + dz * dz;

        /*instantly return 0 if r_c <= norm */
        if(r_c_squared <= scalar_product)
            return {0,0,0};
        else{

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
 *
 * @param p_i Particle \f$ i \f$ for force calculation
 * @param p_j Particle \f$ j \f$ for force calculation
 *
 * @return Three-dimensional vector that corresponds to \f$ f_{ij} \f$
 */
ForceCalculation inline forceSmoothedLennJonesPotentialFunction(std::vector<std::vector<double>>& sigma_mixed,
                                                 std::vector<std::vector<double>>& epsilon_mixed, double cutoff, double r_l){
    return [&sigma_mixed,&epsilon_mixed,cutoff,r_l]
        (const Particle &p_i, const Particle &p_j, const std::array<double,3> &offset) -> std::array<double,3> {
              const auto& x_i = p_i.getX(), x_j = p_j.getX();
    double sigma = sigma_mixed[p_i.getType()][p_j.getType()];
    double epsilon = epsilon_mixed[p_i.getType()][p_j.getType()];
    //make formula more readable, compiler will optimize away
    double r_c = cutoff;

    double dx = x_i[0] - x_j[0] + offset[0];
    double dy = x_i[1] - x_j[1] + offset[1];
    double dz = x_i[2] - x_j[2] + offset[2];

    double r_c_squared = r_c * r_c;
    double scalar_product = dx * dx + dy * dy + dz * dz;

    /*instantly return 0 if r_c <= norm */
    if(r_c_squared <= scalar_product)
        return {0,0,0};

    double norm = std::sqrt(scalar_product);
    //norm = std::max(min_distance, norm);

    /* r_l <= norm  < r_c */
    if(r_l <= norm){
        double norm_pow_6 = std::pow(norm,6), sigma_pow_6  = std::pow(sigma,6);

        double prefactor = -(24 * sigma_pow_6 * epsilon * (r_c - norm)) / 
                            (norm_pow_6 * norm_pow_6 * norm * norm  * std::pow((r_c - r_l),3) );
        
        prefactor *=    r_c_squared * (2 * sigma_pow_6 - norm_pow_6)
                         + r_c * (3 * r_l - norm) * (norm_pow_6 - 2 * sigma_pow_6)
                         + norm * (5 * r_l * sigma_pow_6 - 2 * r_l * norm_pow_6 
                                  - 3 * sigma_pow_6 * norm + norm_pow_6 * norm);
        
        return prefactor * -1 * (x_i - x_j + offset);
    }else /* norm <= r_l */ {
        double prefactor = (-24 * epsilon) / (std::pow(norm, 2));

        prefactor *= (std::pow(sigma / norm, 6) - 2 * std::pow(sigma / norm, 12));

        return prefactor * (x_i - x_j + offset);
    }

    };
}


ForceCalculation inline forceLennJonesPotentialFunction(std::vector<std::vector<double>>& sigma_mixed,
                                                 std::vector<std::vector<double>>& epsilon_mixed, double cutoff){
    return forceSmoothedLennJonesPotentialFunction(sigma_mixed,epsilon_mixed,cutoff,cutoff);
}




