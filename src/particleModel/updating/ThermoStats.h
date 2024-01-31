

#include "particleModel/storage/CellContainer.h"
#include "particleModel/updating/CellCalculator.h"
#include <optional>


class ThermoStats{
    public:

        ThermoStats(CellContainer& container, double delta_t_param,
                std::optional<double> target_temp_param = std::nullopt,
                std::optional<double> max_temp_diff_param = std::nullopt);


        double currentTemp();

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
         * @brief calculates the potential energy of the system according to:
         * for all unique particle pairs i,j:
         * E_pot += 4 * epsilon * (
         *      (sigma/ || particle_i position - particle_j position ||)^12
         *     -(sigma/ || particle_i position - particle_j position ||)^6
         *     )
         *
         * @return potential energy of the system
         */
        double getPotentialEnergy();

        double getPressure();

        /**
         * @brief calculates radial distribution function for the particles currently in CellContainer
         *
         * Iterates over all possible unique((p1,p2) is same as (p2,p1)) particle pairs calculates the distance
         * and stores it (by incrementing a counter in the respective interval, in which the distance falls).
         * Then we have a set of intervals afterwards and for every interval a counter of how much particle
         * distances in the range of that interval, there are. Then for every interval, the radial distribution
         * function is calculated according to the formula:
         * rdf in interval i = (amt particle pair distances in interval i) / (
         *                                  (4pi/3) * ( (upper bound interval of i)^3 - (lower bound interval of i)^3 )
         *                                  )
         *
         * @param interval_size determines the granularity of "buckets", in which different
         *                      particle distances are tracked
         *
         * @return a vector, where every entry represents an the value of the radial distr. function
         *         for an interval
         */
        std::vector<double> getRadialDistributionFunction(double interval_size);

        /**
         *  @brief Calculates the diffusion coefficient by summing up over (current_position - old_position)
         *         for all particles that are currently in the CellContainer and dividing by 1 / (amount of particles).
         *         Where old_position for a particle is retrieved from particle_positions_previous_iteration.
         *         It therefore calculates the movement of the particles with respect to the the reference time,
         *         at which the particle positions were stored the last time. After that this function calls
         *         initDiffusionCoefficient in order to store the positions for the next iteration, at which
         *         the getDiffusionCoefficient will be called again
         *
         * @return diffusion coefficient
         */
        double getDiffusionCoefficient();

        /**
         * @brief stores the particle positions of all particles in CellContainer
         *        as a list of (particle_pointer,position) pairs in the
         *        particle_positions_previous_iteration member of this class
         */
        void initDiffusionCoefficient();

        auto& getParticlePositionsPreviousIteration(){
            return particle_positions_previous_iteration;
        }


private:
        CellContainer& cellContainer;

        double delta_t;

        std::optional<double> target_temp;
        std::optional<double> max_temp_diff;

        std::list<std::pair<Particle*,std::array<double,3>>> particle_positions_previous_iteration = {};

};