

#include "particleModel/storage/CellContainer.h"
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

        double getPotentialEnergy();

        double getPressure();


        std::vector<double> radialDistributionFunction(double interval_size);

        double diffusionCoeff();

        void initDiffCoeff();


    private:
        CellContainer& cellContainer;

        double delta_t;

        std::optional<double> target_temp;
        std::optional<double> max_temp_diff;

        std::list<std::pair<Particle*,std::array<double,3>>> old_particle_positions = {};

};