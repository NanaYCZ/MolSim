#define SPDLOG_ACTIVE_LEVEL SPDLOG_LEVEL_TRACE

#include "particleModel/updating/CellCalculator.h"
#include "Simulation.h"
#include "outputWriter/XYZWriter.h"
#include <spdlog/spdlog.h>
#include <iostream>
#include <string>
#include <sstream>
#include <ostream>

#include <omp.h>
#include <spdlog/sinks/basic_file_sink.h>

void runSimulation(CellContainer &container, CellCalculator& calculator, ThermoStats &thermoStats,
                   const double end_time, const double delta_t, const size_t write_frequency, 
                   std::optional<int> thermostats_frequency, std::optional<int> diffusion_frequency,
                   std::optional<std::pair<double,int>> rdf_interval_and_frequency, bool performance_measurement) {

    outputWriter::VTKWriter writer;
    auto logger = spdlog::get("logger");

    std::chrono::high_resolution_clock::time_point perf_time_start, perf_time_end;

    double current_time = 0;
    int iteration = 0;

    std::string progressBar;
    size_t barWidth, pos = 0;

    SPDLOG_INFO("Starting Simulation");

    #ifdef _OPENMP
    if(calculator.parallelization == concurrency_strategy::first_method
    ||calculator.parallelization ==concurrency_strategy::second_method){
        SPDLOG_INFO("max threads: " + std::to_string(omp_get_max_threads()));
        SPDLOG_INFO("scheduled chunk size: " + std::to_string(chunk_size));
    }
    #endif
    calculator.calculateF();
    calculator.shiftF();

    SPDLOG_LOGGER_DEBUG(logger, "Particles in the simulation:");
    SPDLOG_LOGGER_DEBUG(logger, container.to_string());
    logger->flush();

    auto stat_logger = spdlog::basic_logger_mt("stat_log", "stat.txt");

    std::ostringstream rdf_log;
    std::ostringstream diff_log;
    std::ostringstream pot_log;
    std::ostringstream temp_log;

    // for this loop, we assume: current x, current f and current v are known
    if (performance_measurement)
        perf_time_start = std::chrono::high_resolution_clock::now();

    while (current_time < end_time) {
        SPDLOG_TRACE(std::to_string(current_time));
        SPDLOG_TRACE("Doing a Iteration with CellCalculator");
        //this applies Ghost Particle reflective boundary conditions, only 
        //if any of the boundaries is boundary_condition::ghost_reflective
        int x, f;
#pragma omp parallel if (calculator.parallelization == concurrency_strategy::second_method)
        {
#pragma omp master
            {
#pragma omp task depend(out:x)
                {calculator.calculateX();}
#pragma omp task depend(in:x) depend(out:f)
                {calculator.calculateF();}
#pragma omp task depend(in:f)
                {calculator.calculateV();}

                if (current_time<container.getSpecialTime()) calculator.calculateSpecialForce();
            }
        }
        
        iteration++;

        if (iteration % write_frequency == 0 && !performance_measurement) {
            writer.initializeOutput(container.size());
            container.plotParticles(writer);
            writer.writeFile("out", iteration);
        }

        calculator.shiftF();

        //thermostats_frequency.has_value() will be evaluated first
        if (thermostats_frequency.has_value() &&  iteration % thermostats_frequency.value() == 0) {
            thermoStats.applyThermostats();
        }
        
        if(diffusion_frequency.has_value() && iteration % diffusion_frequency.value() == 0){
            double diffusion = thermoStats.getDiffusionCoefficient();
            diff_log << "(" + std::to_string((iteration)/1000.0) + "," + std::to_string(diffusion) + ")\n";
            pot_log << std::to_string((iteration)/1000.0) + " " + std::to_string(thermoStats.getPotentialEnergy()) + "\n";
            //track the temperature as well
            double temp = thermoStats.currentTemp();
            temp_log << "(" + std::to_string(iteration/1000.0) + "," + std::to_string(temp) + ")\n";
        }

        if(rdf_interval_and_frequency.has_value() && iteration % (rdf_interval_and_frequency.value().second) == 0 ){
            double interval_size = rdf_interval_and_frequency.value().first;
            std::vector<double> stats = thermoStats.getRadialDistributionFunction(interval_size);
            rdf_log << "for time: " + std::to_string(current_time) + " iter:" + std::to_string(iteration) + "\n";
            for(size_t i = 0; i < stats.size();i++){
                rdf_log << "(" + std::to_string(i * interval_size + interval_size/2.0) + "," +
                            std::to_string(stats[i]) + ")";
            }
        }

        /// loading bar
        static int loading_factor = std::max(write_frequency * 5.0, std::ceil(end_time / (delta_t * 100)));

        if (iteration % loading_factor == 0 && !performance_measurement) {
            barWidth = 50;
            pos = static_cast<size_t>(barWidth * (current_time / end_time));
            progressBar = "[" + std::string(pos, '=') + '>'
                          + std::string(barWidth - pos, ' ') + "] "
                          + std::to_string(int((current_time / end_time) * 100.0)) + "%\r";

            spdlog::info(progressBar);
        }

        current_time += delta_t;
    }
    if (performance_measurement) {
        perf_time_end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> perf_duration = perf_time_end - perf_time_start;
        std::cout << "The Computation took: " << perf_duration.count() << " seconds" << std::endl;
    }
    if(diffusion_frequency.has_value()){
        stat_logger->info("diffusion:\n" + diff_log.str());
        stat_logger->info("temp:\n" + temp_log.str());
        stat_logger->info("potential_energy:\n" + pot_log.str());
    }

    if(rdf_interval_and_frequency.has_value())
        stat_logger->info("rdf:\n" + rdf_log.str());

    if(!performance_measurement)
        spdlog::info("[" + std::string(pos, '=') + ">] 100%\r");

    SPDLOG_INFO("output written. Terminating...\r");
}


void plotParticles(CellContainer &container, int iteration) {
  std::string out_name("MD_vtk");
  outputWriter::XYZWriter writer;
  writer.plotParticles(container, out_name, iteration);
}