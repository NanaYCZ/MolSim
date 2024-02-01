#define SPDLOG_ACTIVE_LEVEL SPDLOG_LEVEL_TRACE

#include "particleModel/CellCalculator.h"
#include "Simulation.h"
#include "outputWriter/VTKWriter.h"
#include <spdlog/spdlog.h>
#include <iostream>
#include <string>
#include <sstream>
#include <ostream>
#include <omp.h>
#include <spdlog/sinks/basic_file_sink.h>
#include <optional>



void runSimulation(CellContainer &container, calculator::CellCalculator& calculator, outputWriter::VTKWriter vtkWriter, CellContainer::concurrency_strategy strategy, double t_end,double delta_t,double special_t, size_t write_frequency,bool performance_measurement) {

    auto logger = spdlog::get("logger");

    std::chrono::high_resolution_clock::time_point perf_time_start, perf_time_end;

    double current_time = 0;
    int iteration = 0;

    std::string progressBar;
    size_t barWidth, pos = 0;

    SPDLOG_INFO("Starting Simulation");
    if (strategy != CellContainer::serial) {
        SPDLOG_INFO("max threads: " + std::to_string(omp_get_max_threads()));
    }
    container.setup();
    container.setup();
    calculator.calcF(container);

    SPDLOG_LOGGER_DEBUG(logger, "Particles in the simulation:");

    logger->flush();

    auto stat_logger = spdlog::basic_logger_mt("stat_log", "stat.txt");

    std::string rdf_log = "";
    std::string diff_log = "";
    std::string pot_log = "";
    std::string pres_log = "";
    std::string temp_log = "";


    //size_t before_size = container.size();    

    // for this loop, we assume: current x, current f and current v are known
    if (performance_measurement)
        perf_time_start = std::chrono::high_resolution_clock::now();


    while (current_time < t_end) {
        SPDLOG_TRACE(std::to_string(current_time));

        SPDLOG_TRACE("Doing a Iteration with CellCalculator");
        //this applies Ghost Particle reflective boundary conditions, only 
        //if any of the boundaries is boundary_condition::ghost_reflective

        calculator.calcX(container);
        container.setup();
        calculator.calcF(container);
        calculator.calcV(container);

        if (current_time <= special_t) {
            calculator.applySpecialForces();
        }

        //std::cout << "one iteration done\n";


        iteration++;

        if (iteration % write_frequency == 0 && !performance_measurement) {
            container.cleanup();
            container.setup();
            vtkWriter.write(container, "out", iteration);
        }

        current_time += delta_t;




        /// loading bar
        static int loading_factor = std::max(write_frequency * 5.0, std::ceil(t_end / (delta_t * 100)));

        if (iteration % loading_factor == 0 && !performance_measurement) {
            barWidth = 50;
            pos = static_cast<size_t>(barWidth * (current_time / t_end));
            progressBar = "[" + std::string(pos, '=') + '>'
                          + std::string(barWidth - pos, ' ') + "] "
                          + std::to_string(int((current_time / t_end) * 100.0)) + "%\r";

            spdlog::info(progressBar);
        }


        if (performance_measurement) {
            perf_time_end = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> perf_duration = perf_time_end - perf_time_start;
            std::cout << "The Computation took: " << perf_duration.count() << " seconds" << std::endl;
        }

        if (!performance_measurement) {
            spdlog::info("[" + std::string(pos, '=') + ">] 100%\r");

            SPDLOG_INFO("output written. Terminating...\r");
        }


    }

}