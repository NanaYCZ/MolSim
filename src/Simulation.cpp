#define SPDLOG_ACTIVE_LEVEL SPDLOG_LEVEL_TRACE


#include "particleModel/updating/CellCalculator.h"
#include "Simulation.h"
#include <chrono>
#include <spdlog/spdlog.h>

void runSimulation(SimulationContainer &particleContainer, Calculator& model, const double end_time,
                   const double delta_t, bool performance_measurement) {

    outputWriter::VTKWriter writer;
    auto logger = spdlog::get("logger");

    std::chrono::high_resolution_clock::time_point perf_time_start, perf_time_end;

    double current_time = 0;
    int iteration = 0;

    std::string progressBar;
    size_t barWidth, pos;

    // calculate inital force:
    model.initalize();
    SPDLOG_LOGGER_DEBUG(logger, "Particles in the simulation:");
    SPDLOG_LOGGER_DEBUG(logger, particleContainer.to_string());
    logger->flush();

    // for this loop, we assume: current x, current f and current v are known
    if (performance_measurement)
        perf_time_start = std::chrono::high_resolution_clock::now();
    while (current_time < end_time) {
        SPDLOG_TRACE(std::to_string(current_time));


        model.iteration();

        iteration++;

        if (iteration % 10 == 0 && !performance_measurement) {
            writer.initializeOutput(particleContainer.size());
            particleContainer.plotParticles(writer);
            writer.writeFile("out", iteration);
        }


        /// loading bar
        if (iteration % 50 == 0 && !performance_measurement) {
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

    spdlog::info("[" + std::string(pos, '=') + ">] 100%\r");
    SPDLOG_INFO("output written. Terminating...\r");
}



void plotParticles(ParticleContainer &particleContainer, int iteration) {
  std::string out_name("MD_vtk");
  outputWriter::XYZWriter writer;
  writer.plotParticles(particleContainer, out_name, iteration);
}