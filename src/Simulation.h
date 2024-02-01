#pragma once

#include "particleModel/CellCalculator.h"
#include "outputWriter/VTKWriter.h"

/**
 * @brief start of the particle simulation
 *
 * repeatedly iterates over particles inside particleContainer until
 * end_time is reached with delta_t steps for each iteration
 *
 * \image latex runtime_nice.png "Runtime Measurements"
 *
 * @param container contains all particles to simulate
 * @param calculator is responsible for the particle calculations
 * @param end_time timespan to simulate
 * @param delta_t time step for each iteration
 * @param write_frequency specifies the frequency for the vtk output
 * @param performance_measurement bool to set the performance measuring of the simulation
 */

void runSimulation(CellContainer &container, calculator::CellCalculator& calculator, outputWriter::VTKWriter vtkWriter, CellContainer::concurrency_strategy strategy, double t_end,double delta_t,size_t write_frequency,bool performance_measurement);

/**
* @brief plot the particles to a xyz-file
*/
void plotParticles(CellContainer &container, int iteration);





