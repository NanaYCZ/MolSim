#define SPDLOG_ACTIVE_LEVEL SPDLOG_LEVEL_TRACE

#include "inputHandling/generators/CuboidGeneration.h"
#include "inputHandling/generators/SphereGeneration.h"
#include "inputHandling/generators/MembraneGeneration.h"
#include "inputHandling/FileReaderProgramArgs.h"
#include "inputHandling/FileReader.h"
#include "Simulation.h"
#include <spdlog/sinks/basic_file_sink.h>
#include <spdlog/spdlog.h>
#include <unistd.h>
#include <optional>
#include <omp.h>




/**
 * @brief Parsing the program arguments and starting the Simulation
*/
int main(int argc, char *argsv[])
{
    //initialize default values
    bool performance_measurement = false;
    bool isMembrane = true;
    spdlog::level::level_enum logging_level = spdlog::level::info;

    std::string filename;
    FileReader fileReader;

    auto msg = "Usage ./MolSim [-l<String>] [-p] [-o] -f<String>\n"
               " Info:              See the /input folder for the parameters.xsd schema, in which program arguments should be specified\n"
               " -f<String>:        gives the filename of an .xml file, that has to follow\n"
               "                    the xsd schema defined in input/parameters.xsd.\n"
               "                    from this file all programm arguments / options will be read(see README)\n"
               " -l<String>:        specifies the level of logging, e.g. how fine grained programm logs are.\n"
               "                    can either be \"off\" \"trace\", \"debug\", \"info\", \"error\" or \"critical\".\n"
               "                    The default level is \"debug\".\n"
               " -h                 prints a help message\n"
               " -p                 if the flag is set, the programm will measure the time for the execution.\n"
               "                    therefore no vtk output and no logging will happen (specifing a log level at\n"
               "                    the same time is undefined behaviour)\n"
               "\n"
               "Returns:\n"
               "                  several .vtu files that can be used for visualisation in Paraview\n";

    std::cout <<  "Hello from MolSim for PSE!" << std::endl;

    ///variables for the argument parsing
    int opt;
    std::string log_mode;
    while ((opt = getopt(argc, argsv, "f:l:hpo")) != -1)
    {
        switch (opt)
        {
            case 'l':
                log_mode = std::string(optarg);
                if(log_mode=="off"){
                    logging_level = spdlog::level::off;
                }else if(log_mode=="trace"){
                    logging_level = spdlog::level::trace;   //enables current time logging
                }else if(log_mode=="debug"){
                    logging_level = spdlog::level::debug;   //enables logs.txt writing
                }else if(log_mode=="info"){
                    logging_level = spdlog::level::info;    //enables progress logging
                }else if(log_mode=="error"){
                    logging_level = spdlog::level::err;
                }else if(log_mode=="critical"){
                    logging_level = spdlog::level::critical;
                }
                break;
            case 'f':
                filename = std::string(optarg);
                break;
            case 'p':
                performance_measurement = true;
                logging_level = spdlog::level::off;
                break;
            case 'h':
                std::cout << msg;
                return 0;
            default:
                std::cerr << msg;
                return 1;
        }
    }

    auto logger = spdlog::basic_logger_mt("logger", "logs.txt");
    spdlog::set_level(logging_level);


    FileReader::ProgramArgs args = fileReader.readProgramArguments(filename);
    SPDLOG_INFO("Read:\n" + args.to_string());

    //check if initial velocities need to be initialized according to initialTemp
    // modifies the cuboid / spheres structs such that their mean velocities
    // fit the Temperature and the Maxwell-Boltzmann-Distribution is applied
    // with the correct arguments
    if(args.calculate_thermostats){
       FileReader::initializeCorrectInitialTemp(args);
    }

    if(args.choose_amount_threads.has_value()){
        int number_of_threads = args.choose_amount_threads.value();
        if(number_of_threads < 1){
            throw std::invalid_argument("It is not possible to simulate with the number"
                                        " of threads you provided: " + std::to_string(number_of_threads));
        }
        omp_set_num_threads(number_of_threads);
    }



    auto container = CellContainer(
            args.domain_dimensions[0],
            args.domain_dimensions[1],
            args.domain_dimensions[2],
            args.cut_off_radius,
            args.boundaries,
            std::array<double, 3>{0.0, 0.0, args.gravity_factor},
            args.parallelization_version
    );
    addMembranes(container,args.membranes);
    addCuboids(container,args.cuboids);
    addSpheres(container,args.spheres,3);
    auto io = outputWriter::VTKWriter();
    auto calculator=calculator::CellCalculator(1.0, 1.0, args.cut_off_radius);
    calculator.setDim(3);
    calculator.setMembrane(isMembrane);
    calculator.setRZero(2.2);
    calculator.setDeltaT(args.delta_t);
    calculator.setMembraneForceParameter(300);


    runSimulation(container,calculator,io,args.parallelization_version,args.t_end,args.delta_t,150,args.write_frequency,
                  performance_measurement);


    }

}