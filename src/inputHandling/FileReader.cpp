
#define SPDLOG_ACTIVE_LEVEL SPDLOG_LEVEL_TRACE

#include "FileReader.h"
#include "FileReaderProgramArgs.h"
#include "xmlParsing/parameters.hpp"

#include <spdlog/spdlog.h>
#include <array>
#include <fstream>
#include <iostream>
#include <list>
#include <sstream>
#include <stdexcept>
#include <type_traits>

FileReader::FileReader() = default;

FileReader::~FileReader() = default;





/**
 * @brief Sets a boundary condition based on the specified string.
 *
 * This function sets the boundary condition according to the provided string value.
 * If the specified string matches "reflective" or "outflow", it assigns the corresponding enum value
 * to the 'boundary' parameter.
 *
 * @param boundary Reference to a FileReader::boundary_conditions enum variable where the boundary condition will be set.
 * @param specified_cond A string specifying the desired boundary condition ("reflective" or "outflow"), which normally comes from the xml file.
 *
 * @throws std::invalid_argument If the string is not recognized(not "reflective" or "outflow"), an invalid_argument exception is thrown,
 *                               providing an error message indicating the incorrect string that was given.)
 */
void set_boundary_conditional(boundary_conditions& boundary,std::string specified_cond){
    if(specified_cond == "reflective")
        boundary = boundary_conditions::reflective;
    else if(specified_cond == "outflow")
        boundary = boundary_conditions::outflow;
    else if(specified_cond == "periodic")
        boundary = boundary_conditions::periodic;
    else if(specified_cond == "ghost_reflective")
        boundary = boundary_conditions::ghost_reflective;
    else
        throw std::invalid_argument("The Boundary Conditions were not correctly specified, you gave: " + specified_cond);

}

/**
 * @brief reads xml file and constructs ProgramsArgs struct corresponding to xml file
 * 
 * 
 * Reads an XML file 'filename' and uses Codesynthesis to parse/ validate the xml file then.
 * The information of the object returned from the XML file parser is then writen into an ProgramArgs struct
 * 
 * 
 * @param filename XML file according to parameters.xsd (can be found in input/ folder)
 * 
 * @returns a ProgramArgs struct with the information from the file
 * 
*/
FileReader::ProgramArgs FileReader::readProgramArguments(std::string filename){


    std::unique_ptr<parameters> params = parameters_(filename);

    auto out_params = params->outputParameters();
    auto sim_params = params->simulationParameters();
    auto cuboids = params->cuboids();
    auto spheres = params->spheres();
    auto rdf_params = sim_params.Rdf();
    auto boundary_conditions_xml = sim_params.boundaryConditions();



    boundary_conditions positive_z, negative_z ,positive_x, negative_x , positive_y , negative_y;

    set_boundary_conditional(positive_z,boundary_conditions_xml.boundaryConditionsPositiveZ());
    set_boundary_conditional(negative_z,boundary_conditions_xml.boundaryConditionsNegativeZ());
    set_boundary_conditional(positive_x,boundary_conditions_xml.boundaryConditionsPositiveX());
    set_boundary_conditional(negative_x,boundary_conditions_xml.boundaryConditionsNegativeX());
    set_boundary_conditional(positive_y,boundary_conditions_xml.boundaryConditionsPositiveY());
    set_boundary_conditional(negative_y,boundary_conditions_xml.boundaryConditionsNegativeY());



    ProgramArgs args;
    args.delta_t = sim_params.deltaT();
    args.t_end = sim_params.tEnd();
    args.cut_off_radius = sim_params.cutOffRadius();
    args.cell_size = sim_params.cellSize();
    args.gravity_factor = sim_params.gravityFactor().present() ? sim_params.gravityFactor().get() : 0;
    args.force_type = sim_params.forceType();

    if(sim_params.diffusionStatFrequency().present())
        args.diff_frequency =  sim_params.diffusionStatFrequency().get();

    if(sim_params.Rdf().present()){
        auto rdf = sim_params.Rdf().get();
        args.rdf_interval_and_frequency = std::pair(rdf.rdfTimeInterval(),rdf.rdfStatFrequency());
    }
    if(sim_params.Thermostats().present()){
        auto thermo = sim_params.Thermostats().get();

        if(thermo.maxTempDiff().present())
            args.max_temp_diff = std::make_optional<double>(thermo.maxTempDiff().get());

        if(thermo.targetTemp().present())
            args.target_temp = std::make_optional<double>(thermo.targetTemp().get());

        args.init_temp = thermo.initTemp();
        args.thermo_stat_frequency = thermo.thermoStatFrequency();

        args.calculate_thermostats = true;
    }//if not Thermostats present, struct will have default dummy values

    args.boundaries = {positive_z,negative_z,positive_x,negative_x,positive_y,negative_y};
    
    args.domain_dimensions = {sim_params.domainDimensions().x(),sim_params.domainDimensions().y(),sim_params.domainDimensions().z()};

    args.file_basename = out_params.baseName();
    args.write_frequency = out_params.writeFrequency();

    args.checkpoint_input_file = out_params.checkpointInputFileName().present() ? 
                                  std::optional<std::string>(out_params.checkpointInputFileName().get())
                                : std::nullopt;

    args.checkpoint_output_file = out_params.checkpointOutputFileName().present() ?
                                  std::optional<std::string>(out_params.checkpointOutputFileName().get())
                                : std::nullopt;
    
    for(size_t i = 0; i < cuboids.size() ; i++){
        CuboidData c;
        auto cuboid = cuboids[i];
        c.x = { cuboid.position().x(), cuboid.position().y(), cuboid.position().z() };
        c.v = { cuboid.velocity().x(), cuboid.velocity().y(), cuboid.velocity().z() };
    

        c.N1 = cuboid.dimensions().x();
        c.N2 = cuboid.dimensions().y();
        c.N3 = cuboid.dimensions().z();

        c.m = cuboid.mass();
        c.h = cuboid.meshWidth();
        c.sigma = cuboid.sigma();
        c.epsilon = cuboid.epsilon();

        //zero by default
        c.avg_v =  cuboid.meanVelocity().present() ? std::optional<double>(cuboid.meanVelocity().get()) : std::nullopt;

        args.cuboids.push_back(c);
    }


    for(size_t i = 0; i < spheres.size() ; i++){
        SphereData s;
        auto sphere = spheres[i];
        s.CenterPosition = { sphere.center_position().x(), sphere.center_position().y(), sphere.center_position().z() };
        s.Velocity = { sphere.velocity().x(), sphere.velocity().y(), sphere.velocity().z() };
        s.mass = sphere.mass();
        s.radius = sphere.radius();
        s.meshWidth = sphere.meshWidth();
        s.sigma = sphere.sigma();
        s.epsilon = sphere.epsilon();

        //zero by default
        s.avg_v = sphere.meanVelocity().present() ? std::optional<double>(sphere.meanVelocity().get()) : std::nullopt;
        

        args.spheres.push_back(s);
    }


    return args;    
}


void FileReader::initializeCorrectInitialTemp(FileReader::ProgramArgs& args){
    bool initial_temp_zero = true; 
        for(FileReader::CuboidData& cuboid : args.cuboids){
            if( ! (cuboid.v[0] == 0 && cuboid.v[1] == 0 && cuboid.v[2] == 0)){
                std::string msg = "don't initalize Temp, because of velocity: " +
                            std::to_string(cuboid.v[0]) + " , " + std::to_string(cuboid.v[1])
                            + " , " +  std::to_string(cuboid.v[2]) + "\n";
                SPDLOG_INFO(msg);
                initial_temp_zero = false;
                break;
            }
        }
        for(FileReader::SphereData& sphere : args.spheres){
            if( ! (sphere.Velocity[0] == 0 && sphere.Velocity[1] == 0 && sphere.Velocity[2] == 0)){
                std::string msg = "don't initalize Temp, because of velocity: " +
                            std::to_string(sphere.Velocity[0]) + " , " + std::to_string(sphere.Velocity[1])
                            + " , " +  std::to_string(sphere.Velocity[2]) + "\n";
                SPDLOG_INFO(msg);
                initial_temp_zero = false;
                break;
            }
        }
        if(initial_temp_zero){
            //initalize Particles with Maxwell-Boltzmann to right temperature
            for(FileReader::CuboidData& cuboid : args.cuboids){
                cuboid.avg_v = sqrt(args.init_temp/cuboid.m);
            }
            for(FileReader::SphereData& sphere : args.spheres){
                sphere.avg_v = sqrt(args.init_temp/sphere.mass);
            }
        }else{
            SPDLOG_INFO("Didn't apply inital Temperature of Thermostats,\n because not all initial velocities were zero");
        }
}











