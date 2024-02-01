#include <gtest/gtest.h>
#include "inputHandling/FileReader.h"
#include "inputHandling/FileReaderProgramArgs.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <array>


std::string to_string(boundary_conditions cond){
    if(cond == boundary_conditions::reflective)
        return "reflective";
    else if(cond == boundary_conditions::periodic)
        return "periodic";
    else if(cond == boundary_conditions::outflow)
        return "outflow";
    else    
        throw std::invalid_argument("Not a correct boundary condition");
}

std::string to_string(force_type force,double r_l){
    if(force == force_type::smoothedLJ)
        return "<smoothedLJ>\n"
               "        <r_l>" + std::to_string(r_l) + "</r_l>\n"
               "</smoothedLJ>";
    else if(force == force_type::LJ)
        return "<LJ/>";
    else if(force == force_type::gravitational)
        return "<gravitational/>";
    else
        throw std::invalid_argument("Not a correct force type");
}


std::string to_string(concurrency_strategy strategy,std::optional<int> numThreads){
    if(strategy == concurrency_strategy::first_method){
        std::string string_return = "<first_method>\n";
        if(numThreads)
                   string_return += "                <numThreads>" +  std::to_string(numThreads.value()) + "</numThreads>\n";
        string_return += "</first_method>";
        return  string_return;
    }
    else if(strategy == concurrency_strategy::second_method){
        std::string string_return = "<second_method>\n";
        if(numThreads)
            string_return += "                <numThreads>" +  std::to_string(numThreads.value()) + "</numThreads>\n";
        string_return += "</second_method>";
        return  string_return;

    }else if(strategy == concurrency_strategy::serial)
        return "<serial/>";
    else
        throw std::invalid_argument("Not a correct concurrency strategy");
}


/**
 * @brief writes the struct programArgs into the file 'filename' in xml format
 * 
 * 
 * The function is a helper for the Test afterwards and writes the information
 * from the programArgs struct into an xml file according to a xsd schema
 * defined in parameters.xsd (is copied into build and test folder automatically)
 * 
 * @param programArgs struct containing info for the test
 * @param filename name of the file, in which the info is written
 * 
*/
void writeArgsIntoXMLFile(const FileReader::ProgramArgs &programArgs,std::string filename){
    // Construct the XML string
    std::ostringstream xmlStream;
    xmlStream << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>" << std::endl;
    xmlStream << "<parameters xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" "
            "xsi:noNamespaceSchemaLocation=\"parameters.xsd\">" << std::endl;

    // Add outputParameters node under parameters
    xmlStream << "  <outputParameters>" << std::endl;
    xmlStream << "    <baseName>" << programArgs.file_basename << "</baseName>" << std::endl;
    xmlStream << "    <writeFrequency>" << programArgs.write_frequency << "</writeFrequency>" << std::endl;
    if(programArgs.checkpoint_input_file.has_value())
        xmlStream << "    <checkpointInputFileName>" << programArgs.checkpoint_input_file.value() << "</checkpointInputFileName>" << std::endl;

    if(programArgs.checkpoint_output_file.has_value())
        xmlStream << "    <checkpointOutputFileName>" << programArgs.checkpoint_output_file.value() << "</checkpointOutputFileName>" << std::endl;

    xmlStream << "  </outputParameters>" << std::endl;

    // Add simulationParameters node under parameters
    xmlStream << "  <simulationParameters>" << std::endl;
    xmlStream << "    <tEnd>" << programArgs.t_end << "</tEnd>" << std::endl;
    xmlStream << "    <deltaT>" << programArgs.delta_t << "</deltaT>" << std::endl;
    xmlStream << "    <cutOffRadius>" << programArgs.cut_off_radius << "</cutOffRadius>" << std::endl;
    xmlStream << "    <cellSize>" << programArgs.cell_size << "</cellSize>" << std::endl;
    xmlStream << "    <gravityFactor>" << programArgs.gravity_factor << "</gravityFactor>" << std::endl;
    xmlStream << "     <forceType>" << to_string(programArgs.force_type_param,programArgs.r_l.value_or(1.9)) << "</forceType> " << std::endl;
    if(programArgs.parallelization_version){
        xmlStream << "      <parallelizationVersion>" <<
        to_string(programArgs.parallelization_version.value(),programArgs.choose_amount_threads)
        << "</parallelizationVersion> " << std::endl;
    }
    if(programArgs.rdf_interval_and_frequency){
        xmlStream << "      <Rdf>" <<
                        "<rdfIntervalSize>" + std::to_string(programArgs.rdf_interval_and_frequency->first) + "</rdfIntervalSize>\n"
                        +"<rdfStatFrequency> " + std::to_string(programArgs.rdf_interval_and_frequency->second) +  "  </rdfStatFrequency>\n"
                  << "</Rdf> " << std::endl;
    }
    if(programArgs.diff_frequency){
        xmlStream << "    <diffusionStatFrequency>" << programArgs.diff_frequency.value() << "</diffusionStatFrequency>" << std::endl;
    }
    if(programArgs.calculate_thermostats){
    xmlStream << "    <Thermostats>" << std::endl;
    xmlStream << "      <initTemp>" << programArgs.init_temp << "</initTemp>" << std::endl;                    
    if(programArgs.target_temp.has_value())
        xmlStream << "      <targetTemp>" << programArgs.target_temp.value() << "</targetTemp>" << std::endl;                    

    xmlStream << "      <thermoStatFrequency>" << programArgs.thermo_stat_frequency << "</thermoStatFrequency>" << std::endl;    
    if(programArgs.max_temp_diff.has_value())
        xmlStream << "      <maxTempDiff>" << programArgs.max_temp_diff.value() << "</maxTempDiff>" << std::endl;                    

    xmlStream << "</Thermostats>" << std::endl;
    }
    xmlStream << "    <boundaryConditions>" << std::endl;
    xmlStream << "      <boundaryConditionsPositiveZ>" << (programArgs.boundaries[0] == boundary_conditions::reflective ?  "reflective" : "outflow") << "</boundaryConditionsPositiveZ> " << std::endl;
    xmlStream << "      <boundaryConditionsNegativeZ>" << (programArgs.boundaries[1] == boundary_conditions::reflective ?  "reflective" : "outflow") << "</boundaryConditionsNegativeZ> " << std::endl;
    xmlStream << "      <boundaryConditionsPositiveX>" << (programArgs.boundaries[2] == boundary_conditions::reflective ?  "reflective" : "outflow") << "</boundaryConditionsPositiveX> " << std::endl; 
    xmlStream << "      <boundaryConditionsNegativeX>" << (programArgs.boundaries[3] == boundary_conditions::reflective ?  "reflective" : "outflow") << "</boundaryConditionsNegativeX> " << std::endl; 
    xmlStream << "      <boundaryConditionsPositiveY>" << (programArgs.boundaries[4] == boundary_conditions::reflective ?  "reflective" : "outflow") << "</boundaryConditionsPositiveY> " << std::endl;
    xmlStream << "      <boundaryConditionsNegativeY>" << (programArgs.boundaries[5] == boundary_conditions::reflective ?  "reflective" : "outflow") << "</boundaryConditionsNegativeY> " << std::endl; 
    xmlStream << "</boundaryConditions>" << std::endl;
    xmlStream << "    <domainDimensions>" <<  std::endl;
    xmlStream << "      <x>" << programArgs.domain_dimensions[0] << "</x>" << std::endl;
    xmlStream << "      <y>" << programArgs.domain_dimensions[1] << "</y>" << std::endl;
    xmlStream << "      <z>" << programArgs.domain_dimensions[2] << "</z>" << std::endl;
    xmlStream << "    </domainDimensions>" << std::endl;
    xmlStream << "  </simulationParameters>" << std::endl;

    // Add cuboids under parameters
    for (const auto &cuboid : programArgs.cuboids) {
        xmlStream << "  <cuboids>" << std::endl;
        xmlStream << "    <position>" << std::endl;
        xmlStream << "      <x>" << cuboid.x[0] << "</x>" << std::endl;
        xmlStream << "      <y>" << cuboid.x[1] << "</y>" << std::endl;
        xmlStream << "      <z>" << cuboid.x[2] << "</z>" << std::endl;
        xmlStream << "    </position>" << std::endl;
        xmlStream << "    <velocity>" << std::endl;
        xmlStream << "      <x>" << cuboid.v[0] << "</x>" << std::endl;
        xmlStream << "      <y>" << cuboid.v[1] << "</y>" << std::endl;
        xmlStream << "      <z>" << cuboid.v[2] << "</z>" << std::endl;
        xmlStream << "    </velocity>" << std::endl;
        xmlStream << "    <dimensions>" << std::endl;
        xmlStream << "      <x>" << cuboid.N1 << "</x>" << std::endl;
        xmlStream << "      <y>" << cuboid.N2 << "</y>" << std::endl;
        xmlStream << "      <z>" << cuboid.N3 << "</z>" << std::endl;
        xmlStream << "    </dimensions>" << std::endl;
        if(cuboid.avg_v.has_value())
            xmlStream << "      <meanVelocity>" << cuboid.avg_v.value() << "</meanVelocity>" << std::endl;
        xmlStream << "    <mass>" << cuboid.m << "</mass>" << std::endl;
        xmlStream << "    <meshWidth>" << cuboid.h << "</meshWidth>" << std::endl;
        xmlStream << "    <sigma>" << cuboid.sigma << "</sigma>" << std::endl;
        xmlStream << "    <epsilon>" << cuboid.epsilon << "</epsilon>" << std::endl;
        xmlStream << "  </cuboids>" << std::endl;
    }

    for (const auto &sphere : programArgs.spheres) {
        xmlStream << "  <spheres>" << std::endl;
        xmlStream << "    <center_position>" << std::endl;
        xmlStream << "      <x>" << sphere.CenterPosition[0] << "</x>" << std::endl;
        xmlStream << "      <y>" << sphere.CenterPosition[1] << "</y>" << std::endl;
        xmlStream << "      <z>" << sphere.CenterPosition[2] << "</z>" << std::endl;
        xmlStream << "    </center_position>" << std::endl;
        xmlStream << "    <velocity>" << std::endl;
        xmlStream << "      <x>" << sphere.Velocity[0] << "</x>" << std::endl;
        xmlStream << "      <y>" << sphere.Velocity[1] << "</y>" << std::endl;
        xmlStream << "      <z>" << sphere.Velocity[2] << "</z>" << std::endl;
        xmlStream << "    </velocity>" << std::endl;
        if(sphere.avg_v.has_value())
            xmlStream << "      <meanVelocity>" << sphere.avg_v.value() << "</meanVelocity>" << std::endl;
        xmlStream << "    <mass>" << sphere.mass << "</mass>" << std::endl;
        xmlStream << "    <radius>" << sphere.radius << "</radius>" << std::endl;
        xmlStream << "    <meshWidth>" << sphere.meshWidth << "</meshWidth>" << std::endl;
        xmlStream << "    <sigma>" << sphere.sigma << "</sigma>" << std::endl;
        xmlStream << "    <epsilon>" << sphere.epsilon << "</epsilon>" << std::endl;
        xmlStream << "  </spheres>" << std::endl;
    }

    xmlStream << "</parameters>" << std::endl;

    // Write the XML string to a file
    std::ofstream outputFile(filename);
    if (outputFile.is_open()) {
        outputFile << xmlStream.str();
        outputFile.close();
        std::cout << "XML file generated successfully." << std::endl;
    } else {
        std::cerr << "Unable to open file for writing." << std::endl;
    }

}
    




/**
 * @brief Creates a dummy programArgs struct and writes it into a file "args_test.xml".
 *        Then the method from our project is used to again extract the information from the 
 *        xml file (the corresponding .xsd file has to be within the test folder as well,
 *        parameters.xsd should be copied into the /tests folder within the build folder 
 *        automatically). In the end the pre-defined struct and the struct gained 
 *        through writing + reading are compared for equality.
*/
TEST(test_readProgArgs,test_big){
    std::string filename = "args_test.xml";

    /*
    std::array<double, 3> CenterPosition;
    std::array<double, 3> Velocity;
    double mass;
    double radius;
    double meshWidth;
    double sigma;
    double epsilon;

    //by default no maxwell-boltzmann is applied
    std::optional<double> avg_v;
    */
    FileReader::SphereData sphere = {
        {1.5, 2.0, 3.0}, 
        {5.12, 3.343, 7.8}, 
        17.0,            
        2.5,             
        1.65,             
        0.5,             
        0.8,
        0.1          
    };

    /*
    std::array<double, 3> x, v;



    /// N1: amount of particles along dimension 1
    /// N2: amount of particles along dimension 2
    /// N3: amount of particles along dimension 3
    uint64_t N1, N2, N3;

    /// Mass m of the particles in the cuboid
    /// Mesh width h
    double m, h;

    /// sigma and epsilon parameters for the force calculation
    /// between particles of this cuboid
    double sigma, epsilon;

    /// Average velocity default 0 means by default no Maxwell-boltzmann is applied
    std::optional<double> avg_v;
    */

    // Initialize CuboidData with specific values
    FileReader::CuboidData cuboid = {
        {0.0, 23.2324, 9.9}, 
        {15.23, 1.435, 7.7}, 
        5,      //N1            
        9,      //N2     
        5,      //N3     
        20.0,   //m        
        2.0,    //h      
        0.3,             
        0.6,             
        0.1     //avg_v         
    };
    
    /*

    bool calculate_thermostats = false;

    double delta_t;
    double t_end;
    double cut_off_radius;
    double cell_size;
    double gravity_factor = 0;
    double init_temp = 0;
    std::optional<double> max_temp_diff = std::nullopt;
    std::optional<double> target_temp = std::nullopt;
    int thermo_stat_frequency = 0;
    std::array<boundary_conditions,6> boundaries;
    std::array<double,3> domain_dimensions;

    std::optional<std::string> checkpoint_input_file;
    std::optional<std::string> checkpoint_output_file;



    std::string file_basename = "out";
    size_t write_frequency = 10;

    */


    FileReader::ProgramArgs programArgs = {
        false,
        0.054,           // delta_t
        1422.0,          // t_end
        2.0,            //cut off radius
        1.0,            //cell size
        -12.44,     //gravity factor
        10.0,       //initial temp
        force_type::LJ,
        1.9,
        concurrency_strategy::serial,
        16,
        20.0,       //max temp diff
        std::pair<double,int>(1.0,100),
        0.5,
        30.0,       //target temp
        50,          //thermostat write frequency
        {boundary_conditions::reflective,boundary_conditions::reflective,
        boundary_conditions::reflective,boundary_conditions::reflective,
        boundary_conditions::reflective,boundary_conditions::reflective
        },    //boundary conditions
        {5,5,5},        //domain size
        "in_checkp",    //in checkpoint file
        "out_checkp",   //out checkpoint file
        "out",          // file_basename
        10,             // write_frequency
        {cuboid},      // spheres
        {sphere}       // cuboids
    };

    programArgs.calculate_thermostats = true;


    FileReader fileReader;

    // write the information from the struct into
    // a xml file according to the defined xml format
    writeArgsIntoXMLFile(programArgs,filename);

    // let the method from the fileReader
    // read out the data from the previously created file
    FileReader::ProgramArgs programArgs_read = fileReader.readProgramArguments(filename);


    std::cout << programArgs.to_string() << std::endl;
    std::cout << "Read:\n" ;
    std::cout << programArgs_read.to_string() << std::endl;

    // check if the structs are equal
    // comparison operation is overloaded for the defined structs
    ASSERT_EQ(programArgs_read,programArgs);
}




