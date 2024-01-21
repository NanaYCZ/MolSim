#include <gtest/gtest.h>
#include "inputHandling/FileReader.h"
#include "inputHandling/FileReaderProgramArgs.h"
#include "inputHandling/generators/CuboidGeneration.h"
#include "inputHandling/generators/SphereGeneration.h"
#include "particleModel/storage/CellContainer.h"
#include "particleModel/storage/CellContainerIterators.h"
#include "particleModel/updating/CellCalculator.h"
#include "particleModel/updating/ThermoStats.h"
#include "utils/ArrayUtils.h"

/**
 * @brief calculates the temperature of the system(all particles in the container, that are 
 *        within the domain boundary)
 * 
 * @return temperature of the system
*/
double getTemp(CellContainer& container){

    double temp = 0;

    for(auto iter = container.begin(); iter != container.end(); ++iter){
        for(Particle* particle_ptr : *iter){   
            //std::cout << "encountered particle"  << std::endl;
            auto v = particle_ptr->getV();
            // if(ArrayUtils::L2Norm(v) > 50)
            //     std::cout << "high veloctiy for: " << particle_ptr->toString() << "\n" << std::endl;

            temp += ((v[0]* v[0] + v[1] * v[1] + v[2] * v[2] ) * particle_ptr->getM());
        }
    }

    std::cout << "Afterwards have kinetic energy(times 2): " << temp << std::endl;

    //dimension should be 2 and boltzman constant is 1
    temp = temp / ( 2.0 * container.size() *  1);

    return temp;
}

/**
 * @brief Tests if the Thermostat applied to a CellContainer with particles changes the 
 *        temperature of the system(particles in the CellContainer) correctly. First some 
 *        particles are added to a CellContainer. Then without anything else, the Thermostat
 *        is applied to the CellContainer by calling it on a CellCalculator, that wraps our
 *        CellContainer. The CellCalculator was initialized with  
 *        target_temp = 30 and no max_temp_diff. Therefore calling `applyThermostats` should
 *        directly change the temperature of the system to 30.
 *
 */
TEST(test_Thermo_Stat,test_basic){
    CellContainer container(15,15,0,3.0,3.0);
    CellCalculator calculator(container,0.0014,3.0,1.9,
        {boundary_conditions::outflow,boundary_conditions::outflow,
        boundary_conditions::outflow,boundary_conditions::outflow,
        boundary_conditions::outflow,boundary_conditions::outflow},"LJ");
    ThermoStats thermoStats(container,0.0014,30.0);

    container.addParticle({1,1,0},{2,2,2},3);
    container.addParticle({6,5,0},{2,2,2},4);
    container.addParticle({7,12,0},{2,2,2},4);
    container.addParticle({7,4,0},{3,4,5},7);

    container.createPointers();


    
    //after this call the Temperature of the system should be 30 (because it's the target temperature)
    thermoStats.applyThermostats();

    //particles afterwards

    double temp = getTemp(container);

    //due to rounding errors etc. we can't expect to get the exact double temperature again
    ASSERT_NEAR(30,temp,0.00001);

}

/**
 * @brief Tests if the Thermostat correctly heats up a system(particles in the CellContainer).
 *        Frist some particles are added to a CellContainer. The system will then have the inital
 *        temperature=40. Then a CellCalculator, that wraps the CellContainer is initialized 
 *        with target_temp = 100 and max_temp_diff = 5.0 and we do 20 iterations of the simulation
 *        and apply the the Thermostat in every step. Since we start with temp=40 and then can increase
 *        the temp by 5 every step, after 20 steps our system should have the target temp=100
 */
TEST(test_Thermo_Stat,test_heating){
    CellContainer container(50,50,0,3.0,3.0);
    CellCalculator calculator(container,0.0014,3.0,1.9,
        {boundary_conditions::reflective,boundary_conditions::reflective,
        boundary_conditions::reflective,boundary_conditions::reflective,
        boundary_conditions::reflective,boundary_conditions::reflective
        },"LJ");  
    ThermoStats thermoStats(container,0.0014,100.0,5.0);


    //have some particles to simulate
    container.addParticle({1,1,0},{2,2,0},3);
    container.addParticle({3,3,0},{2,2,0},4);
    container.addParticle({7,12,0},{2,2,0},4);
    container.addParticle({7,7,0},{3,4,0},7);
    container.addParticle({12,8,0},{2,4,0},3);
    container.addParticle({30,23,0},{3,4,0},8);
    container.addParticle({10,8,0},{1,1,0},9);

    container.createPointers();

    double temp = getTemp(container);

    std::cout << "The Temperature before the simulation is: " << temp << std::endl; 
    //This will be 40

    calculator.calculateF();
    calculator.shiftF();


    for(int i = 0; i < 20; i++){
        calculator.calculateX();
        calculator.calculateF();
        calculator.calculateV();
        calculator.shiftF();
        thermoStats.applyThermostats();
        temp = getTemp(container);
        std::cout << "The current Temperature is: " << temp << std::endl; 
    }

    temp = getTemp(container);

    std::cout << "The Temperature after the simulation is: " << temp << std::endl;
    //due to rounding errors etc. we can't expect to get the exact double temperature 
    ASSERT_NEAR(100,temp,0.00001);
}

/**
 * @brief Tests if the Thermostat correctly cools down a system(particles in the CellContainer).
 *        Frist some particles are added to a CellContainer. The system will then have the inital
 *        temperature=40. Then a CellCalculator, that wraps the CellContainer is initialized 
 *        with target_temp = 20 and max_temp_diff = 1.0 and we do 30 iterations of the simulation
 *        and apply the the Thermostat in every step. Since we start with temp=40 and then can decrease
 *        the temp by 1.0 every step, after 30 steps our system should have the target temp=20
 * 
 */
TEST(test_Thermo_Stat,test_cooling){
    CellContainer container(50,50,0,3.0,3.0);
    CellCalculator calculator(container,0.0014,3.0,1.9,
        {boundary_conditions::reflective,boundary_conditions::reflective,
        boundary_conditions::reflective,boundary_conditions::reflective,
        boundary_conditions::reflective,boundary_conditions::reflective
        },"LJ");  

    ThermoStats thermoStats(container,0.0014,20.0,1.0);
    //max_temp_diff is 5 and target_temp is 20
    //so in every Thermostat iteration, the temperature is increased by one maximum


    //have some particles to simulate
    container.addParticle({1,1,0},{2,2,0},3);
    container.addParticle({6,5,0},{2,2,0},4);
    container.addParticle({7,12,0},{2,2,0},4);
    container.addParticle({7,4,0},{3,4,0},7);
    container.addParticle({20,30,0},{2,4,0},3);
    container.addParticle({3,5,0},{3,4,0},8);
    container.addParticle({24,8,0},{1,1,0},9);

    container.createPointers();

    double temp = getTemp(container);

    std::cout << "The Temperature before the simulation is: " << temp << std::endl; 
    //This will be 40

    calculator.calculateF();
    calculator.shiftF();


    for(int i = 0; i < 20; i++){
        calculator.calculateX();
        calculator.calculateF();
        calculator.calculateV();
        calculator.shiftF();
        thermoStats.applyThermostats();
        temp = getTemp(container);
        std::cout << "The current Temperature is: " << temp << std::endl; 
    }

    temp = getTemp(container);

    std::cout << "The Temperature after the simulation is: " << temp << std::endl;
    //due to rounding errors etc. we can't expect to get the exact double temperature again
    ASSERT_NEAR(20,temp,0.00001);
}


/**
 * @brief Just an informal check to see if the inital temp is correctly applied an the system
 *        has roughly the expected value of temp=30. Due to the non determinism, we don't check
 *        that the system actually has the exact inital temp = 30.
 * 
 * 
*/
TEST(test_Thermo_Stat,test_initial_Temp){
     FileReader::SphereData sphere = {
        {100, 100, 100}, 
        {0, 0, 0}, 
        1.0,            
        50,             
        1.0,             
        1,             
        1              
    };

    // Initialize CuboidData with specific values
    FileReader::CuboidData cuboid = {
        {20, 20, 0}, 
        {0, 0, 0}, 
        50,               
        50,               
        0,               
        1,            
        1,             
        0.3,             
        0.6,             
        0              
    };
    
    FileReader::ProgramArgs args = {
        false,              
        0.054,           // delta_t
        50.0,          // t_end
        2.0,            //cut of radius
        2.0,            //cell size
        0.0,            //gravity factor
        30.0,       //initial temp
        "LJ",
        std::nullopt,
        std::nullopt,
        std::nullopt,       //max temp diff
        std::nullopt,       //target temp
        50,          //thermostat write frequency
        {boundary_conditions::reflective,boundary_conditions::reflective,
        boundary_conditions::reflective,boundary_conditions::reflective,
        boundary_conditions::reflective,boundary_conditions::reflective
        },    //boundary conditions
        {500,500,0},        //domain size
        std::nullopt,       
        std::nullopt,
        "out",          // file_basename
        10,             // write_frequency
        {cuboid},      // spheres
        {sphere}       // cuboids
    };

    std::cout << args.to_string() << std::endl;

    CellContainer cellContainer(args.domain_dimensions[0],args.domain_dimensions[1],args.domain_dimensions[2],args.cut_off_radius,args.cell_size);

    FileReader::initializeCorrectInitialTemp(args);
    
    std::cout << args.to_string() << std::endl;

    addCuboids(cellContainer,args.cuboids);
    addSpheres(cellContainer,args.spheres,2);

    cellContainer.createPointers();

    //std::cout << cellContainer.to_string() << std::endl;
    double temp = getTemp(cellContainer);


    // no comparison possible, because nondeterministic initialization
    // with boltzmann distribution
    std::cout << "temp is: " << temp << std::endl;
    std::cout << "temp should be roughly 30" << std::endl;

}









