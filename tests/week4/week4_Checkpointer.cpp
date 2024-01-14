#include <gtest/gtest.h>
#include "inputHandling/FileReader.h"
#include "inputHandling/generators/CuboidGeneration.h"
#include "inputHandling/Checkpointer.h"
#include "inputHandling/FileReaderProgramArgs.h"





/**
 * @brief simple test, that writes all particles of a CellContainer into a checkpoint file,
 *        then read them out again and compare if exactly the particles were
 *        read, that were written before
 * 
*/
TEST(test_Checkpointing,test_basic){
    FileReader::CuboidData cuboid = {
        {20, 20, 20}, 
        {0, 0, 0}, 
        10,               
        10,               
        10,               
        1,            
        1,             
        0.3,             
        0.6,             
        5              
    };

    CellContainer container(100,100,100,5.0,5.0);
    
    addCuboids(container,{cuboid});

    container.createPointers();

    std::list<Particle>& particles = container.getInstances();

    

    std::list<Particle> particles_list(particles.begin(),particles.end());



    Checkpointer::writeCheckpoint(particles_list,"checkpoint_test.txt");

    std::list<std::tuple<Particle,double,double>> particles_read;

    Checkpointer::readCheckpoint(particles_read,"checkpoint_test.txt");


    std::cout << "Checking list sizes" << std::endl;
    ASSERT_TRUE(particles_read.size() == particles_list.size());

    auto it1 = particles_list.begin();
    auto it2 = particles_read.begin();

    for(;it1 != particles_list.end() && it2 != particles_read.end();it1++ , it2++){
        Particle other_particle = std::get<0>(*it2);
        ASSERT_TRUE(*it1 == other_particle);
    }

}



/**
 * @brief advanced test, that writes all particles of a CellContainer into a checkpoint file,
 *        where the particles have different sigmas and epsilons, then read them out again and
 *        compare if exactly the particles were read, that were written before
 * 
*/
TEST(test_Checkpointing,test_advanced){
    FileReader::CuboidData cuboid = {
        {20, 20, 20}, 
        {0, 0, 0}, 
        10,               
        10,               
        10,               
        1,            
        1,             
        0.3,             
        0.6,             
        5              
    };

    FileReader::CuboidData cuboid_different_sigma_epsilon = {
        {20, 20, 20}, 
        {0, 0, 0}, 
        10,               
        10,               
        10,               
        1,            
        1,             
        4,             
        5,             
        5              
    };

    CellContainer container(100,100,100,5.0,5.0);

    CellContainer container_for_reading(100,100,100,5.0,5.0);
    
    addCuboids(container,{cuboid,cuboid_different_sigma_epsilon});

    container.createPointers();

    Checkpointer::storeCheckpointparticles(container,"checkpoint_test.txt");


    Checkpointer::addCheckpointparticles(container_for_reading,"checkpoint_test.txt");

    container_for_reading.createPointers();

    std::list<Particle>& particles_written = container.getInstances();

    std::list<Particle>& particles_read = container_for_reading.getInstances();

    std::cout << "Checking list sizes" << std::endl;
    ASSERT_TRUE(particles_read.size() == particles_written.size());

    auto it1 = particles_written.begin();
    auto it2 = particles_read.begin();

    for(;it1 != particles_written.end() && it2 != particles_read.end();it1++ , it2++){
        ASSERT_TRUE(*it1 == *it1);
    }

}










