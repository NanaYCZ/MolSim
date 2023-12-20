#include <gtest/gtest.h>
#include "inputHandling/FileReader.h"
#include "inputHandling/generators/CuboidGeneration.h"
#include "inputHandling/Checkpointer.h"
#include "inputHandling/FileReaderProgramArgs.h"






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

    std::vector<Particle>& particles = container.getInstances();

    

    std::list<Particle> particles_list(particles.begin(),particles.end());



    Checkpointer::writeCheckpoint(particles_list,"checkpoint_test.txt");

    std::list<std::tuple<Particle,double,double>> particles_read;

    Checkpointer::readCheckpoint(particles_read,"checkpoint_test.txt");

    //particles_read.reverse();


    std::cout << "Checking list sizes" << std::endl;
    ASSERT_TRUE(particles_read.size() == particles_list.size());

    auto it1 = particles_list.begin();
    auto it2 = particles_read.begin();

    for(;it1 != particles_list.end() && it2 != particles_read.end();it1++ , it2++){
        Particle other_particle = std::get<0>(*it2);
        ASSERT_TRUE(*it1 == other_particle);
    }

}













