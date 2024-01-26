#include <gtest/gtest.h>
#include "particleModel/updating/ThermoStats.h"
#include "particleModel/storage/CellContainer.h"

#include "inputHandling/FileReaderProgramArgs.h"


TEST(test_Rdf,test_basic){
    CellContainer container(15,15,0,3.0,3.0);
    ThermoStats thermoStats(container,0.1,0,0);

    // basic equal-sided pyramid (every side has length 2 * sqrt(2) ~ 2.828 )
    container.addParticle({2,2,2},{0,0,0},3);
    container.addParticle({2,0,0},{0,0,0},4);
    container.addParticle({0,2,0},{0,0,0},4);
    container.addParticle({0,0,2},{0,0,0},7);

    std::vector<double> rdf = thermoStats.getRadialDistributionFunction(0.1);
    //-> there should be 6 entries be counted in the 28th "backet"
    //because there are 6 pairs with the distance 2.828 -> they belong in 2.8 - 2.9 bucket

    double expected = 6 / ( (4.0 * M_PI / 3.0) * 
                      ( std::pow(2.9,3) - std::pow(2.8,3)) );
    
    std::cout << ArrayUtils::to_string(rdf) << "\n";

    for(size_t i = 0; i < rdf.size();i++){
        if(i == 28){ // 28 * 0.1(step size) = 2.8 -> bucket 2.8 - 2.9
            ASSERT_EQ(expected,rdf[i]);
        }else{
            ASSERT_EQ(0,rdf[i]);
        }
    }

}



