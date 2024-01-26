#include <gtest/gtest.h>
#include "particleModel/updating/ThermoStats.h"
#include "particleModel/storage/CellContainer.h"
#include "particleModel/updating/CellCalculator.h"

#include "inputHandling/FileReaderProgramArgs.h"


TEST(test_Misc,test_force){
    CellContainer container(15,15,0,3.0,3.0);
    CellCalculator cellCalculator(container,0.001,3.6,1.9,
                          {boundary_conditions::periodic,boundary_conditions::periodic,
                           boundary_conditions::periodic,boundary_conditions::periodic,
                           boundary_conditions::periodic,boundary_conditions::periodic},"LJ");
    Particle a{{4.2,1.8,1.8},{0,0,0},1};
    Particle b{{6.6,3,4.2},{0,0,0},1};
    std::array<double, 3> F = cellCalculator.force(a, b, {0, 0, 0});
    std::cout << ArrayUtils::to_string(F) << "\n";

}



