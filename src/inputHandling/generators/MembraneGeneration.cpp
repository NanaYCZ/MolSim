#include "MembraneGeneration.h"
#include "inputHandling/FileReaderProgramArgs.h"

void generateMembrane(FileReader::MembraneData& membrane, CellContainer& container) {

    for(int z = 0; z < membrane.N3; z++) {

        for (int y = 0; y < membrane.N2; y++) {

            for (int x = 0; x < membrane.N1; x++) {
                std::array<double, 3> cords(membrane.x);
                std::array<double, 3> vel(membrane.v);
                std::array<int, 3> grid({x,y,z});

                cords[0] += x * membrane.h;
                cords[1] += y * membrane.h;
                cords[2] += z * membrane.h;


                container.addParticle(cords, vel, grid,membrane.m);

            }
        }
    }
}


void addMembranes(CellContainer &container, std::list<FileReader::MembraneData> membranes) {
    for (auto &membrane : membranes) {
        generateMembrane(membrane, container);
    }
}

