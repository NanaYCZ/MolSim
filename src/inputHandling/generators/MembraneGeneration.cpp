#include "MembraneGeneration.h"
#include "inputHandling/FileReaderProgramArgs.h"

void generateMembrane(FileReader::MembraneData& membrane, FileReader::SpecialForcesData& specialForces,CellContainer& container) {
    for(int z = 0; z < membrane.N3; z++) {

        for (int y = 0; y < membrane.N2; y++) {

            for (int x = 0; x < membrane.N1; x++) {
                std::array<double, 3> special = {0,0,0};
                std::array<double, 3> cords(membrane.x);
                std::array<double, 3> vel(membrane.v);
                std::array<int, 3> grid({x,y,z});

                cords[0] += x * membrane.h;
                cords[1] += y * membrane.h;
                cords[2] += z * membrane.h;

                if (int(specialForces.position[0])== x && int(specialForces.position[1])==y && int(specialForces.position[2])==z){
                    special=specialForces.f;
                }

                container.addParticle(cords, vel, grid,membrane.a, membrane.h, membrane.m, {0,0,7});

            }
        }
    }
}


void addMembranes(CellContainer &container, std::list<FileReader::SpecialForcesData> specialForces, std::list<FileReader::MembraneData> membranes) {
    for (auto &membrane : membranes) {
        for (auto &specialForce : specialForces){
            generateMembrane(membrane, specialForce,container);
        }

    }
}

