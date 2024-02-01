#include "MembraneGeneration.h"
#include "inputHandling/FileReaderProgramArgs.h"

void generateMembrane(FileReader::MembraneData& membraneData, CellContainer& container, size_t dim) {

    for(uint64_t z = 0; z < membraneData.N3; z++) {

        for (uint64_t y = 0; y < membraneData.N2; y++) {

            for (uint64_t x = 0; x < membraneData.N1; x++) {
                std::array<double, 3> cords(membraneData.x);
                std::array<double, 3> vel(membraneData.v);
                std::array<uint64_t, 3> grids={x, y, z};

                cords[0] += x * membraneData.h;
                cords[1] += y * membraneData.h;
                cords[2] += z * membraneData.h;

                container.addParticle(grids, cords, vel, membraneData.m);

            }
        }
    }
}
