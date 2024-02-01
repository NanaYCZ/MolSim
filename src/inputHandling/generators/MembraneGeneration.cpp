#include "MembraneGeneration.h"
#include "inputHandling/FileReaderProgramArgs.h"

void generateMembrane(FileReader::MembraneData& membrane, ParticleContainer &particles) {
    int total_count = membrane.N1 * membrane.N2;
    particles.reserve(total_count);


    for (int y = 0; y < membrane.N2; y++) {

        for (int x = 0; x < membrane.N1; x++) {

            int z=1;



                std::array<double, 3> cords(membrane.x);
                std::array<double, 3> vel(membrane.v);
                cords[0] += x * membrane.h;
                cords[1] += y * membrane.h;
                cords[2] += 0;

                Particle part = Particle{cords, vel, membrane.m};
                auto index = std::array<int, 3>{x, y, z};
                part.setGridIndex(index);
                part.setMembrane(true);
                part.setGrav({0,0,-0.001});


                if( x==17&&y==24 || x==17&&y==25 || x==18&&y==24 || x==18&&y==25 ){
                    std::array<double, 3> FZ_UP = {0, 0, 0.8};
                    part.setBaseForce(FZ_UP);
                }
                part.applyBaseForceAndGrav();
                particles.emplace_back(part);


        }
    }
}

void addMembranes(ParticleContainer &particles, std::list<FileReader::MembraneData> membranes) {
    for (auto &membrane : membranes) {
        generateMembrane(membrane, particles);
    }
}
