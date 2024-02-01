#include "CuboidGeneration.h"
#include "utils/MaxwellBoltzmannDistribution.h"
#include "inputHandling/FileReaderProgramArgs.h"

void generateCuboid(FileReader::CuboidData& cuboid, ParticleContainer &particles, size_t dim) {
    uint64_t total_count = cuboid.N1 * cuboid.N2 * cuboid.N3;
    particles.reserve(total_count);
    for(uint64_t z = 0; z < cuboid.N3; z++) {

        for (uint64_t y = 0; y < cuboid.N2; y++) {

            for (uint64_t x = 0; x < cuboid.N1; x++) {
                std::array<double, 3> cords(cuboid.x);
                std::array<double, 3> vel(cuboid.v);
                //add boltzman distributed velocity only if needed
                std::array<double, 3> dist = cuboid.avg_v.has_value()?
                                                maxwellBoltzmannDistributedVelocity(cuboid.avg_v.value(), dim)
                                                : std::array<double, 3>({0,0,0});

                cords[0] += x * cuboid.h;
                cords[1] += y * cuboid.h;
                cords[2] += z * cuboid.h;

                vel[0] += dist[0];
                vel[1] += dist[1];
                vel[2] += dist[2];

                Particle part = Particle{cords, vel, cuboid.m};


                particles.emplace_back(part);

            }
        }
    }
}

void addCuboids(ParticleContainer &particles, std::list<FileReader::CuboidData> cuboids) {
    size_t dim{2};
    double z = cuboids.front().x[2];
    for (auto &cube : cuboids) {
        if(1 < cube.N3 || cube.x[2] != z || cube.v[2] != 0) dim = 3;
    }

    for (auto &cube : cuboids) {
        generateCuboid(cube, particles, dim);
    }
}