#include "SphereGeneration.h"
#include "utils/ArrayUtils.h"
#include "utils/MaxwellBoltzmannDistribution.h"
#include "inputHandling/FileReaderProgramArgs.h"
#include "cmath"

void addSpheres(ParticleContainer& particles, std::list<FileReader::SphereData> spheres, size_t dim) {
    if (dim == 2 ) {
        for (auto &sphere : spheres) {
            generateSpheresMethod2D(sphere, particles);
        }
    }
    if (dim == 3 ) {
        for (auto &sphere : spheres) {
            generateSpheresMethod3D(sphere, particles);
        }
    }

}

void generateSpheresMethod2D(FileReader::SphereData& sphere, ParticleContainer& particles){
    
    std::array<double, 3> cords{0,0,0};

    double radiuslength = sphere.radius * sphere.meshWidth;
    std::array<double,2> BottomLeftCorner {sphere.CenterPosition[0]-radiuslength, sphere.CenterPosition[1]-radiuslength};
    for (auto x=0; x<=2*sphere.radius;x++){
        for (auto y=0; y<=2*sphere.radius;y++){
            std::array<double,2> position {x*sphere.meshWidth+BottomLeftCorner[0],y*sphere.meshWidth+BottomLeftCorner[1]};
            if (pow(position[0]-sphere.CenterPosition[0],2)+ pow(position[1]-sphere.CenterPosition[1],2)<pow(radiuslength,2)+0.1){
                std::array<double, 3> dist = sphere.avg_v.has_value()? 
                                                maxwellBoltzmannDistributedVelocity(sphere.avg_v.value(), 2)
                                                : std::array<double, 3>({0,0,0});
                std::array<double, 3> vel(sphere.Velocity + dist);
                cords[0]=position[0];
                cords[1]=position[1];
                Particle part = Particle{cords, vel, sphere.mass};


                particles.emplace_back(part);
            }
        }
    }
}

void generateSpheresMethod3D(FileReader::SphereData& sphere, ParticleContainer& particles){
    std::array<double, 3> cords{0,0,0};

    double radiuslength = sphere.radius * sphere.meshWidth;
    std::array<double,3> BottomLeftBehindCorner {sphere.CenterPosition[0]-radiuslength, sphere.CenterPosition[1]-radiuslength, sphere.CenterPosition[2]-radiuslength};
    for (auto x=0; x<=2*sphere.radius;x++){
        for (auto y=0; y<=2*sphere.radius;y++){
            for (auto z=0; z<=2*sphere.radius;z++) {
                std::array<double, 3> position{x * sphere.meshWidth + BottomLeftBehindCorner[0],
                                               y * sphere.meshWidth + BottomLeftBehindCorner[1],
                                               z * sphere.meshWidth + BottomLeftBehindCorner[2]};
                if (pow(position[0] - sphere.CenterPosition[0], 2) + pow(position[1] - sphere.CenterPosition[1], 2)
                    + pow(position[2] - sphere.CenterPosition[2], 2) <pow(radiuslength, 2) + 0.1) {
                    std::array<double, 3> dist = sphere.avg_v.has_value()? 
                                                maxwellBoltzmannDistributedVelocity(sphere.avg_v.value(), 2)
                                                : std::array<double, 3>({0,0,0});
                    std::array<double, 3> vel(sphere.Velocity + dist);
                    cords[0] = position[0];
                    cords[1] = position[1];
                    cords[2] = position[2];

                    Particle part = Particle{cords, vel, sphere.mass};


                    particles.emplace_back(part);
                }
            }
        }
    }
}