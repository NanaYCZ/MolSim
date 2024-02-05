#include <iostream>
#include "CellContainer.h"
#include "particleModel/updating/CellContainerIterators.h"
#include "particleModel/updating/CellCalculator.h"
#include <cmath>
#include <spdlog/spdlog.h>

CellContainer::CellContainer(double d_width, double d_height, double d_depth, double r_cutoff, double cell_size)
            : cell_size(cell_size), domain_bounds({d_width, d_height, d_depth})
              {

    domain_max_dim = {static_cast<dim_t>(d_width / cell_size + 1),
                      static_cast<dim_t>(d_height / cell_size + 1),
                      static_cast<dim_t>(d_depth / cell_size + 1)};
    //check if modulo would be 0
    if(isApproximatelyEqual(std::fmod(d_width, cell_size), 0.0)) {
        --domain_max_dim[0];
    }
    if(isApproximatelyEqual(std::fmod(d_height, cell_size), 0.0)) {
        --domain_max_dim[1];
    }
    if(isApproximatelyEqual(std::fmod(d_depth, cell_size), 0.0)) {
        --domain_max_dim[2];
    }
    particles = {};
    particles.resize(static_cast<dim_t>(domain_max_dim[0] + 2),
                     std::vector<std::vector<std::vector<Particle*>>>(
                             static_cast<dim_t>(domain_max_dim[1] + 2),
                             std::vector<std::vector<Particle*>>(
                                     static_cast<dim_t>(domain_max_dim[2] + 2)
                             )
                     )
    );

    for (auto cells = begin_CellIterator(); cells != end_CellIterator(); ++cells) {
        auto &current_cell = *cells;
        current_cell = {};
    }
    patterns.clear();

    if (cell_size < r_cutoff) {
        comparing_depth = std::ceil(r_cutoff / cell_size);
    }

    if(domain_max_dim[2] <= 1) {
        domain_max_dim[2] = 1;

        if(domain_max_dim[0] <= comparing_depth) {
            throw std::invalid_argument("Domain width is too small for the r_cutoff to cell size ratio."
                                        "Consider increasing the cell size or using ParticleContainer instead.");

        } else if(domain_max_dim[1] < comparing_depth) {
            throw std::invalid_argument("Domain height is too small for the r_cutoff to cell size ratio."
                                        "Consider increasing the cell size.");
        }

        three_dimensions = false;

    } else {
        if(domain_max_dim[0] <= comparing_depth) {
            throw std::invalid_argument("Domain width is too small for the r_cutoff to cell size ratio."
                                        "Consider increasing the cell size.");

        } else if(domain_max_dim[1] <= comparing_depth) {
            throw std::invalid_argument("Domain height is too small for the r_cutoff to cell size ratio."
                                        "Consider increasing the cell size.");

        }
        else if(domain_max_dim[2] < comparing_depth) {
            throw std::invalid_argument("Domain depth is too small for the r_cutoff to cell size ratio."
                                        "Consider increasing the cell size.");
        }

        three_dimensions = true;
    }

    //precalculate all patterns
    std::array<dim_t, 3> pattern{};

    if(three_dimensions) {
        while(setNext3dPattern(pattern)) patterns.emplace_back(pattern);

    } else {
        while(setNext2dPattern(pattern)) patterns.emplace_back(pattern);
    }

    patterns.shrink_to_fit();

    std::ostringstream out_str;
    out_str << "Domain cells x: 1 - " << (domain_max_dim[0]) <<  " y: 1 - " << (domain_max_dim[1]) << " z: 1 - "  << (domain_max_dim[2]) << std::endl;
    SPDLOG_INFO(out_str.str());
}

CellContainer::~CellContainer() {}

enum direction_status {
    first_subset, second_subset, third_subset
};

bool CellContainer::setNext3dPattern(std::array<dim_t, 3> &pattern) {
    static direction_status status = first_subset;

    switch (status) {
        case (first_subset):
            // (0, 0, 1 to depth)
            if (pattern[2] < comparing_depth) {
                ++pattern[2];
            } else {
                pattern[0] = 1;
                pattern[2] = -comparing_depth;
                status = second_subset;
            }
            break;

        case (second_subset):
            // (1 to depth, 0, -depth to depth)
            if (pattern[0] < comparing_depth) {
                ++pattern[0];

            } else if (pattern[2] < comparing_depth) {
                pattern[0] = 1;
                ++pattern[2];

            } else {
                pattern[0] = -comparing_depth;
                pattern[1] = 1;
                pattern[2] = -comparing_depth;
                status = third_subset;
            }
            break;

        case (third_subset):
            // (-depth to depth, 1 to depth, -depth to depth)
            if (pattern[0] < comparing_depth) {
                ++pattern[0];

            } else if (pattern[1] < comparing_depth) {
                pattern[0] = -comparing_depth;
                ++pattern[1];

            } else if (pattern[2] < comparing_depth) {
                pattern[0] = -comparing_depth;
                pattern[1] = 1;
                ++pattern[2];
            } else {
                //finished
                status = first_subset;
                return false;
            }
            break;
    }
    return true;
}

bool CellContainer::setNext2dPattern(std::array<dim_t, 3> &pattern) {
    static direction_status status = first_subset;

    switch (status) {
        case (first_subset):
            // (1 to depth, 0, 0)
            if (pattern[0] < comparing_depth) {
                ++pattern[0];

            } else {
                pattern[0] = -comparing_depth;
                pattern[1] = 1;
                status = second_subset;
            }
            break;

        case (second_subset):
            // (-depth to depth, 1 to depth, 0)
            if (pattern[0] < comparing_depth) {
                ++pattern[0];

            } else if (pattern[1] < comparing_depth) {
                pattern[0] = -comparing_depth;
                ++pattern[1];

            } else {
                //finished
                status = first_subset;
                return false;
            }
            break;
    }
    return true;
}

void CellContainer::addParticle(std::array<double, 3> x_arg, std::array<double, 3> v_arg, double m_arg) {
    addParticle(x_arg, v_arg, m_arg, default_grid_index, 1,1,sigma_mixed[0][0], epsilon_mixed[0][0]);
}

void CellContainer::addParticle(const Particle& particle,double sigma, double epsilon){
    addParticle(particle.getX(), particle.getV(), particle.getM(), default_grid_index, 1,1,sigma, epsilon);
}

void CellContainer::addParticle(std::array<double, 3> x_arg, std::array<double, 3> v_arg, std::array<int,3> grid, double rz, double fp, double m_arg){
    addParticle(x_arg, v_arg, m_arg, grid,rz,fp,sigma_mixed[0][0], epsilon_mixed[0][0]);
}

void CellContainer::addParticle(std::array<double, 3> x_arg, std::array<double, 3> v_arg, double m_arg, double sigma, double epsilon)
{
    addParticle(x_arg, v_arg, m_arg, default_grid_index,1,1,sigma, epsilon);
}

void CellContainer::addParticle(std::array<double, 3> x_arg, std::array<double, 3> v_arg,
                                double m_arg, std::array<int,3> grid, double rz, double fp,double sigma, double epsilon) {
    if(domain_bounds[0] <= x_arg[0] || domain_bounds[1] <= x_arg[1] || (domain_bounds[2] <= x_arg[2] && three_dimensions) ||
      x_arg[0] < 0 || x_arg[1] < 0 || x_arg[2] < 0)
    {
        std::cerr << "Wanted to add at: " << x_arg[0] << " , " << x_arg[1] << " , " << x_arg[2] << "\n";
        throw std::invalid_argument("The provided coordinates are outside the domain borders.");
    }

    int type = 0;
    static std::vector<std::array<double,2>> types;

    //find the existing type/sigma-epsilon-pair
    auto it = std::find_if(types.begin(), types.end(), [sigma, epsilon](std::array<double,2> pair) {
       return pair[0] == sigma && pair[1] == epsilon;
    });

    if(it != types.end()) {
        type = std::distance(types.begin(), it);
    } else {
        //add new type and update Lorentz-Berthelot mixing matrices
        type = types.size();
        types.emplace_back(std::array<double,2>{sigma, epsilon});

        sigma_mixed.clear();
        sigma_mixed.resize(types.size(), std::vector<double>(types.size()));
        epsilon_mixed.clear();
        epsilon_mixed.resize(types.size(), std::vector<double>(types.size()));

        for (int i = 0; i < types.size(); ++i) {
            for (int j = 0; j < types.size(); ++j) {

                sigma_mixed[i][j] = (types[i][0] + types[j][0]) / 2;
                epsilon_mixed[i][j] = std::sqrt(types[i][1] * types[j][1]);
            }
        }
    }
    Particle newp = Particle(x_arg, v_arg, m_arg, type, grid, rz, fp);
//    std::cout<<"grid:"<<newp.getGrid()<<std::endl;
//    newp.setGrid(grid);
//    std::cout<<"grid after:"<<newp.getGrid()<<std::endl;
    particle_instances.push_back(newp);
//    std::cout << "SPECIAL:" << newp.getSpecial() << std::endl;
//    newp.setSpecial(specialF);
//    std::cout << "newSPECIAL:" << newp.getSpecial() << std::endl;
//    particle_instances.emplace_back(x_arg, v_arg, m_arg, type, grid, rz, fp, specialF);
//    std::cout<<particle_instances<<std::endl;
}


bool CellContainer::isApproximatelyEqual(double a, double b, double epsilon) {
    return std::abs(a - b) < epsilon;
}


void CellContainer::createPointers(){
    particle_instances.shrink_to_fit();

    for(Particle& particle : particle_instances){
        std::array<dim_t , 3> pos;
        std::array<double,3> x_arg = particle.getX();
        allocateCellFromPosition(x_arg, pos);
        particles.at(pos[0]).at(pos[1]).at(pos[2]).push_back(&particle);
        particle_amount++;
    }
}


void CellContainer::plotParticles(outputWriter::VTKWriter &writer) {
    for(Particle& particle : particle_instances){
        writer.plotParticle(particle);
    }
}


std::string CellContainer::to_string() {
    std::ostringstream out_str;

    size_t amt = 0;

    out_str << "The actual domain has  \n" << particles.size() << " cells in x dir. \n" << particles[0].size()
            << " cells in y dir. \n" << particles[0][0].size() << " cells in z dir." << std::endl;
    out_str << "The actual domain is from  \nx: 1 - " << (domain_max_dim[0]) << "(domain_max_dim[0])\ny: 1 - "
            << (domain_max_dim[1]) << "(domain_max_dim[1])\nz: 1 - " << (domain_max_dim[2]) << "(domain_max_dim[2])"
            << std::endl;
    out_str << "Are we in 3d?: " << (three_dimensions ? "Yes" : "No") << std::endl;
    out_str << "cell_size: " << cell_size << std::endl;
    out_str << "comparing_depth: " << comparing_depth << std::endl;
    out_str << "domain_bounds [0]:" << domain_bounds[0] << " [1]:" << domain_bounds[1] << " [2]:" << domain_bounds[2]
            << std::endl;

    for (CellIterator it = begin_CellIterator(); it != end_CellIterator(); ++it) {
        out_str << "The cell with index x=" << it.position()[0] << " y=" << it.position()[1] << " z="
                << it.position()[2] << std::endl;
        out_str << "Has the following Particles: " << std::endl;

        for (auto *particle: *it) {
            out_str << (*particle).toString() << std::endl;
            amt++;
        }
        out_str << "\n\n";
    }

    out_str << "in total amt: " << size() << std::endl;

    return out_str.str();
}


std::list<Particle> CellContainer::to_list(){
    return {particle_instances.begin(),particle_instances.end()};
}


size_t CellContainer::size() {
    return particle_instances.size();
}

std::array<dim_t, 3> CellContainer::domain_max_dim{};
std::vector<std::array<dim_t,3>> CellContainer::patterns{};
std::vector<std::vector<std::vector<std::vector<Particle*>>>> CellContainer::particles{};

std::vector<Particle> CellContainer::getParticleInstances(){
    return particle_instances;
}

std::vector<std::array<double,3>> CellContainer::getSpecialPositions(){
    return specialPositions;
}
int CellContainer::getSpecialTime(){
    return spetialTime;
}
std::vector<std::array<double,3>> CellContainer::getSpecialStrength() {
    return specialStrength;
}
void CellContainer::pushbackSpecialPosition(std::array<double,3> sp){
    specialPositions.push_back(sp);
}
void CellContainer::setSpecialTime(int st){
    spetialTime=st;
}
void CellContainer::pushbackSpecialForce(std::array<double,3> sf){
    specialStrength.push_back(sf);
}










