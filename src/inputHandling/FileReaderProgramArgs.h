 #include "FileReader.h"
 

/**
 * @struct SphereData
 * 
 * contains relevant information for one sphere
*/
 struct FileReader::SphereData {
    std::array<double, 3> CenterPosition;
    std::array<double, 3> Velocity;
    double mass;
    double radius;
    double meshWidth;
    double sigma;
    double epsilon;

    //by default no maxwell-boltzmann is applied
    std::optional<double> avg_v;

    std::string to_string() const{
     auto sphereData = (*this);
    std::ostringstream oss;

    oss << "  SphereData:" << std::endl;
    oss << "  center: (" << sphereData.CenterPosition[0] << ", "
        << sphereData.CenterPosition[1] << ", " << sphereData.CenterPosition[2] << ")" << std::endl;
    oss << "  velocity: (" << sphereData.Velocity[0] << ", "
        << sphereData.Velocity[1] << ", " << sphereData.Velocity[2] << ")" << std::endl;
    if(avg_v.has_value())
      oss << "  v_avg: " << avg_v.value() << std::endl;
    oss << "  mass: " << sphereData.mass << std::endl;
    oss << "  radius: " << sphereData.radius << std::endl;
    oss << "  mesh width: " << sphereData.meshWidth << std::endl;
    oss << "  sigma: " << sphereData.sigma << std::endl;
    oss << "  epsilon: " << sphereData.epsilon << std::endl;

    return oss.str();
    }

    bool operator==(const SphereData& other) const {
        return (CenterPosition == other.CenterPosition &&
                Velocity == other.Velocity &&
                mass == other.mass &&
                radius == other.radius &&
                meshWidth == other.meshWidth &&
                sigma == other.sigma &&
                epsilon == other.epsilon);
    }


  };


/**
 * @struct CuboidData
 * 
 * contains relevant information for one cuboid
 * 
*/
struct FileReader::CuboidData {
    /// initial velocity and position vectors
    std::array<double, 3> x, v;


    /// N1: amount of particles along dimension 1
    /// N2: amount of particles along dimension 2
    /// N3: amount of particles along dimension 3
    uint64_t N1, N2, N3;

    /// Mass m of the particles in the cuboid
    /// Mesh width h
    double m, h;

    /// sigma and epsilon parameters for the force calculation
    /// between particles of this cuboid
    double sigma, epsilon;

    /// Average velocity default 0 means by default no Maxwell-boltzmann is applied
    std::optional<double> avg_v;

    /**
     * @brief Convert CuboidData to a string
     *
     * @return String representation of CuboidData
     */
    std::string to_string() const {
        std::stringstream ss;

        ss << "CuboidData:" << std::endl;
        ss << "  x: (" << x[0] << ", " << x[1] << ", " << x[2] << ")"
           << std::endl;
        ss << "  v: (" << v[0] << ", " << v[1] << ", " << v[2] << ")"
           << std::endl;
        if (avg_v.has_value())
            ss << "  v_avg: " << avg_v.value() << std::endl;
        ss << "  N1: " << N1 << std::endl;
        ss << "  N2: " << N2 << std::endl;
        ss << "  N3: " << N3 << std::endl;
        ss << "  m: " << m << std::endl;
        ss << "  h: " << h << std::endl;
        ss << "  sigma: " << sigma << std::endl;
        ss << "  epsilon: " << epsilon << std::endl;

        return ss.str();
    }

    /**
    * @brief Compare two CuboidData structs. Used for Testing
    *
    * @return True if this CuboidData struct has the same content as other
    */
    bool operator==(const CuboidData &other) const {
        return (x == other.x && v == other.v && N1 == other.N1 &&
                N2 == other.N2 && N3 == other.N3 && m == other.m &&
                h == other.h && sigma == other.sigma &&
                epsilon == other.epsilon);
    }

};

    /**
 * @struct MembraneData
 *
 * contains relevant information for one Membrane
 *
*/
struct FileReader::MembraneData {
        /// initial velocity and position vectors
        std::array<double, 3> x, v;

        double a, f;

        /// N1: amount of particles along dimension 1
        /// N2: amount of particles along dimension 2
        /// N3: amount of particles along dimension 3
        int N1, N2, N3;

        /// Mass m of the particles in the Membrane
        /// Mesh width h
        double m, h;

        /// sigma and epsilon parameters for the force calculation
        /// between particles of this Membrane
        double sigma, epsilon;

        /**
         * @brief Convert MembraneData to a string
         *
         * @return String representation of MembraneData
         */
        std::string to_string() const {
            std::stringstream ss;

            ss << "MembraneData:" << std::endl;
            ss << "  x: (" << x[0] << ", " << x[1] << ", " << x[2] << ")"
               << std::endl;
            ss << "  v: (" << v[0] << ", " << v[1] << ", " << v[2] << ")"
               << std::endl;
            ss << "  N1: " << N1 << std::endl;
            ss << "  N2: " << N2 << std::endl;
            ss << "  N3: " << N3 << std::endl;
            ss << "  m: " << m << std::endl;
            ss << "  h: " << h << std::endl;
            ss << "  sigma: " << sigma << std::endl;
            ss << "  epsilon: " << epsilon << std::endl;

            return ss.str();
        }

    /**
     * @brief Compare two MembraneData structs. Used for Testing
     *
     * @return True if this MembraneData struct has the same content as other
     */
    bool operator==(const MembraneData &other) const {
      return (a == other.a && f == other.f && x == other.x && v == other.v && N1 == other.N1 &&
              N2 == other.N2 && N3 == other.N3 && m == other.m &&
              h == other.h && sigma == other.sigma &&
              epsilon == other.epsilon);
    }
};


/**
 * @struct ProgramArgs 
 * @brief contains all arguments relevant for starting the program
 * 
*/
 struct FileReader::ProgramArgs {
    // arguments of the simulation

    bool calculate_thermostats = false;

    double delta_t;
    double t_end;
    double cut_off_radius;
    double cell_size;
    double gravity_factor;
    double init_temp = 0;
    std::string force_type = "LJ";
    CellContainer::concurrency_strategy parallelization_version = CellContainer::serial;
    std::optional<int> choose_amount_threads = std::nullopt;
    std::optional<int> diff_frequency = std::nullopt;
    std::optional<std::pair<double,int>> rdf_interval_and_frequency = std::nullopt;
    std::optional<double> max_temp_diff = std::nullopt;
    std::optional<double> target_temp = std::nullopt;
    int thermo_stat_frequency = 0;
    std::array<CellContainer::boundary_conditions,6> boundaries;
    std::array<double,3> domain_dimensions;

    std::optional<std::string> checkpoint_input_file;
    std::optional<std::string> checkpoint_output_file;



    std::string file_basename = "out";
    size_t write_frequency = 10;

    // spheres and cuboids to simulate
    std::list<CuboidData> cuboids;
    std::list<SphereData> spheres;
    std::list<MembraneData> membranes;
  

    std::string to_string() const{
    
    std::ostringstream oss;

    oss << "Delta_t: " << delta_t << std::endl;
    oss << "T_end: " << t_end << std::endl;
    oss << "File basename: " << file_basename << std::endl;
    oss << "Write frequency: " << write_frequency << std::endl;
    oss << "Domain bounds: [" << std::endl;
    oss << "width: "<< domain_dimensions[0] << std::endl;
    oss << "Height: "<< domain_dimensions[1] << std::endl;
    oss << "Depth: "<< domain_dimensions[2] << std::endl;
    oss << "]" << std::endl;
    oss << "cut_of_radius: " << cut_off_radius << std::endl;
    oss << "cell_size: " << cell_size << std::endl;
    oss << "gravity_factor: " << gravity_factor << std::endl;
    oss << "init_temp: " << init_temp << std::endl;
    oss << "max_temp_diff: " << (max_temp_diff.has_value() ? std::to_string(*max_temp_diff) : "nullopt") << std::endl;
    oss << "target_temp: " << (target_temp.has_value() ? std::to_string(*target_temp) : "nullopt") << std::endl;
    oss << "thermo_stat_frequency: " << thermo_stat_frequency << std::endl;

    oss << "Boundary conditions: [" << std::endl;
    int side = 0;
    for (const auto condition : boundaries) {
        std::string condition_name;
        std::string side_name;
        switch(condition) {
            case CellContainer::boundary_conditions::reflective:
                condition_name = "reflective";
                break;
            case CellContainer::boundary_conditions::outflow:
                condition_name = "outflow";
                break;
            case CellContainer::boundary_conditions::periodic:
                condition_name = "periodic";
                break;
            default:
                condition_name = "undefined";

        }
        switch(side) {
            case(0):
                side_name = "Positive Z: ";
                break;
            case(1):
                side_name = "Negative Z: ";
                break;
            case(2):
                side_name = "Positive X: ";
                break;
            case(3):
                side_name = "Negative X: ";
                break;
            case(4):
                side_name = "Positive Y: ";
                break;
            case(5):
                side_name = "Negative Y: ";
                break;
        }

        oss << side_name << condition_name << std::endl;
        side++;
    }
    oss << "]" << std::endl;

    oss << "Spheres: [" << std::endl;
    for (const auto& sphere : spheres) {
        oss <<  sphere.to_string() << std::endl;
    }
    oss << "]" << std::endl;

    oss << "Cuboids: [" << std::endl;
    for (const auto& cuboid : cuboids) {
        oss << "  " << cuboid.to_string() << std::endl;
    }
    oss << "]" << std::endl;

    return oss.str();
    }

    bool operator==(const ProgramArgs& other) const {
        if(cuboids.size() == other.cuboids.size() && spheres.size() == other.spheres.size()){
          auto cuboids_it = cuboids.begin();
          auto cuboids_it2 = other.cuboids.begin();
          

          for(; cuboids_it != cuboids.end() && cuboids_it2 != other.cuboids.end() ; ++cuboids_it , ++cuboids_it2){
            if( !(*cuboids_it ==  * cuboids_it2)){
                std::cout << "Comp failed bc these two were not equal:" << (*cuboids_it).to_string() << (*cuboids_it2).to_string() << std::endl;
                return false;
            }
          }

          auto spheres_it = spheres.begin();
          auto spheres_it2 = other.spheres.begin();

          for(;spheres_it != spheres.end() && spheres_it2 != other.spheres.end(); ++spheres_it , ++spheres_it2){
            if(!(*spheres_it == *spheres_it2)){
              std::cout << "Comp failed bc these two were not equal:" << (*spheres_it).to_string() << (*spheres_it2).to_string() << std::endl;
                return false;
            }
          }
        }else{
          std::cout << "Comp failed bc the structs had a different amount of cuboids or spheres" << std::endl;
          return false;
        }

    if (!(delta_t == other.delta_t)) {
        std::cout << "Comparison failed: delta_t" << std::endl;
        return false;
    }
    if (!(t_end == other.t_end)) {
        std::cout << "Comparison failed: t_end" << std::endl;
        return false;
    }
    if (!(cut_off_radius == other.cut_off_radius)) {
        std::cout << "Comparison failed: cut_off_radius" << std::endl;
        return false;
    }
    if (!(cell_size == other.cell_size)) {
        std::cout << "Comparison failed: cell_size" << std::endl;
        return false;
    }
    if (!(gravity_factor == other.gravity_factor)) {
        std::cout << "Comparison failed: gravity_factor" << std::endl;
        return false;
    }
    if (!(init_temp == other.init_temp)) {
        std::cout << "Comparison failed: init_temp" << std::endl;
        return false;
    }

    // Seperate if clauses for all the other comparisons using ==
    if (!(target_temp == other.target_temp)) {
        std::cout << "Comparison failed: target_temp" << std::endl;
        return false;
    }
    if (!(max_temp_diff == other.max_temp_diff)) {
        std::cout << "Comparison failed: max_temp_diff" << std::endl;
        return false;
    }
    if (!(thermo_stat_frequency == other.thermo_stat_frequency)) {
        std::cout << "Comparison failed: thermo_stat_frequency" << std::endl;
        return false;
    }
    if (!(boundaries == other.boundaries)) {
        std::cout << "Comparison failed: boundaries" << std::endl;
        return false;
    }
    if (!(domain_dimensions == other.domain_dimensions)) {
        std::cout << "Comparison failed: domain_dimensions" << std::endl;
        return false;
    }
    if (!(checkpoint_input_file == other.checkpoint_input_file)) {
        std::cout << "Comparison failed: checkpoint_input_file" << std::endl;
        return false;
    }
    if (!(checkpoint_output_file == other.checkpoint_output_file)) {
        std::cout << "Comparison failed: checkpoint_output_file" << std::endl;
        return false;
    }
    if (!(file_basename == other.file_basename)) {
        std::cout << "Comparison failed: file_basename" << std::endl;
        return false;
    }
    if (!(write_frequency == other.write_frequency)) {
        std::cout << "Comparison failed: write_frequency" << std::endl;
        return false;
    }
    if (!(cuboids == other.cuboids)) {
        std::cout << "Comparison failed: cuboids" << std::endl;
        return false;
    }
    if (!(spheres == other.spheres)) {
        std::cout << "Comparison failed: spheres" << std::endl;
        return false;
    }

    return true;        
    }



 };