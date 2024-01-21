#include "ThermoStats.h"
#include "particleModel/storage/CellContainerIterators.h"
#include "utils/ArrayUtils.h"
#include <cmath>




ThermoStats::ThermoStats(CellContainer& container, double delta_t_param,
        std::optional<double> target_temp_param,
        std::optional<double> max_temp_diff_param) 
        : cellContainer(container) , delta_t(delta_t_param), target_temp(target_temp_param) ,
         max_temp_diff(max_temp_diff_param) { }


double ThermoStats::getPotentialEnergy(){
  double sum = 0;
  auto& instances = cellContainer.getInstances();

  for(auto p1 = instances.begin(); p1 != instances.end(); p1++){
        for(auto p2 = std::next(p1);p2 != instances.end();  p2++){
            double dist = ArrayUtils::L2Norm(p1->getX() - p2->getX());
            sum +=  4 * (std::pow(1/dist,12) - std::pow(1/dist,6));
        }
    }

  return sum;
}

double square(std::array<double,3> x1, std::array<double,3> x2){
  return (x1[0] * x2[0] + x1[1] * x2[1] + x1[2] * x2[2]);
}

double ThermoStats::getPressure(){
  double sum = 0;
  auto& instances = cellContainer.getInstances();

  for(auto p1 = instances.begin(); p1 != instances.end(); p1++){
      sum += (p1->getM() * square(p1->getV(),p1->getV()) +  square(p1->getF(),p1->getX()));
  }

  return sum;
}

void ThermoStats::initDiffCoeff(){
  //std::cout << "Start init diffcoeff\n";
  old_particle_positions = {};
  std::list<Particle>& particles = cellContainer.getInstances();
  for(Particle& particle : particles){
    old_particle_positions.push_back(std::make_pair(&particle,particle.getX()));
  }
  //std::cout << "End init diffcoeff\n";
}

double ThermoStats::diffusionCoeff(){
  //std::cout << "Start diffcoeff\n";
    double sum = 0;
    size_t amt = 0;

    for(auto iter = begin_CI(); iter != end_CI(); ++iter){
        for(Particle* particle_ptr : *iter){
            auto it = std::find_if(old_particle_positions.begin(),old_particle_positions.end(),
                                  [particle_ptr](std::pair<Particle*,std::array<double,3>> ptr_pos){return ptr_pos.first == particle_ptr;});
            if(it == old_particle_positions.end()){
              throw std::runtime_error("position of a particle was not correctly stored for diffusion coefficient");
            }else{
              std::array<double,3> old_pos,current_pos,remove_periodic,diff;
              old_pos = it->second;
              current_pos = particle_ptr->getX();
              remove_periodic = cellContainer.getDomainBounds();
              std::array<int,3> crossed_periodic = particle_ptr->getBoundariesCrossed();
              remove_periodic[0] = remove_periodic[0] * static_cast<double>(crossed_periodic[0]);
              remove_periodic[1] = remove_periodic[1] * static_cast<double>(crossed_periodic[1]);
              remove_periodic[2] = remove_periodic[2] * static_cast<double>(crossed_periodic[2]);
              diff = current_pos - (remove_periodic + old_pos);
              sum += std::pow(ArrayUtils::L2Norm(diff),2);
              amt++;
              particle_ptr->setPeriodicZero();
            }
        }
    } 
    initDiffCoeff();
    //std::cout << "end diffcoeff\n";
    return sum / static_cast<double>(amt);
}


std::vector<double> ThermoStats::radialDistributionFunction(double interval_size){
    double max_dist = std::sqrt(std::pow(cellContainer.domain_bounds[0],2) + std::pow(cellContainer.domain_bounds[1],2));
    if(cellContainer.hasThreeDimensions())
        max_dist = std::sqrt(std::pow(max_dist,2) + std::pow(cellContainer.domain_bounds[2],2));

    std::list<Particle>& instances = cellContainer.particle_instances;
    //first the amount of particle pairs in that interval will be calculated
    //afterwards the rdf statistic, therefore it uses double
    std::vector<size_t> pairs_per_interval(max_dist / interval_size + 1,0.0);
    std::vector<double> rdf_per_interval(max_dist / interval_size + 1,0.0);

    for(auto p1 = instances.begin(); p1 != instances.end(); p1++){
        for(auto p2 = std::next(p1);p2 != instances.end();  p2++){
            double dist = ArrayUtils::L2Norm(p1->getX() - p2->getX());
            size_t interval_index = static_cast<size_t>(dist / interval_size);
            if(interval_index < pairs_per_interval.size())
                //there is particle pair that has a distance that
                //lies within intervals[interval_index]
                pairs_per_interval[interval_index] +=1; 
            else    
                std::cout << "Err\n";
        }
    }



    for(size_t i = 0; i < rdf_per_interval.size();i++ ){
        //calculate rdf statistic for this interval
        rdf_per_interval[i] = static_cast<double>(pairs_per_interval[i]) / 
                    ( (4.0 * M_PI / 3.0) * 
                      ( std::pow((i+1) * interval_size,3) - std::pow(i * interval_size,3)) );
        
    }

    return rdf_per_interval;
}




double ThermoStats::currentTemp(){
  double kinetic_energy = 0;
  size_t amt = 0;

  for(auto iter = begin_CI(); iter != end_CI(); ++iter){
    for(Particle* particle_ptr : *iter){
      const std::array<double,3> &v = particle_ptr->getV();
      double v_squared = v[0] * v[0]  + v[1] * v[1] + v[2] * v[2];
      double m = particle_ptr->getM();
      kinetic_energy += v_squared * m;
      amt++;
    }
  }
  double k_boltzman = 1;
  //calculate the current temperatur from the current kinetic energy in the system
  //assuming we only have two kinds of dimensions namely 2 or 3
  double current_temp = kinetic_energy/((cellContainer.hasThreeDimensions() ? 3 : 2) * amt * k_boltzman);
  return current_temp;
}

void ThermoStats::applyThermostats(){
  double current_temp = currentTemp();
  double next_temp;
  if(target_temp.has_value())
    next_temp = target_temp.value();  
  else  
    throw std::invalid_argument("applyThermostats was called, altough target temp was not provided ");
  

  //if the temperatur diffference would be too big cap it
  double temp_diff = next_temp - current_temp;
  if(max_temp_diff.has_value() && std::abs(temp_diff) > max_temp_diff){
    next_temp = (std::signbit(temp_diff) ? -1 : 1) * max_temp_diff.value() + current_temp;
  }


  // the scaling factor to reach the target temperature
  // assuming this is never negative because we only calculate with kelvin
  double temp_scaling = sqrt(next_temp/current_temp);

  //apply scaling
  for(auto iter = begin_CI(); iter != end_CI(); ++iter){
      for(Particle* particle_ptr : *iter){
        std::array<double,3> v = particle_ptr->getV();
        particle_ptr->setV(temp_scaling * v);
      }
  }
}
