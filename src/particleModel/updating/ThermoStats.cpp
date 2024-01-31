#include "ThermoStats.h"
#include "CellContainerIterators.h"
#include "utils/ArrayUtils.h"
#include <cmath>
#include <spdlog/spdlog.h>




ThermoStats::ThermoStats(CellContainer& container, double delta_t_param,
        std::optional<double> target_temp_param,
        std::optional<double> max_temp_diff_param) 
        : cellContainer(container) , delta_t(delta_t_param), target_temp(target_temp_param) ,
         max_temp_diff(max_temp_diff_param) { }


double ThermoStats::getPotentialEnergy(){
  double sum = 0;
  auto& particle_instances = cellContainer.getInstances();

  for(auto particle1 = particle_instances.begin(); particle1 != particle_instances.end(); particle1++){
        for(auto particle2 = std::next(particle1); particle2 != particle_instances.end(); particle2++){
            double distance = ArrayUtils::L2Norm(particle1->getX() - particle2->getX());
            double mixed_sigma = sigma_mixed[particle1->getType()][particle2->getType()];
            double mixed_epsilon = epsilon_mixed[particle1->getType()][particle2->getType()];
            sum +=  4 * mixed_epsilon * (std::pow(mixed_sigma / distance, 12) - std::pow(mixed_sigma / distance, 6));
        }
    }

  return sum;
}



double ThermoStats::getPressure(){
  double sum = 0;
  auto& instances = cellContainer.getInstances();

  for(auto p1 = instances.begin(); p1 != instances.end(); p1++){
      sum += (p1->getM() * ArrayUtils::scalarProduct(p1->getV(),p1->getV()) +
                        ArrayUtils::scalarProduct(p1->getF(),p1->getX())  );
  }

  return sum;
}

void ThermoStats::initDiffusionCoefficient(){
  particle_positions_previous_iteration = {};
  std::list<Particle>& particles = cellContainer.getInstances();
  for(Particle& particle : particles){
    particle_positions_previous_iteration.emplace_back(&particle, particle.getX());
  }
}

double ThermoStats::getDiffusionCoefficient(){
    double sum = 0;
    size_t amount = 0;

    for(auto iter = begin_CI(); iter != end_CI(); ++iter){
        for(Particle* particle_ptr : *iter){
            auto current_particle_old_position = std::find_if(
                                                              particle_positions_previous_iteration.begin(),
                                                              particle_positions_previous_iteration.end(),
                                                              [particle_ptr](std::pair<Particle*,std::array<double,3>> ptr_and_position)
                                   {return ptr_and_position.first == particle_ptr;});
            if(current_particle_old_position == particle_positions_previous_iteration.end()){
              throw std::runtime_error("position of a particle was not correctly stored for diffusion coefficient");
            }else{
              std::array<double,3> old_position,current_position,domain_size,remove_periodic_boundaries,diff;
              old_position = current_particle_old_position->second;
                current_position = particle_ptr->getX();
              //we have to remove the displacement that happened to
              //the particle's position, because of the periodic boundaries
              domain_size = cellContainer.getDomainBounds();
              std::array<int,3> how_often_periodic_boundaries_crossed = particle_ptr->getBoundariesCrossed();
              //in every direction remove the domain_bounds which were added to the current particle's
              //position to place the particle back into the domain within the simulation
              // e.g. in direction i how_often_periodic_boundaries_crossed[i] = 3
              //then we know that the boundary in direction i was crossed three times in positive
              //direction -> therefore domain_size[i] was subtracted three times from the original position
                remove_periodic_boundaries[0] = -domain_size[0] *
                                static_cast<double>(how_often_periodic_boundaries_crossed[0]);
                remove_periodic_boundaries[1] = -domain_size[1] *
                                static_cast<double>(how_often_periodic_boundaries_crossed[1]);
                remove_periodic_boundaries[2] = -domain_size[2] *
                                static_cast<double>(how_often_periodic_boundaries_crossed[2]);
              //remove_periodic_boundaries now contains exactly what was added/subtracted
              //to the particles position due to periodic boundaries
              //-> subtract remove_periodic_boundaries from the current position
              diff = current_position - remove_periodic_boundaries - old_position;
              sum += ArrayUtils::scalarProduct(diff,diff); //no need to calculate L2norm and then square
              amount++;
              particle_ptr->setBoundariesCrossedZero();
            }
        }
    }
    initDiffusionCoefficient();
    return sum / static_cast<double>(amount);
}


std::vector<double> ThermoStats::getRadialDistributionFunction(double interval_size){
    //simple calculation(pythagoras) of the maximal distance two particles can have in our current system
    double max_distance_of_two_particles = std::sqrt(std::pow(cellContainer.domain_bounds[0], 2) +
                                                        std::pow(cellContainer.domain_bounds[1], 2));
    if(cellContainer.hasThreeDimensions())
        max_distance_of_two_particles = std::sqrt(std::pow(max_distance_of_two_particles, 2) + std::pow(cellContainer.domain_bounds[2], 2));

    std::list<Particle>& particle_instances = cellContainer.particle_instances;
    //first the amount of particle pairs in that interval will be calculated
    //afterwards the rdf statistic, therefore it uses double
    std::vector<size_t> particle_pairs_per_interval(static_cast<size_t>(max_distance_of_two_particles / interval_size + 1), 0.0);
    std::vector<double> rdf_per_interval(static_cast<size_t>(max_distance_of_two_particles / interval_size + 1), 0.0);

    for(auto particle1 = particle_instances.begin(); particle1 != particle_instances.end(); particle1++){
        for(auto particle2 = std::next(particle1); particle2 != particle_instances.end(); particle2++){
            double dist = ArrayUtils::L2Norm(particle1->getX() - particle2->getX());
            size_t interval_index = static_cast<size_t>(dist / interval_size);
            if(interval_index < particle_pairs_per_interval.size())
                //there is particle pair that has a distance that
                //lies within intervals[interval_index]
                particle_pairs_per_interval[interval_index] +=1;
            else
                spdlog::error("rdf Array was to small");
        }
    }



    for(size_t i = 0; i < rdf_per_interval.size();i++ ){
        //calculate rdf statistic for this interval
        rdf_per_interval[i] = static_cast<double>(particle_pairs_per_interval[i]) /
                ( (4.0 * M_PI / 3.0) * ( std::pow((i+1) * interval_size,3) - std::pow(i * interval_size,3)) );
        
    }

    return rdf_per_interval;
}




double ThermoStats::currentTemp(){
  double kinetic_energy = 0;
  size_t amount = 0;

  for(auto iter = begin_CI(); iter != end_CI(); ++iter){
    for(Particle* particle_ptr : *iter){
      const std::array<double,3> &v = particle_ptr->getV();
      double v_squared = v[0] * v[0]  + v[1] * v[1] + v[2] * v[2];
      double m = particle_ptr->getM();
      kinetic_energy += v_squared * m;
      amount++;
    }
  }
  double k_boltzman = 1;
  //calculate the current temperature from the current kinetic energy in the system
  //assuming we only have two kinds of dimensions namely 2 or 3
  double current_temp = kinetic_energy/(static_cast<double>(cellContainer.hasThreeDimensions() ? 3 : 2) * amount * k_boltzman);
  return current_temp;
}

void ThermoStats::applyThermostats(){
  double current_temp = currentTemp();
  double next_temp;
  if(target_temp.has_value())
    next_temp = target_temp.value();  
  else  
   throw std::invalid_argument("applyThermostats was called, altough target temp was not provided ");
  

  //if the temperatur difference would be too big cap it
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
