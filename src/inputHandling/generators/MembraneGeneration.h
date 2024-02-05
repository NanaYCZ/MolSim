#pragma once

#include "inputHandling/FileReader.h"

/**
 * @brief creates all particles of the membrane and adds them into the ParticleContainer
 *
 * iterates over each dimension N of the membrane and creates particles with given
 * dimension and the cuboids initial velocity, mass and coordinates. The coordinates
 * offset is based on the particles location within the cuboid and a distance h
 * to the particles around it. Grid of the particles is based on their initial index.
 *
 * @param membrane data for the generation
 * @param container reference to add the cuboids particles to
 */
void generateMembrane(FileReader::MembraneData& membrane, CellContainer& container);

/**
 * @brief determines global dimension of the membranes particles and iteratively generates membranes
 *
 * @param container reference to add particles to
 * @param membranes all the membranes data to iterate over
 */
void addMembranes(CellContainer &container, std::list<FileReader::MembraneData> membranes);


/**
 * @brief receives a list of special forces and iteratively generates special forces
 * including their position, strength adn last time
 *
 * @param container reference to add special forces to
 * @param specialForces all the special Forces data to iterate over
 */
void addSpecialForces(CellContainer &container, std::list<FileReader::SpecialForcesData> specialForces);