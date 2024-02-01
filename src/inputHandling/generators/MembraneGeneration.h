#pragma once

#include "inputHandling/FileReader.h"

void generateMembrane(FileReader::MembraneData& membrane, ParticleContainer &particles);

void addMembranes(ParticleContainer &particles, std::list<FileReader::MembraneData> membranes);