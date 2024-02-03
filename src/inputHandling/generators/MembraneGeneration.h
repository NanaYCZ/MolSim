#pragma once

#include "inputHandling/FileReader.h"

void generateMembrane(FileReader::MembraneData& membrane, FileReader::SpecialForcesData& specialForces,CellContainer& container);
void addMembranes(CellContainer &container, std::list<FileReader::SpecialForcesData>  specialForces,std::list<FileReader::MembraneData> membranes);