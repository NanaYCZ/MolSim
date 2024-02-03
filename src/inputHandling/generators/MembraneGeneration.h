#pragma once

#include "inputHandling/FileReader.h"

void generateMembrane(FileReader::MembraneData& membrane, CellContainer& container);
void addMembranes(CellContainer &container, std::list<FileReader::MembraneData> membranes);