#pragma once

#include "inputHandling/FileReader.h"

void generateMembrane(FileReader::MembraneData& membraneData, CellContainer& container, size_t dim);

void addSpecialForcesToTheMembrane(CellContainer& container, size_t dim);