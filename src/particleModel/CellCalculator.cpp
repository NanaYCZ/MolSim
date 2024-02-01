#include <iostream>
#include "CellCalculator.h"

namespace calculator {

    void CellCalculator::calcFWithinCell(CellContainerIterators &cell) {
        for (auto it = cell.begin(); it != cell.end(); ++it) {
            for (auto it2 = it; it2 != cell.end(); ++it2) {
                if (*it != *it2) {
                    calcForce(*it, *it2, true);
                }
            }
        }
    }

    void CellCalculator::calcX(ParticleContainer &container) const {
        auto &gridLC = static_cast<CellContainer &>(container);

        if (gridLC.getStrategy() == CellContainer::serial) {
            for (auto &curCell: gridLC.getGrid()) { // Loop through every cell
                // If a cell is empty, skip it
                if (curCell.getParticles().empty()) continue;
                calcNewX(gridLC, curCell);
            }
        } else {
#ifdef _OPENMP
#pragma omp parallel for default(none) schedule(guided) shared(gridLC)
#endif //_OPENMP
            for (auto &curCell: gridLC.getGrid()) { // Loop through every cell
                // If a cell is empty, skip it
                if (curCell.getParticles().empty()) continue;
                calcNewX(gridLC, curCell);
            }
        }
    }

    void CellCalculator::moveParticles(CellContainer &grid) {
        std::array<int, 3> currentIndexes{};
        for (auto &curCell: grid.grid) {
            currentIndexes = curCell.getIndex();

            // checks if it is a border cell
            if (currentIndexes[0] == 0 || currentIndexes[0] == grid.getDim()[0] - 1 ||
                currentIndexes[1] == 0 || currentIndexes[1] == grid.getDim()[1] - 1 ||
                currentIndexes[2] == 0 || currentIndexes[2] == grid.getDim()[2] - 1) {

                for (auto &it: curCell) {

                    // Checks whether any particle has crossed the boundaries
                    for (int d = 0; d < 3; ++d) {
                        if (it->getX()[d] < 0) {
                            // outflow, removing the particle
                            if (std::get<0>(grid.getBorders(currentIndexes, d)) ==
                                CellContainer::outflow) {
                                spdlog::info("Removing Particle");
                                it->valid = false;
                                break;
                            }
                                // periodic
                            else if (std::get<0>(grid.getBorders(currentIndexes, d)) ==
                                     CellContainer::periodic) {
                                // set X to the opposite site
                                spdlog::info("Particle was at d: {} and position {} {} {} now at {}", d,
                                             it->getX()[0], it->getX()[1], it->getX()[2],
                                             grid.getLenDim()[d] + it->getX()[d]);
                                it->setX(d, grid.getLenDim()[d] + it->getX()[d]);
                            }
                        } else if (it->getX()[d] >= grid.getLenDim()[d]) {
                            // outflow, removing the particle
                            if (std::get<0>(grid.getBorders(currentIndexes, d)) ==
                                CellContainer::outflow) {
                                spdlog::info("Removing Particle");
                                it->valid = false;
                                break;
                            }
                                // periodic
                            else if (std::get<0>(grid.getBorders(currentIndexes, d)) ==
                                     CellContainer::periodic) {
                                // set X to the opposite site
                                spdlog::info("Particle was at d: {} and position {} {} {} now at {}", d,
                                             it->getX()[0], it->getX()[1], it->getX()[2],
                                             it->getX()[d] - grid.getLenDim()[d]);
                                it->setX(d, it->getX()[d] - grid.getLenDim()[d]);
                            }
                        }
                    }
                }
            }
        }
    }

    void CellCalculator::calcFSingleNeighbor(CellContainer &grid, const std::array<int, 3> &neighbor,
                                             Particle *p, bool newton) {
        // Loops through every particle of the neighbor
        for (auto &p_other: grid.grid[grid.index(neighbor)]) {
            calcForce(p, p_other, newton);
        }
    }


    void CellCalculator::calcPerNeighbors(CellContainer &grid, const std::array<int, 3> &neighbors, Particle *p,
                                          const std::array<double, 3> &mirror, bool newton) const {
        // Loop through the neighbors
        for (auto &p_other: grid.grid[grid.index(neighbors)]) {
            // This if is important if the domain only has one cell
            if (p != p_other) {
                auto mirroredX = p_other->getX() + mirror;
                double sqrd_dist = 0;
                for (int i = 0; i < DIM; i++) {
                    sqrd_dist += CellCalculator::sqr(mirroredX[i] - p->getX()[i]);
                }
                if (sqrd_dist <= CellCalculator::sqr(rCut)) {
                    std::array<double, 3> force{};

                    double s = sqr(sigma) / sqrd_dist;
                    s = s * s * s; // s = sigma⁶/dist⁶
                    double f = 24 * epsilon / sqrd_dist * (s - 2 * s * s);
                    force = f * (mirroredX - p->getX());


                    p->setF(p->getF() + force);
                    if (newton) {
                        p_other->setF(p_other->getF() - force);
                    }
                }
            }
        }
    }

    void CellCalculator::calcFWithNeighbors(CellContainer &grid,
                                            Particle *p,
                                            const std::array<int, 3> &neighbor,
                                            bool newton) {

        if (neighbor[0] < grid.getDim()[0] && neighbor[1] < grid.getDim()[1] &&
            neighbor[2] < grid.getDim()[2] && neighbor[0] >= 0 && neighbor[1] >= 0 &&
            neighbor[2] >= 0) {

            calcFSingleNeighbor(grid, neighbor, p, newton);
            //}
        } else if (grid.isPeriodic(neighbor)) {
            std::array<int, 3> neigh;
            for (int i = 0; i < 3; ++i) {
                neigh[i] = (neighbor[i] + grid.getDim()[i]) % grid.getDim()[i];
            }

            // the mirror we are adding so that the particle gets mirrored
            std::array<double, 3> mirror{};
            mirror = {
                    neighbor[0] == -1 ? -grid.getLenDim()[0] : neighbor[0] == grid.getDim()[0] && grid.getDim()[0] != 1
                                                               ? grid.getLenDim()[0] : 0.0,

                    neighbor[1] == -1 ? -grid.getLenDim()[1] : neighbor[1] == grid.getDim()[1] && grid.getDim()[1] != 1
                                                               ? grid.getLenDim()[1] : 0.0,

                    neighbor[2] == -1 ? -grid.getLenDim()[2] : neighbor[2] == grid.getDim()[2] && grid.getDim()[2] != 1
                                                               ? grid.getLenDim()[2] : 0.0};

            CellCalculator::calcPerNeighbors(grid, neigh, p, mirror, newton);
        }
    }

    void CellCalculator::calcFCellSubdomain(CellContainerIterators &curCell,
                                            CellContainer &grid,
                                            const CellContainer::SubDomain &subDomain) {
        for (auto &p: curCell) {
            // get all the neighbors
            for (const std::array<int, 3> &neighbors: grid.is2D() ? curCell.getNeighbors2D()
                                                                  : curCell.getNeighbors3D()) {
                // The neighbor is in the subdomain
                if (subDomain.isInSubdomain(neighbors)) {
                    calcFWithNeighbors(grid, p, neighbors, true);
                } else {
                    calcFWithNeighbors(grid, p, neighbors, false);
                }
            }

            for (const std::array<int, 3> &neighbors: grid.is2D() ? curCell.getRemainingNeighbors2D()
                                                                  : curCell.getRemainingNeighbors3D()) {
                if (!subDomain.isInSubdomain(neighbors)) {
                    calcFWithNeighbors(grid, p, neighbors, false);
                }
            }
        }

        // Calculates the forces within a cell
        calcFWithinCell(curCell);
        auto currentIndexes = curCell.getIndex();
        // checks if it is a border cell, if yes also calculate border forces
        if (curCell.isBorderCell1()) {
            reflectiveBoundary(grid, currentIndexes);
        }
    }

    void CellCalculator::calcFCell(CellContainerIterators &curCell, CellContainer &grid) {
        // Loop through every particle in this cell
        for (auto &p: curCell) {
            // Get all the neighbors of this cell
            for (const std::array<int, 3> &neighbors:
                    (grid.is2D() ? curCell.getNeighbors2D() : curCell.getNeighbors3D())) {
                calcFWithNeighbors(grid, p, neighbors, true);
            }
        }

        // Calculates the forces within a cell
        calcFWithinCell(curCell);
        auto currentIndexes = curCell.getIndex();
        // checks if it is a border cell, if yes also calculate border forces
        if (curCell.isBorderCell1()) {
            reflectiveBoundary(grid, currentIndexes);
        }
    }

    void CellCalculator::calcF(ParticleContainer &container) {
        auto &grid = static_cast<CellContainer &>(container);
        // Switches between parallelization technique
        switch (grid.getStrategy()) {
            case CellContainer::primitiveX:
            case CellContainer::primitiveY:
            case CellContainer::primitiveZ:
#ifdef _OPENMP
#pragma omp parallel shared(grid) default(none)
#endif
            {
#ifdef _OPENMP
#pragma omp for schedule(dynamic)
#endif
                // First part
                for (size_t i = 0; i < grid.getIndicesThreadVector().size(); ++i) {
                    for (int pos: grid.getIndicesThreadVector()[i]) {
                        if (!grid.grid[pos].getParticles().empty()) {
                            calcFCell(grid.grid[pos], grid);
                        }
                    }
                }
#ifdef _OPENMP
#pragma omp for schedule(dynamic) nowait // no wait since next line doesn't have to wait
#endif
                // Second part
                for (size_t i = 0; i < grid.getIndicesThreadVector().size(); ++i) {
                    for (int pos: grid.getIndicesThreadVector()[i]) {
                        if (!grid.grid[pos + grid.getThreadOffset()].getParticles().empty()) {
                            calcFCell(grid.grid[pos + grid.getThreadOffset()], grid);
                        }
                    }
                }
#ifdef _OPENMP
#pragma omp single // single last remaining line
#endif
                {
                    for (int pos: grid.getResidualThreadVector()) {
                        if (!grid.grid[pos].getParticles().empty()) {
                            calcFCell(grid.grid[pos], grid);
                        }
                    }
                };
            }; // end parallelization area
                break;

            case CellContainer::subDomain:
#ifdef _OPENMP
#pragma omp parallel shared(grid) default(none)
#endif
            {
#ifdef _OPENMP
#pragma omp for schedule(dynamic)
#endif
                for (auto &subDomain: grid.getSubDomainVector()) {
                    // Normal cells
                    for (int pos: subDomain.getCellIndices()) {
                        if (!grid.grid[pos].getParticles().empty()) {
                            calcFCell(grid.grid[pos], grid);
                        }
                    }
                    // Border cells
                    for (int pos: subDomain.getBorderCellIndices()) {
                        if (!grid.grid[pos].getParticles().empty()) {
                            calcFCellSubdomain(grid.grid[pos], grid, subDomain);
                        }
                    }
                }
            }; // end parallelization area
                break;

            case CellContainer::serial:
                for (auto &curCell: grid.grid) {
                    calcFCell(curCell, grid);
                }
                break;
        }
    }

    void CellCalculator::harmonic_potential(Particle *p1, Particle *p2, double sqrd_dist, bool newton) const {
        double factor = SQR_ROOT_OF_TWO;
        double distance = sqrt(sqrd_dist);
        // if orthogonal neighbors
        if (abs(p1->getGridIndex()[0] - p2->getGridIndex()[0]) + abs(p1->getGridIndex()[1] - p2->getGridIndex()[1])
            + abs(p1->getGridIndex()[2] - p2->getGridIndex()[2]) == 1) {
            factor = 1;
        }
        auto f = membraneForceParameter * (distance - factor * rZero) * (1. / distance) * (p2->getX() - p1->getX());
        p1->setF(p1->getF() + f);
        if (newton) {
            p2->setF(p2->getF() - f);
        }
    }


    void CellCalculator::reflectiveBoundary(CellContainer &grid, const std::array<int, 3> &currentIndex) const {
        // saves the reflective borders of the current cell
        std::vector<int> reflBorder{};
        // go through all three or two axis and acquire the borders of currentIndex that are reflective
        for (int d = 0; d < (grid.is2D() ? 2 : 3); ++d) {
            const auto [bordType, bord] = grid.getBorders(currentIndex, d);
            if (bordType == CellContainer::reflective) {
                reflBorder.push_back(bord);
            }
        }

        /**
         * Border type overview:
         * 0: LEFT
         * 1: RIGHT
         * 2: UPPER
         * 3: LOWER
         * 4: FRONT
         * 5: BACK
         */
        if (!reflBorder.empty()) {
            for (auto &p: grid.grid[grid.index(currentIndex)]) {

                for (int bord: reflBorder) {
                    double r = grid.getDistance(p->getX(), bord);
                    // reflect distance depending on different sigma for each particle;
                    double reflectDistance = SIXTH_ROOT_OF_TWO * sigma;
                    if (r <= reflectDistance) {
                        double s = (sigma * sigma) / (r * r);
                        s = s * s * s;
                        auto force = -24 * epsilon / r * (s - 2 * s * s);

                        auto newF{p->getF()};
                        switch (bord) {
                            case 0:
                                newF[0] += force;
                                break;
                            case 1:
                                newF[0] -= force;
                                break;
                            case 2:
                                newF[1] += force;
                                break;
                            case 3:
                                newF[1] -= force;
                                break;
                            case 4:
                                newF[2] += force;
                                break;
                            case 5:
                                newF[2] -= force;
                                break;
                            default:
                                spdlog::critical("DEFAULT CASE SOMETHING WRONG ALARM");
                        }
                        p->setF(newF);
                    }
                }
            }
        }
    }

    bool CellCalculator::gridNeighbors(const std::array<int, 3> &index1, const std::array<int, 3> &index2) {
        return abs(index1[0] - index2[0]) <= 1 && abs(index1[1] - index2[1]) <= 1 && abs(index1[2] - index2[2]) <= 1;
    }


    void CellCalculator::setRZero(double rZ) {
        CellCalculator::rZero = rZ;
    }

    void CellCalculator::setMembraneForceParameter(double kFP) {
        CellCalculator::membraneForceParameter = kFP;
    }

    void CellCalculator::setMembrane(bool mem) {
        CellCalculator::membrane = mem;
    }

    std::string CellCalculator::toString() {
        return "CellCalculator";
    }
}

}