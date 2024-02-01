#include <iostream>
#include "CellContainerIterators.h"

std::vector<Particle *>::iterator CellContainerIterators::begin() {
    return particles.begin();
}

std::vector<Particle *>::iterator CellContainerIterators::end() {
    return particles.end();
}

std::vector<Particle *>::const_iterator CellContainerIterators::begin() const {
    return particles.begin();
}

std::vector<Particle *>::const_iterator CellContainerIterators::end() const {
    return particles.end();
}

std::vector<Particle *>::iterator CellContainerIterators::erase(std::vector<Particle *>::const_iterator position) {
    return particles.erase(position);
}

void CellContainerIterators::add_particle(Particle &p) {
    Particle *pp = &p;
    particles.emplace_back(pp);
}


void CellContainerIterators::clear() {
    particles.clear();
}

const std::array<int, 3> &CellContainerIterators::getIndex() const {
    return index;
}

void CellContainerIterators::setIndex(const std::array<int, 3> &indexV) {
    CellContainerIterators::index = indexV;
}

CellContainerIterators::CellContainerIterators() : particles{std::vector<Particle *>{}} {}

void CellContainerIterators::setNeighbors2DInX() {
    CellContainerIterators::neighbors2D = {
            std::array<int, 3>{index[0] + 1, index[1], index[2]}, // right
            std::array<int, 3>{index[0], index[1] + 1, index[2]}, // down
            std::array<int, 3>{index[0] + 1, index[1] + 1, index[2]}, // down right
            std::array<int, 3>{index[0] + 1, index[1] - 1, index[2]} // up right
    };
}

void CellContainerIterators::setNeighbors3DInX() {
    CellContainerIterators::neighbors3D = {
            std::array<int, 3>{index[0], index[1] + 1, index[2]}, // down, center, center
            std::array<int, 3>{index[0] + 1, index[1] + 1, index[2]}, // down, right, center
            std::array<int, 3>{index[0] + 1, index[1], index[2]}, // center, right, center
            std::array<int, 3>{index[0] + 1, index[1] - 1, index[2]}, // up, right, center
            std::array<int, 3>{index[0], index[1] + 1, index[2] - 1}, // down, center, front
            std::array<int, 3>{index[0], index[1], index[2] - 1}, // front
            std::array<int, 3>{index[0], index[1] - 1, index[2] - 1}, // up, center, front
            std::array<int, 3>{index[0] + 1, index[1] - 1, index[2] + 1}, // up, right, back
            std::array<int, 3>{index[0] + 1, index[1], index[2] + 1}, // center, right, back
            std::array<int, 3>{index[0] + 1, index[1] + 1, index[2] + 1}, // down, right, back
            std::array<int, 3>{index[0] + 1, index[1] - 1, index[2] - 1}, // up, right, front
            std::array<int, 3>{index[0] + 1, index[1], index[2] - 1}, // right, front
            std::array<int, 3>{index[0] + 1, index[1] + 1, index[2] - 1}, // down, right, front
    };
}

void CellContainerIterators::setNeighbors2DInY() {
    CellContainerIterators::neighbors2D = {
            std::array<int, 3>{index[0], index[1] + 1, index[2]},
            std::array<int, 3>{index[0] + 1, index[1], index[2]},
            std::array<int, 3>{index[0] + 1, index[1] + 1, index[2]},
            std::array<int, 3>{index[0] - 1, index[1] + 1, index[2]}
    };
}

void CellContainerIterators::setNeighbors3DInY() {
    CellContainerIterators::neighbors3D = {
            std::array<int, 3>{index[0] + 1, index[1], index[2]}, // right
            std::array<int, 3>{index[0] + 1, index[1] + 1, index[2]}, // right, down
            std::array<int, 3>{index[0], index[1] + 1, index[2]}, // down
            std::array<int, 3>{index[0] - 1, index[1] + 1, index[2]}, // left down
            std::array<int, 3>{index[0] + 1, index[1], index[2] - 1}, // front right
            std::array<int, 3>{index[0], index[1], index[2] - 1}, // front
            std::array<int, 3>{index[0] - 1, index[1], index[2] - 1}, // front left
            std::array<int, 3>{index[0] - 1, index[1] + 1, index[2] + 1}, // down back left
            std::array<int, 3>{index[0], index[1] + 1, index[2] + 1}, // down, back
            std::array<int, 3>{index[0] + 1, index[1] + 1, index[2] + 1}, // right, down, back
            std::array<int, 3>{index[0] - 1, index[1] + 1, index[2] - 1}, // left, down, front
            std::array<int, 3>{index[0], index[1] + 1, index[2] - 1}, // down, front
            std::array<int, 3>{index[0] + 1, index[1] + 1, index[2] - 1}, // right, down, front
    };
}

void CellContainerIterators::setNeighbors3DInZ() {
    CellContainerIterators::neighbors3D = {
            std::array<int, 3>{index[0], index[1] + 1, index[2]},
            std::array<int, 3>{index[0] + 1, index[1], index[2]},
            std::array<int, 3>{index[0] + 1, index[1] + 1, index[2]},
            std::array<int, 3>{index[0] - 1, index[1] + 1, index[2]},
            std::array<int, 3>{index[0], index[1], index[2] + 1},
            std::array<int, 3>{index[0], index[1] + 1, index[2] + 1},
            std::array<int, 3>{index[0] + 1, index[1], index[2] + 1},
            std::array<int, 3>{index[0] + 1, index[1] + 1, index[2] + 1},
            std::array<int, 3>{index[0] - 1, index[1], index[2] + 1},
            std::array<int, 3>{index[0], index[1] - 1, index[2] + 1},
            std::array<int, 3>{index[0] - 1, index[1] - 1, index[2] + 1},
            std::array<int, 3>{index[0] - 1, index[1] + 1, index[2] + 1},
            std::array<int, 3>{index[0] + 1, index[1] - 1, index[2] + 1}
    };
}

const std::vector<std::array<int, 3>> &CellContainerIterators::getNeighbors2D() const {
    return neighbors2D;
}

void CellContainerIterators::setNeighbors2D1(const std::vector<std::array<int, 3>> &neighbors_2_d) {
    neighbors2D = neighbors_2_d;
}

const std::vector<std::array<int, 3>> &CellContainerIterators::getNeighbors3D() const {
    return neighbors3D;
}

void CellContainerIterators::setNeighbors3D1(const std::vector<std::array<int, 3>> &neighbors_3_d) {
    neighbors3D = neighbors_3_d;
}

bool CellContainerIterators::isBorderCell1() const {
    return isBorderCell;
}

void CellContainerIterators::setIsBorderCell(bool is_border_cell) {
    isBorderCell = is_border_cell;
}

void CellContainerIterators::setRemainingNeighbors2D() {
    CellContainerIterators::remainingNeighbors2D = {
            std::array<int, 3>{index[0], index[1] - 1, index[2]},
            std::array<int, 3>{index[0] - 1, index[1], index[2]},
            std::array<int, 3>{index[0] - 1, index[1] - 1, index[2]},
            std::array<int, 3>{index[0] + 1, index[1] - 1, index[2]},
    };
}

void CellContainerIterators::setRemainingNeighbors3D() {
    CellContainerIterators::remainingNeighbors3D = {
            std::array<int, 3>{index[0] - 1, index[1], index[2]},
            std::array<int, 3>{index[0] - 1, index[1] - 1, index[2]},
            std::array<int, 3>{index[0], index[1] - 1, index[2]},
            std::array<int, 3>{index[0] + 1, index[1] - 1, index[2]},
            std::array<int, 3>{index[0] - 1, index[1], index[2] + 1},
            std::array<int, 3>{index[0], index[1], index[2] + 1},
            std::array<int, 3>{index[0] + 1, index[1], index[2] + 1},
            std::array<int, 3>{index[0] + 1, index[1] - 1, index[2] - 1},
            std::array<int, 3>{index[0], index[1] - 1, index[2] - 1},
            std::array<int, 3>{index[0] - 1, index[1] - 1, index[2] - 1},
            std::array<int, 3>{index[0] + 1, index[1] - 1, index[2] + 1},
            std::array<int, 3>{index[0], index[1] - 1, index[2] + 1},
            std::array<int, 3>{index[0] - 1, index[1] - 1, index[2] + 1},
    };
}

const std::vector<std::array<int, 3>> &CellContainerIterators::getRemainingNeighbors2D() const {
    return remainingNeighbors2D;
}

const std::vector<std::array<int, 3>> &CellContainerIterators::getRemainingNeighbors3D() const {
    return remainingNeighbors3D;
}
