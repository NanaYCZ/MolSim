#include "CellContainerIterators.h"
#include <stdexcept>

CellIterator::CellIterator(dim_t x, dim_t y, dim_t z) : x(x), y(y), z(z) {
    remaining = CellContainer::domain_max_dim[0] * CellContainer::domain_max_dim[1] * CellContainer::domain_max_dim[2];
}

CellIterator &CellIterator::operator+=(int amount) {
    for (int i = 0; i < amount; ++i) {
        next_index();
    }
    return *this;
}

CellIterator &CellIterator::operator++() {
    next_index();
    return *this;
}

std::array<dim_t,3> CellIterator::position() {
    return {x,y,z};
}

std::vector<Particle*>& CellIterator::operator*() const {
    return CellContainer::particles[x][y][z];
}

bool CellIterator::operator!=(CellIterator other) const {
    return !(x == other.x && y == other.y && z == other.z);
}

void CellIterator::next_index() {
    if(x == -1) {
        throw std::runtime_error("too many iterations");
    }

    --remaining;

    if (x < CellContainer::domain_max_dim[0]) {
        x++;
    } else if (y < CellContainer::domain_max_dim[1]) {
        x = 1;
        y++;
    } else if (z < CellContainer::domain_max_dim[2]) {
        x = 1;
        y = 1;
        z++;
    } else {
        x=-1;
        y=-1;
        z=-1;
    }
}

int operator-(CellIterator a, CellIterator b) {
    return a.remaining;
}

CellIterator begin_CellIterator() {
    return CellIterator{1,1,1};
}

CellIterator end_CellIterator() {
    return CellIterator{-1,-1,-1};
}


StartPointIterator::StartPointIterator(std::array<dim_t, 3> pattern) : progress(pattern){
    next_index();
    remaining = std::abs(pattern[0]) * CellContainer::domain_max_dim[1] * CellContainer::domain_max_dim[2] +
            std::abs(pattern[1]) * (CellContainer::domain_max_dim[0] - std::abs(pattern[0])) *
            CellContainer::domain_max_dim[2] +
            std::abs(pattern[2]) * (CellContainer::domain_max_dim[0] - std::abs(pattern[0])) *
            (CellContainer::domain_max_dim[1] - std::abs(pattern[1]));
}

StartPointIterator::StartPointIterator() : progress({0, 0, 0}), plane_axis(finished), remaining(0) {}

StartPointIterator &StartPointIterator::operator+=(int p) {
    for (int i = 0; i < p; ++i) {
        next_index();
    }
    return *this;
}

StartPointIterator StartPointIterator::operator++() {
    next_index();
    return *this;
}

bool StartPointIterator::operator!=(StartPointIterator other) {
    return progress[0] != other.progress[0] ||
           progress[1] != other.progress[1] ||
           progress[2] != other.progress[2] ||
           plane_axis != other.plane_axis;
}

std::array<dim_t, 3> StartPointIterator::outside() {
    return {current[0] + mapping[0], current[1] + mapping[1], current[2] + mapping[2]};
}

std::array<dim_t, 3> StartPointIterator::operator*() {
    return current;
}

//"iterates over every point of every plane"
void StartPointIterator::next_index() {
    switch (plane_axis) {
        case x_axis: {
            next_point_on_plane(1, 2);
            break;
        }
        case y_axis: {
            next_point_on_plane(0, 2);
            break;
        }
        case z_axis: {
            next_point_on_plane(0, 1);
            break;
        }
        case reset: {
            next_plane_corner();
            break;
        }
        default:
            throw std::runtime_error("too many iterations");
    }

    --remaining;
}

//"iterates over every point of a plane"
void StartPointIterator::next_point_on_plane(short a, short b) {
    if (current[a] < max[a]) {
        ++current[a];

    } else if (current[b] < max[b]) {
        current[a] = min[a];
        ++current[b];
    }

    if (current[a] == max[a] && current[b] == max[b]) {
        //indicator that last point on current plane is reached
        plane_axis = reset;
    }
}

//"iterates over every plane"
void StartPointIterator::next_plane_corner() {
    if (progress[0] != 0) {
        plane_axis = x_axis;

        set_axis();
        current[1] = min[1];
        current[2] = min[2];

        mapping[1] = progress[1];
        mapping[2] = progress[2];

    } else if (progress[1] != 0) {
        plane_axis = y_axis;

        current[0] = min[0];
        set_axis();
        current[2] = min[2];

        mapping[0] = 0;
        mapping[2] = progress[2];

    } else if (progress[2] != 0) {
        plane_axis = z_axis;

        current[0] = min[0];
        current[1] = min[1];
        set_axis();

        mapping[0] = 0;
        mapping[1] = 0;

    } else {
        //indicator that all points got iterated
        plane_axis = finished;
    }
}

void StartPointIterator::set_axis() {
    if (progress[plane_axis] < 0) {
        //position of the next plane
        current[plane_axis] = CellContainer::domain_max_dim[plane_axis] + progress[plane_axis] + 1;
        //mapping to receive domain exit point from current
        mapping[plane_axis] = -CellContainer::domain_max_dim[plane_axis];
        //update limit to avoid overlap with other planes
        max[plane_axis] = std::min(max[plane_axis], current[plane_axis] - 1);
        //update progress
        ++progress[plane_axis];

    } else {
        //position of the next plane
        current[plane_axis] = progress[plane_axis];
        //mapping to receive domain exit point from current
        mapping[plane_axis] = CellContainer::domain_max_dim[plane_axis];
        //update limit to avoid overlap with other planes
        min[plane_axis] = std::max(min[plane_axis], current[plane_axis] + 1);
        //update progress
        --progress[plane_axis];
    }
}


int operator-(StartPointIterator a, StartPointIterator b) { return b.remaining; }

StartPointIterator begin_StartIterator(std::array<dim_t,3> pattern) {
    return {pattern};
}

StartPointIterator end_StartIterator() {
    return {};
}
