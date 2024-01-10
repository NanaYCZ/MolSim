

#include <iostream>
#include "CellContainerIterators.h"
#include <stdexcept>

void CellContainer::CustomIterator::next_correct_index_custom(dim_t& x, dim_t& y, dim_t& z){
  if(x ==-1){
    throw new std::invalid_argument("Trying to increment an end pointer");
  }
  if(x < upper_x){
    x++;
  }else if(y < upper_y){
    x=lower_x;
    y++;
  }else if(z < upper_z){
    x=lower_x;
    y=lower_y;
    z++;
  }else{
    x = -1;
    y = 1;
    z = 1;
  }
}


CellContainer::CustomIterator &CellContainer::CustomIterator::operator++() {
  //std::cout << "Domain max: [" << cell.domain_max_dim[0] << " " <<  cell.domain_max_dim[1] << " " <<  cell.domain_max_dim[2] << "]";
  next_correct_index_custom(x,y,z);
  return *this;
}


std::vector<Particle*>& CellContainer::CustomIterator::operator*(){
  return cell.particles[x][y][z];
}

bool CellContainer::CustomIterator::operator==(const CellContainer::CustomIterator& other){
    return (x == other.x && y == other.y && z == other.z);
}

bool CellContainer::CustomIterator::operator!=(const CellContainer::CustomIterator& other){
    return (!(x == other.x && y == other.y && z == other.z));
}





void CellContainer::Iterator::next_correct_index(dim_t& x, dim_t& y, dim_t& z){
  if(x ==-1){
    throw new std::invalid_argument("Trying to increment an end pointer");
  }
  if(x < cell.domain_max_dim[0]){
    x++;
  }else if(y < cell.domain_max_dim[1]){
    x=1;
    y++;
  }else if(z < cell.domain_max_dim[2]){
    x=1;
    y=1;
    z++;
  }else{
    x = -1;
    y = 1;
    z = 1;
  }

}


CellContainer::Iterator &CellContainer::Iterator::operator++() {
  //std::cout << "Domain max: [" << cell.domain_max_dim[0] << " " <<  cell.domain_max_dim[1] << " " <<  cell.domain_max_dim[2] << "]";
  next_correct_index(x,y,z);
  return *this;
}


std::vector<Particle*>& CellContainer::Iterator::operator*(){
  return cell.particles[x][y][z];
}

bool CellContainer::Iterator::operator==(const CellContainer::Iterator& other){
    return (x == other.x && y == other.y && z == other.z);
}

bool CellContainer::Iterator::operator!=(const CellContainer::Iterator& other){
    return (!(x == other.x && y == other.y && z == other.z));
}


CellIterator::CellIterator(dim_t x, dim_t y, dim_t z) : x(x), y(y), z(z) {
    total = CellContainer::domain_max_dim[0] * CellContainer::domain_max_dim[1] * CellContainer::domain_max_dim[2];
}

CellIterator &CellIterator::operator+=(int p) {
    for (int i = 0; i < p; ++i) {
        next_index();
    }
    return *this;
}

CellIterator &CellIterator::operator++() {
    next_index();
    return *this;
}

std::vector<Particle*>& CellIterator::operator*() const {
    std::vector<Particle*>& ret = CellContainer::particles[x][y][z];
    return ret;
}

bool CellIterator::operator!=(CellIterator other) const {
    return !(x == other.x && y == other.y && z == other.z);
}

void CellIterator::next_index() {
    if(x == -1) {
        throw std::runtime_error("too many iterations");
    }

    --total;

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
    return a.total;
}

CellIterator begin_CI() {
    return CellIterator{1,1,1};
}

CellIterator end_CI() {
    return CellIterator{-1,-1,-1};
}


PeriodIterator::PeriodIterator(std::array<dim_t, 3> pattern) : progress(pattern){
    next_index();
    total = std::abs(pattern[0]) * CellContainer::domain_max_dim[1] * CellContainer::domain_max_dim[2] +
            std::abs(pattern[1]) * (CellContainer::domain_max_dim[0] - std::abs(pattern[0])) *
            CellContainer::domain_max_dim[2] +
            std::abs(pattern[2]) * (CellContainer::domain_max_dim[0] - std::abs(pattern[0])) *
            (CellContainer::domain_max_dim[1] - std::abs(pattern[1]));
}

PeriodIterator::PeriodIterator() : progress({0,0,0}), plane_axis(finished), total(0) {}

PeriodIterator &PeriodIterator::operator+=(int p) {
    for (int i = 0; i < p; ++i) {
        next_index();
    }
    return *this;
}

PeriodIterator PeriodIterator::operator++() {
    next_index();
    return *this;
}

bool PeriodIterator::operator!=(PeriodIterator other) {
    return progress[0] != other.progress[0] ||
           progress[1] != other.progress[1] ||
           progress[2] != other.progress[2] ||
           plane_axis != other.plane_axis;
}

std::array<dim_t, 3> PeriodIterator::operator*() {
    return {current[0] + mapping[0], current[1] + mapping[1], current[2] + mapping[2]};
}

//"iterates over every point of every plane"
void PeriodIterator::next_index() {
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

    --total;
}

//"iterates over every point of a plane"
void PeriodIterator::next_point_on_plane(short a, short b) {
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
void PeriodIterator::next_plane_corner() {
    if (progress[0] != 0) {
        plane_axis = x_axis;

        set_plane_position_on_axis();
        current[1] = min[1];
        current[2] = min[2];

        mapping[1] = progress[1];
        mapping[2] = progress[2];

    } else if (progress[1] != 0) {
        plane_axis = y_axis;

        current[0] = min[0];
        set_plane_position_on_axis();
        current[2] = min[2];

        mapping[0] = 0;
        mapping[2] = progress[2];

    } else if (progress[2] != 0) {
        plane_axis = z_axis;

        current[0] = min[0];
        current[1] = min[1];
        set_plane_position_on_axis();

        mapping[0] = 0;
        mapping[1] = 0;

    } else {
        //indicator that all points got iterated
        plane_axis = finished;
    }
}

void PeriodIterator::set_plane_position_on_axis() {
    if (progress[plane_axis] < 0) {
        //position of the next plane on axis
        current[plane_axis] = CellContainer::domain_max_dim[plane_axis] + progress[plane_axis] + 1;
        //mapping to receive domain exit point from current
        mapping[plane_axis] = -CellContainer::domain_max_dim[plane_axis];
        //update limit to avoid overlap with other planes
        max[plane_axis] = std::min(max[plane_axis], current[plane_axis] - 1);
        //update progress
        ++progress[plane_axis];

    } else {
        current[plane_axis] = progress[plane_axis];
        mapping[plane_axis] = CellContainer::domain_max_dim[plane_axis];
        min[plane_axis] = std::max(min[plane_axis], current[plane_axis] + 1);
        --progress[plane_axis];
    }
}


int operator-(PeriodIterator a, PeriodIterator b) { return b.total; }

PeriodIterator begin_PI(std::array<dim_t,3> pattern) {
    return {pattern};
}

PeriodIterator end_PI() {
    return {};
}
