/*
 * Particle.cpp
 *
 *  Created on: 23.02.2010
 *      Author: eckhardw
 */

#define SPDLOG_ACTIVE_LEVEL SPDLOG_LEVEL_TRACE

#include "Particle.h"
#include "utils/ArrayUtils.h"

#include <iostream>
#include "spdlog/spdlog.h"

Particle::Particle(int type_arg) {
  type = type_arg;
  SPDLOG_TRACE("Particle generated!");
  f_1 = {0., 0., 0.};
  f_2 = {0., 0., 0.};
  boundaries_crossed =  {0,0,0};
}

Particle::Particle(const Particle &other) {
  x = other.x;
  v = other.v;
  m = other.m;
  type = other.type;
  f_1 = other.getF();
  f_2 = other.getOldF();
  boundaries_crossed = other.boundaries_crossed;
  grid = other.getGrid();
  RZ = other.RZ;
  FP = other.FP;
  special = other.getSpecial();
  SPDLOG_TRACE("Particle generated by copy!");
}


Particle::Particle(std::array<double, 3> x_arg, std::array<double, 3> v_arg,
                   double m_arg, int type_arg)
    : x(x_arg),
      v(v_arg),
      m(m_arg),
      type(type_arg),
      f_1({0., 0., 0.}),
      f_2({0., 0., 0.}),
      boundaries_crossed({0,0,0})
//      RZ(),
//      FP(),
//      special(),
//      grid()
      {
    SPDLOG_TRACE("Particle generated!");
}


Particle::Particle(std::array<double, 3> x_arg, std::array<double, 3> v_arg, double m_arg, int type_arg,
                   std::array<int, 3> grid, double rz, double fp, std::array<double, 3> special)
        : x(x_arg),
          v(v_arg),
          m(m_arg),
          type(type_arg),
          f_1({0., 0., 0.}),
          f_2({0., 0., 0.}),
          boundaries_crossed({0,0,0}),
          grid(grid),
          RZ(rz),
          FP(fp),
          special(special)

{
    SPDLOG_TRACE("Particle generated!");
}


Particle::~Particle() {
    SPDLOG_TRACE("Particle destructed!");
}

const std::array<double, 3> &Particle::getX() const { return x; }

const std::array<double, 3> &Particle::getV() const { return v; }

const std::array<int, 3> &Particle::getGrid()  const { return grid; }

const std::array<double, 3> &Particle::getF() const {
  if (secondIsOld) {
    return f_1;
  }
  return f_2;
}

const std::array<double, 3> &Particle::getOldF() const {
  if (secondIsOld) {
    return f_2;
  }
  return f_1;
}

double Particle::getM() const { return m; }

int Particle::getType() const { return type; }

double Particle::getRZ() const { return RZ; }

double Particle::getFP() const { return FP; }

std::array<double, 3> Particle::getSpecial() const {return special;}

void Particle::setX(int index, double value) { x[index] = value; }

void Particle::setX(std::array<double,3> new_x) {x = new_x;};

void Particle::setV(int index, double value) { v[index] = value; }

void Particle::setV(std::array<double,3> new_v){v = new_v;};

void Particle::setGrid(std::array<int,3> index){grid = index;};

void Particle::setRZ(double value) { RZ = value; }

void Particle::setFP(double value) { FP = value; }

void Particle::setSpecial(std::array<double, 3> value) { special = value; }


void Particle::addX(int index, double value){
  x[index] += value;
}

void Particle::addX(std::array<double,3> &x_add) {
    x[0] += x_add[0];
    x[1] += x_add[1];
    x[2] += x_add[2];
}

void Particle::addF(int index, double value) {
  if (secondIsOld) {
    f_1[index] += value;
  } else {
    f_2[index] += value;
  }
}

void Particle::addF(std::array<double,3> add_f){
  if (secondIsOld) {
    f_1 = f_1 + add_f;
  } else {
    f_2 = f_2 + add_f;
  }
}

void Particle::shiftF() {
  if (secondIsOld) {
    secondIsOld = false;
    f_2 = {0.0, 0.0, 0.0};
  } else {
    secondIsOld = true;
    f_1 = {0.0, 0.0, 0.0};
  }
}

std::array<int,3>& Particle::getBoundariesCrossed(){
    return boundaries_crossed;
}

int& Particle::getBoundariesCrossed(int i){
    return boundaries_crossed[i];
}

void  Particle::incBoundariesCrossedI(int i){
    boundaries_crossed[i]+=1;
}

void Particle::decBoundariesCrossedI(int i){
    boundaries_crossed[i]-=1;
};

void Particle::setBoundariesCrossedZero(){
    boundaries_crossed = {0, 0, 0};
}

std::string Particle::toString() const {
  std::stringstream stream;
  stream << "Particle: X:" << x << " v: " << v << " f: " << getF() << " crossed: " << boundaries_crossed 
         << " old_f: " << getOldF() << " type: " << type << " mass: " << m;
  return stream.str();
}

bool Particle::operator==(Particle &other) {
  auto arr_eq = [](std::array<double,3> arr1,std::array<double,3> arr2) -> bool{
    for(int i = 0; i < 3 ; i++){
      //allow for error
      if (std::fabs(arr1[i] - arr2[i]) > 1e-4) {
              return false; // Arrays are not close
      }
    }
    return true;
  };

  return arr_eq(x,other.x) and arr_eq(v,other.v) and arr_eq(getF(),other.getF()) and
         (type == other.type) and (m == other.m) and
         arr_eq(getOldF(),other.getOldF());
}

bool Particle::operator==(const Particle& other) const
{
    if (this == &other) return true;
    else return false;
}

