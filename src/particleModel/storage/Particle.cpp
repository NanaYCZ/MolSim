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
}

Particle::Particle(const Particle &other) {
  x = other.x;
  v = other.v;
  m = other.m;
  type = other.type;
  f_1 = other.getF();
  f_2 = other.getOldF();
  SPDLOG_TRACE("Particle generated by copy!");
}


Particle::Particle(std::array<double, 3> x_arg, std::array<double, 3> v_arg,
                   double m_arg, int type_arg)
    : x(x_arg),
      v(v_arg),
      m(m_arg),
      type(type_arg),
      f_1({0., 0., 0.}),
      f_2({0., 0., 0.}) {
    SPDLOG_TRACE("Particle generated!");
}

Particle::~Particle() {
    SPDLOG_TRACE("Particle destructed!");
}

const std::array<double, 3> &Particle::getX() const { return x; }

const std::array<double, 3> &Particle::getV() const { return v; }

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

void Particle::setX(int index, double value) { x[index] = value; }

void Particle::setV(int index, double value) { v[index] = value; }

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

void Particle::shiftF() {
  if (secondIsOld) {
    secondIsOld = false;
    f_2 = {0.0, 0.0, 0.0};
  } else {
    secondIsOld = true;
    f_1 = {0.0, 0.0, 0.0};
  }
}

std::string Particle::toString() const {
  std::stringstream stream;
  stream << "Particle: X:" << x << " v: " << v << " f: " << getF()
         << " old_f: " << getOldF() << " type: " << type << " mass: " << m;
  return stream.str();
}

bool Particle::operator==(Particle &other) {
  return (x == other.x) and (v == other.v) and (getF() == other.getF()) and
         (type == other.type) and (m == other.m) and
         (getOldF() == other.getOldF());
}

bool Particle::operator==(const Particle& other) const
{
    if (this == &other) return true;
    else return false;
}

std::ostream &operator<<(std::ostream &stream, Particle &p) {
  stream << p.toString();
  return stream;
}
