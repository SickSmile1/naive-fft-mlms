#ifndef BOUSSINESQ1_HH
#define BOUSSINESQ1_HH

#include "BoussinesqMatrix1.hh"
#include <cmath>

inline double boussinesq(double x, double y, double a, double b) {
  double plus_sqrt = std::sqrt((x + a) * (x + a) + (y + b) * (y + b));
  double minu_sqrt = std::sqrt((x - a) * (x - a) + (y - b) * (y - b));
  double pm_sqrt = std::sqrt((x + a) * (x + a) + (y - b) * (y - b));
  double mp_sqrt = std::sqrt((x - a) * (x - a) + (y + b) * (y + b));
  return 1. / M_PI *
         ((x + a) * std::log((y + b + plus_sqrt) / (y - b + pm_sqrt)) +
          (y + b) * std::log((x + a + plus_sqrt) / (x - a + mp_sqrt)) +
          (x - a) * std::log((y - b + minu_sqrt) / (y + b + mp_sqrt)) +
          (y - b) * std::log((x - a + minu_sqrt) / (x + a + pm_sqrt)));
}

inline double boussinesq(int i, int j, const vec2d& metric, const vec2d& pixel) {
  return boussinesq(i * metric[0], j * metric[1], pixel[0] / 2, pixel[1] / 2);
}

double deflection(const matrix& pressure, vec2d source, vec2d metric, vec2d pixel);

void naive(const matrix& pressure, matrix& displacement, vec2d metric, vec2d pixel);
  

#endif // BOUSSINESQ1_HH
