#include <iostream>
#include "Boussinesq.hh"

double deflection(const matrix &pressure, vec2d source, vec2d metric, vec2d pixel) {
  double result = 0;
  double x = source[0], y = source[1];
  double a = pixel[0] / 2, b = pixel[1] / 2;
  double dx = metric[0], dy = metric[1];

  for (std::size_t i = 0; i < pressure.shape[0]; ++i) {
    for (std::size_t j = 0; j < pressure.shape[1]; ++j) {
      double x_p = dx * i, y_p = dy * j;
      auto u = boussinesq(x - x_p, y - y_p, a, b);
      result += pressure(i, j) * u;
    }
  }

  return result;
}

void naive(const matrix& pressure, matrix& displacement, vec2d metric, vec2d pixel) {
  // compute deflection
  for (std::size_t i = 0; i < displacement.shape[0]; ++i) {
    for (std::size_t j = 0; j < displacement.shape[1]; ++j) {
      displacement(i, j) = deflection(pressure, {i*metric[0], j*metric[1]}, metric, pixel);
    }
  }
}
