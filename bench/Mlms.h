/* Copyright 2023 <Ilia Fedotov @ Uni Freiburg> */
#ifndef MLMS_H_
#define MLMS_H_

// #include <gtest/gtest.h>
#include <vector>
#include <array>

struct matrix{
  std::array<std::size_t, 2> shape;
  std::vector<double> data;

  explicit matrix(std::array<std::size_t, 2> shape_) : shape(shape_) {
    data.resize(shape[0]*shape[1]);
  }
  double& operator() (std::size_t i, std::size_t j) {
    return data[j+(i*shape[1])];
  }

  const double& operator() (std::size_t i, std::size_t j) const {
    return data[j+(i*shape[1])];
  }
};

void initializePressureArray(matrix &Pa, double lower_b,
                             double upper_b, double pressure);

void initializeRandomPressureArray(matrix &Pa);

void initializeDisplacementArray(matrix &Ic);

void calculation_loop(matrix &Ic, const matrix &Pa, double cell_size,
                 double v, double E);

double calc_displacement(const matrix &pressure,
                const matrix &Ic,
                double y, double x,
                double a, double b,
                double v, double E,
                double cell);

double calculate(double a, double b, double x, double y);

void printarray(const matrix &array);

#endif  // MLMS_H_
