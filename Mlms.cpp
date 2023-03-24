/* Copyright 2023 <Ilia Fedotov @ Uni Freiburg> */
#include <fstream>
#include <iostream>
#include <cmath>
#include <vector>
#include <array>
#include "./Mlms.h"

const double PI = 3.141592653589793238463;

// __________________________________________________________________
void printarray(const matrix &array) {
  for (std::size_t i = 0; i < array.shape[0]; i++) {
    for (std::size_t j = 0; j < array.shape[1]; j++) {
      printf("%.1f\t", array(i, j));
    }
    printf("\n");
  }
  return;
}

// __________________________________________________________________
void initializePressureArray(matrix &Pa, double lower_b,
                             double upper_b, double pressure) {
  for (std::size_t i = 0; i < Pa.shape[0]; i++) {
    for (std::size_t j = 0; j < Pa.shape[1]; j++) {
      if (j >= lower_b && i >= lower_b && j < upper_b && i < upper_b) {
        Pa(i, j) = 1;
      } else {
        Pa(i, j) = 0;
      }
    }
  }
}

// __________________________________________________________________
void initializeDisplacementArray(matrix &Ic) {
  for (std::size_t i = 0; i < Ic.shape[0]; i++) {
    for (std::size_t j = 0; j < Ic.shape[1]; j++) {
      Ic(i, j) = 0;
    }
  }
}
// __________________________________________________________________
inline double calculate(double a, double b, double x, double y) {
  double first, second, third, fourth;
  // naive formula implementation for calculating a pressure patch
  first = ((y + b) + sqrt(pow((y + b), 2) + pow((x + a), 2))) /
          ((y - b) + sqrt(pow((y - b), 2) + pow((x + a), 2)));
  first = (x + a) * log(first);

  second = ((x + a) + sqrt(pow((y + b), 2) + pow((x + a), 2))) /
           ((x - a) + sqrt(pow((y + b), 2) + pow((x - a), 2)));
  second = (y + b) * log(second);

  third = ((y - b) + sqrt(pow((y - b), 2) + pow((x - a), 2))) /
          ((y + b) + sqrt(pow((y + b), 2) + pow((x - a), 2)));
  third = (x - a) * log(third);

  fourth = ((x - a) + sqrt(pow((y - b), 2) + pow((x - a), 2))) /
            ((x + a) + sqrt(pow((y - b), 2) + pow((x + a), 2)));
  fourth = (y - b) * log(fourth);
  double res = (first+second+third+fourth);
  return res;
}

/* // __________________________________________________________________
inline double calculate(double a, double b, double x, double y) {
  double first, second, third, fourth;
  double y_min_b = pow((y - b), 2);
  double y_plus_b = pow((y + b), 2);
  double x_min_a = pow((x - a), 2);
  double x_plus_a = pow((x + a), 2);
  // naive formula implementation for calculating a pressure patch
  first = ((y + b) + sqrt(y_plus_b + x_plus_a)) /
          ((y - b) + sqrt(y_min_b + x_plus_a));
  first = (x + a) * log(first);

  second = ((x + a) + sqrt(y_plus_b + x_plus_a)) /
           ((x - a) + sqrt(y_plus_b + x_min_a));
  second = (y + b) * log(second);

  third = ((y - b) + sqrt(y_min_b + x_min_a)) /
          ((y + b) + sqrt(y_plus_b + x_min_a));
  third = (x - a) * log(third);

  fourth = ((x - a) + sqrt(y_min_b + x_min_a)) /
            ((x + a) + sqrt(y_min_b + x_plus_a));
  fourth = (y - b) * log(fourth);
  double res = (first+second+third+fourth);
  return res;
} */


// __________________________________________________________________
void calculation_loop(matrix &Ic, const matrix &Pa,
                double cell_size,
                double v, double E) {
  // outer loop over grid to call displacement calculation
  // for parallelisation uncomment, compiler option needed: -fopen
  // #pragma omp parallel for
  for (int i = 0; i < static_cast<int>(Ic.shape[0]); i++) {
    for (int j = 0; j < static_cast<int>(Ic.shape[1]); j++) {
      Ic(i, j) = calc_displacement(Pa, Ic, i, j, cell_size/2,
                          cell_size/2, v, E, cell_size);
    }
  }
}

// __________________________________________________________________
double calc_displacement(const matrix &pressure,
                const matrix &Ic,
                double y, double x,
                double a, double b , double v,
                double E, double cell) {
  double res = 0;

  for (int i = 0; i < static_cast<int>(Ic.shape[0]); i+=1) {
    for (int j = 0; j < static_cast<int>(Ic.shape[1]); j+=1) {
      res += calculate(a, b, (x-j)*cell, (y-i)*cell) *
                            ((1-v)/(PI*E)) * pressure(i, j);
    }
  }
  return res;
}
