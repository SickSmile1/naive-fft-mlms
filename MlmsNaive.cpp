/* Copyright 2023 <Ilia Fedotov @ Uni Freiburg> */
#include <fstream>
#include <iostream>
#include <cmath>
#include <vector>
#include "./Mlms.h"

const double PI = 3.141592653589793238463;

// __________________________________________________________________
void printarray(const matrix &array) {
  for (std::size_t i = 0; i < array.shape[0]; i++) {
    for (std::size_t j = 0; j < array.shape[1]; j++) {
      printf("%.4f\t", array(i, j));
    }
    printf("\n");
  }
  return;
}

// __________________________________________________________________
void initializePressureArray(matrix &Pa, double lower_b,
                             double upper_b, double pressure) {
  initializeDisplacementArray(Pa);
  for (std::size_t i = lower_b; i < upper_b; i++) {
    for (std::size_t j = lower_b; j < upper_b; j++) {
      Pa(i, j) = pressure;
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


// __________________________________________________________________
void calculation_loop(matrix &Ic, const matrix &Pa,
                double cell_size,
                double v, double E) {
  // outer loop over grid to call displacement calculation
  // for parallelisation uncomment, compiler option needed: -fopen
  #pragma omp parallel for simd
  for (std::size_t i = 0; i < Ic.shape[0]; i++) {
    for (std::size_t j = 0; j < Ic.shape[1]; j++) {
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

  for (std::size_t i = 0; i < Ic.shape[0]; i+=1) {
    for (std::size_t j = 0; j < Ic.shape[1]; j+=1) {
      double xj = (x*cell)-(j*cell);
      double yi = (y*cell)-(i*cell);
      // std::cout << "xi: " << xj << " yi: " << yi << std::endl;
      res += calculate(a, b, xj, yi) *
                            ((1-v)/(PI*E)) * pressure(i, j);
    }
  }
  return res;
}
