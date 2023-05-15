/* Copyright 2023 <Ilia Fedotov @ Uni Freiburg> */
#include <fftw3.h>
#include <fstream>
#include <iostream>
#include <cmath>
// #include <vector>
// #include <complex>
#include "./Boussinesq.h"

// __________________________________________________________________
void initializePressureArray(matrix &Pa, double lower_b, // NOLINT
                             double upper_b, double pressure) {
  initializeDisplacementArray(Pa);
  for (int i = lower_b; i < upper_b; i++) {
    for (int j = lower_b; j < upper_b; j++) {
      Pa(i, j) = pressure;
    }
  }
}

// __________________________________________________________________
void initializeDisplacementArray(matrix &Ic) { // NOLINT
  for (int i = 0; i < Ic.shape[0]; i++) {
    for (int j = 0; j < Ic.shape[1]; j++) {
      Ic(i, j) = 0;
    }
  }
}

// __________________________________________________________________
double calcBoussinesq(double a, double b, double x, double y) {
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
  double res = (first+second+third+fourth)/M_PI;
  return res;
}

// __________________________________________________________________
double calcBoussinesq(int i, int j, double dxc, double dyc,
    double dxf, double dyf) {
  return calcBoussinesq(dxf/2, dyf/2, i*dxc, j*dyc);
}
// __________________________________________________________________
void naiveCalculation(matrix &Ic, const matrix &Pa, // NOLINT
                double cell_size) {
  // outer loop over grid to call displacement calculation
  // for parallelisation uncomment, compiler option needed: -fopen
  #pragma omp parallel for simd
  for (int i = 0; i < Ic.shape[0]; i++) {
    for (int j = 0; j < Ic.shape[1]; j++) {
      Ic(i, j) = calc_displacement(Pa, Ic, i, j, cell_size/2,
                          cell_size/2, cell_size);
    }
  }
}

// __________________________________________________________________
double calc_displacement(const matrix &pressure,
              const matrix &Ic,
              double y, double x,
              double a, double b, double cell) {
  double res = 0;

  for (int i = 0; i < Ic.shape[0]; i+=1) {
    for (int j = 0; j < Ic.shape[1]; j+=1) {
      double xj = (x*cell)-(j*cell);
      double yi = (y*cell)-(i*cell);
      res += calcBoussinesq(a, b, xj, yi) * pressure(i, j);
    }
  }
  return res;
}

void calc_displacement(const matrix &pF, double cS, double fS,
                         matrix &cD) { // NOLINT
  int shape = pF.shape[0];
  // #pragma omp parallel for simd
  for (int i = 0; i < shape; ++i) {
    for (int j = 0; j < shape; ++j) {
      // inner loop
      for (int k = 0; k < shape; ++k) {
        for (int l = 0; l < shape; ++l) {
          cD(i, j) += calcBoussinesq(i-k, j-l,  cS, cS,
                          fS, fS)*pF(k, l);
        }
      }
    }
  }
}

