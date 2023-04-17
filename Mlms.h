/* Copyright 2023 <Ilia Fedotov @ Uni Freiburg> */
#ifndef MLMS_H_
#define MLMS_H_

// #include <gtest/gtest.h>
#include <vector>
#include <array>
#include <cassert>

struct matrix{
  // data container for pressure arrays, grid and result
  std::array<std::size_t, 2> shape;
  std::vector<double> data;

  // constructor, to initialize simulated 2d array in a 1d shape
  explicit matrix(std::array<std::size_t, 2> shape_) : shape(shape_) {
    data.resize(shape[0]*shape[1]);
  }

  // overload () for aquiring (conversion) i,j idexes in 1d array
  double& operator() (std::size_t i, std::size_t j) {
    assert(i < shape[0] && j < shape[1]);
    return data[j+(i*shape[1])];
  }

  // overload () for aquiring (conversion) i,j indexes in 1d array
  // in static form
  const double& operator() (std::size_t i, std::size_t j) const {
    assert(i < shape[0] && j < shape[1]);
    return data[j+(i*shape[1])];
  }
};

// initialization of pressure array
void initializePressureArray(matrix &Pa, double lower_b,
                             double upper_b, double pressure);

// initialization with 0 for all elements
void initializeDisplacementArray(matrix &Ic);



// most outer loop, for we need a n*n loop for the algorithm
void calculation_loop(matrix &Ic, const matrix &Pa, double cell_size,
                 double v, double E);

void calculation_loop2(matrix &Ic, const matrix &Pa, double cell_size,
                 double v, double E);

// inner loop calling the calulation for every n
double calculate(double a, double b, double x, double y);

// actual calculation
double calc_displacement(const matrix &pressure,
                const matrix &Ic,
                double y, double x,
                double a, double b,
                double v, double E,
                double cell);

// printing function for the data container
void printarray(const matrix &array);

// coarse pressure array calculation
void calc_coarse_pressure(const matrix &fP, const matrix &st,
                          matrix &cP, std::size_t t);

// initialization with 1 for all elements
void initializeMultiplicationArray(matrix &Ic);

// calculate the sizes for the coarse grids
void prepareCoarseSizes(std::vector<std::size_t> &gridLen1,
                std::vector<std::size_t> &gridLen2,
                std::vector<std::size_t> &coarseGrid,
                const std::size_t shape,
                const std::size_t t);


void deflectionCorrection(matrix &kM, const double mc,
                const matrix &cC, const matrix &cP, const std::size_t t);

void calcCorrMatrix(matrix &correctionCoefficients, const double mc,
                const matrix &st, const std::size_t t,
                const double fineSizeA, const double fineSizeB);


void coarseToFineGrid(const matrix &cM, matrix &kM, const matrix &st,
                const std::size_t t);

void fineGridCorrection(matrix &kM, const matrix &st,
               const double mc, const std::size_t ts,
               const matrix fP);
void myBreakpoint();
#endif  // MLMS_H_
