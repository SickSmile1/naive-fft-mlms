/* Copyright 2023 <Ilia Fedotov @ Uni Freiburg> */
#ifndef MLMS_H_
#define MLMS_H_

// #include <gtest/gtest.h>
#include <vector>
#include <array>
#include <cassert>

struct matrix_i{
  // data container for pressure arrays, grid and result
  std::array<int, 2> shape;
  std::vector<double> data;

  // constructor, to initialize simulated 2d array in a 1d shape
  explicit matrix_i(std::array<int, 2> shape_) : shape(shape_) {
    data.resize(shape[0]*shape[1]);
  }
  auto begin() {
    return data.begin();
  }

  auto end() {
    return data.end();
  }
  // overload () for aquiring (conversion) i,j idexes in 1d array
  double& operator() (int i, int j) {
    // assert(i < shape[0] && j < shape[1]);
    return data[j+(i*shape[1])];
  }

  // overload () for aquiring (conversion) i,j indexes in 1d array
  // in static form
  const double& operator() (int i, int j) const {
    // assert(i < shape[0] && j < shape[1]);
    return data[j+(i*shape[1])];
  }
};

// initialization of pressure array
void initializePressureArray(matrix_i &Pa, double lower_b, // NOLINT
                             double upper_b, double pressure); // NOLINT

// initialization with 0 for all elements
void initializeDisplacementArray(matrix_i &Ic); // NOLINT

void initializeStylusArray(matrix_i &st, int t); // NOLINT

// most outer loop, for we need a n*n loop for the algorithm
void calculation_loop(matrix_i &Ic, const matrix_i &Pa, double cell_size, // NOLINT
                 double v, double E);

// inner loop calling the calulation for every n
inline double calculate(double a, double b, double x, double y);

inline double calculate(int i, int j, double dxc, double dyc,
    double dxf, double dyf);

// actual calculation
double calc_displacement(const matrix_i &pressure,
                const matrix_i &Ic,
                double y, double x,
                double a, double b,
                double v, double E,
                double cell);

void calc_displacement(const matrix_i &pF, double cS, double fS,
                         matrix_i &cD);

void calcCoarsePressure(const std::vector<int>& qs, std::vector<matrix_i> &pFVec, // NOLINT
                        std::vector<matrix_i> &cDVec, // NOLINT
                        int t, const matrix_i& st);
// printing function for the data container
void printarray(const matrix_i &array);

// coarse pressure array calculation
void calc_coarse_pressure(const matrix_i &fP, const matrix_i &st,
                          matrix_i &cP, std::size_t t); // NOLINT

// initialization with 1 for all elements
void initializeMultiplicationArray(matrix_i &Ic); // NOLINT

// calculate the sizes for the coarse grids
void prepareCoarseSizes(std::vector<std::size_t> &gridLen1, // NOLINT
                std::vector<std::size_t> &gridLen2, // NOLINT
                std::vector<std::size_t> &coarseGrid, // NOLINT
                const std::size_t shape,
                const std::size_t t);


bool boundaryCheck(matrix_i &m, int i, int j); // NOLINT

void correctionSteps(matrix_i& cC, const matrix_i& st, int mc, // NOLINT
                    int t, double fineSizeA,
                    double fineSizeB, double halfSize);


void applyCorrection(matrix_i &coarseDisplacement, const matrix_i cC, // NOLINT
                     const matrix_i Ip, int t);

double correctionHelper(const matrix_i& cC, const matrix_i& Ip, int t,
                        int i, int j);

void interpolateGrid(matrix_i &nextDisplacement, // NOLINT
                     const matrix_i coarseDisplacement, // NOLINT
                     const matrix_i st);

void secondCorrectionStep(int mc, const matrix_i& st, double fineSizeA,
                          double fineSizeB, double hS, const matrix_i& // NOLINT
                          pF, matrix_i& cD, const std::vector<matrix_i> &cCVec); // NOLINT


void createCorrectionArrays(std::vector<matrix_i> &cCVec, // NOLINT
                            const matrix_i &st, double hS,
                            double fineSizeA, double fineSizeB,
                            int mc);
#endif  // MLMS_H_
