/* Copyright 2023 <Ilia Fedotov @ Uni Freiburg> */
#ifndef BOUSSINESQ_H_
#define BOUSSINESQ_H_

// #include <gtest/gtest.h>
#include <vector>
#include <array>
#include <cassert>
#include <complex>
#include <fstream>
#include <string>

template <class T>
struct matrixTemplate{
  // data container for pressure arrays, grid and result
  std::array<int, 2> shape;
  std::vector<T> data;

  // constructor, to initialize simulated 2d array in a 1d shape
  explicit matrixTemplate(std::array<int, 2> shape_) : shape(shape_) {
    data.resize(shape[0]*shape[1]);
  }
  auto begin() {
    return data.begin();
  }

  auto end() {
    return data.end();
  }
  // overload () for aquiring (conversion) i,j idexes in 1d array
  T& operator() (int i, int j) {
    assert(i < shape[0] && j < shape[1]);
    return data[j+(i*shape[1])];
  }

  // overload () for aquiring (conversion) i,j indexes in 1d array
  // in static form
  const T& operator() (int i, int j) const {
    assert(i < shape[0] && j < shape[1]);
    return data[j+(i*shape[1])];
  }
};

using matrix = matrixTemplate<double>;
using cMatrix = matrixTemplate<std::complex<double>>;

// initialization of pressure array
void initializePressureArray(matrix &Pa, double lower_b, // NOLINT
                             double upper_b, double pressure); // NOLINT

// initialization with 0 for all elements
void initializeDisplacementArray(matrix &Ic); // NOLINT

void initializeStylusArray(matrix &st, int t); // NOLINT

// most outer loop, for we need a n*n loop for the algorithm
void naiveCalculation(matrix &Ic, const matrix &Pa, double cell_size); // NOLINT

// inner loop calling the calulation for every n
inline double calcBoussinesq(double a, double b, double x, double y);

inline double calcBoussinesq(int i, int j, double dxc, double dyc,
    double dxf, double dyf);

// actual calculation
double calc_displacement(const matrix &pressure,
                const matrix &Ic,
                double y, double x,
                double a, double b,
                double cell);

void calc_displacement(const matrix &pF, double cS, double fS,
                         matrix &cD);  // NOLINT

void calcCoarsePressure(const std::vector<int>& qs, std::vector<matrix> &pFVec, // NOLINT
                        std::vector<matrix> &cDVec, // NOLINT
                        int t, const matrix& st);
// printing function for the data container
// void printarray(const matrixTemplate<T> &array);

template <class T>
// __________________________________________________________________
void printarray(const matrixTemplate<T> &array) {
  for (int i = 0; i < array.shape[0]; i++) {
    for (int j = 0; j < array.shape[1]; j++) {
      printf("%.4f\t", array(i, j));
    }
    printf("\n");
  }
  return;
}

template <class T>
// __________________________________________________________________
void writeToFile(const matrixTemplate<T> &arr, const std::string name) {
  std::ofstream file;
  file.open(name, std::ofstream::out | std::ofstream::trunc);
  for (int i = 0; i < arr.shape[0]; i++) {
    for (int j = 0; j < arr.shape[1]; j++) {
      file << arr(i, j);
      file << "\t";
    }
    file << std::endl;
  }
  file.close();
}

// coarse pressure array calculation
void calc_coarse_pressure(const matrix &fP, const matrix &st,
                          matrix &cP, std::size_t t); // NOLINT

// initialization with 1 for all elements
void initializeMultiplicationArray(matrix &Ic); // NOLINT

// calculate the sizes for the coarse grids
void prepareCoarseSizes(std::vector<std::size_t> &gridLen1, // NOLINT
                std::vector<std::size_t> &gridLen2, // NOLINT
                std::vector<std::size_t> &coarseGrid, // NOLINT
                const std::size_t shape,
                const std::size_t t);


bool boundaryCheck(matrix &m, int i, int j); // NOLINT

void correctionSteps(matrix& cC, const matrix& st, int mc, // NOLINT
                    int t, double fineSizeA,
                    double fineSizeB, double halfSize);


void applyCorrection(matrix &coarseDisplacement, const matrix cC, // NOLINT
                     const matrix Ip, int t);

double correctionHelper(const matrix& cC, const matrix& Ip, int t,
                        int i, int j);

void interpolateGrid(matrix &nextDisplacement, // NOLINT
                     const matrix coarseDisplacement, // NOLINT
                     const matrix st);

void secondCorrectionStep(int mc, const matrix& st, double fineSizeA,
                          double fineSizeB, double hS, const matrix&  // NOLINT
                          pF, matrix& cD, const std::vector<matrix> &cCVec);  // NOLINT


void createCorrectionArrays(std::vector<matrix> &cCVec, // NOLINT
                            const matrix &st, double hS,
                            double fineSizeA, double fineSizeB,
                            int mc);

void copyPressureArray(matrix& p, const matrix& tempP);  // NOLINT

void calculateGmn(matrix &Gmn, double dx, double dy);  // NOLINT

void transformGmnP(int Nx, int Ny, matrix& Gmn, cMatrix& Gmn_tild,  // NOLINT
                  matrix& p, cMatrix& p_tild);  // NOLINT


void multiplyTransformed(cMatrix& Gmn_tild, cMatrix& Umn_tild,
                        cMatrix& p_tild); // NOLINT

void transformToReal(cMatrix& Umn_tild, matrix& Umn, int Nx, int Ny); // NOLINT

void writeToResultArray(const matrix& Umn, matrix& Umn_res,  // NOLINT
                        int Nx, int Ny);

#endif  // BOUSSINESQ_H_
