/* Copyright 2023 <Ilia Fedotov @ Uni Freiburg> */
#ifndef BOUSSINESQ_H_
#define BOUSSINESQ_H_

// #include <gtest/gtest.h>
#include <vector>
#include <array>
#include <cassert>
#include <complex>
#include <sstream>
#include <iostream>
#include <fstream>
#include <string>
// #include <eigen3/Eigen/Dense>
// #include<eigen3/Eigen/StdVector>

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
    // assert(i < shape[0] && j < shape[1]);
    return data[j+(i*shape[1])];
  }

  // overload () for aquiring (conversion) i,j indexes in 1d array
  // in static form
  const T& operator() (int i, int j) const {
    // assert(i < shape[0] && j < shape[1]);
    return data[j+(i*shape[1])];
  }
};

using matrix = matrixTemplate<double>;
using smatrix = matrixTemplate<std::string>;
using cMatrix = matrixTemplate<std::complex<double>>;

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

// initialization of pressure array
void initializePressureArray(matrix &Pa, double lower_b, // NOLINT
                             double upper_b, double pressure); // NOLINT

// initialization with 0 for all elements
void initializeDisplacementArray(matrix &Ic); // NOLINT

double calcBoussinesq(double a, double b, double x, double y);

double calcBoussinesq(int i, int j, double dxc, double dyc,
    double dxf, double dyf);

void naiveCalculation(matrix &Ic, const matrix &Pa, // NOLINT
                double cell_size);

double calc_displacement(const matrix &pressure,
              const matrix &Ic,
              double y, double x,
              double a, double b, double cell);

void calc_displacement(const matrix &pF, double cS, double fS,
                         matrix &cD);
#endif  // BOUSSINESQ_H_
