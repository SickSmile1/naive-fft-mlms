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
#include <iomanip>
#include <Eigen/Core>

using matrix = Eigen::MatrixXd;
using cMatrix = Eigen::MatrixXcd;

// printing function for the data container
template <class T>
// __________________________________________________________________
void printarray(const T &array) {
  std::cout << array << std::endl;
  return;
}

template <class T>
// __________________________________________________________________
void writeToFile(const T &arr, const std::string name) {
  std::ofstream file;
  file.open(name, std::ofstream::out | std::ofstream::trunc);
  file << std::setprecision (25);
  file << arr << std::endl;
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
                double fS, double cS);

double calc_displacement(const matrix &pressure,
              const matrix &Ic,
              double y, double x,
              double a, double b, double cell);

void calc_displacement(const matrix &pF, double cS, double fS,
                         matrix &cD); // NOLINT
#endif  // BOUSSINESQ_H_
