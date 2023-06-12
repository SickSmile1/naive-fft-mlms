/* Copyright 2023 <Ilia Fedotov @ Uni Freiburg> */
#ifndef BOUSSINESQFFT_H_
#define BOUSSINESQFFT_H_

#include "Boussinesq.h"

void copyPressureArray(matrix& p, const matrix& tempP); // NOLINT

void copyPressureArray(matrix& p, const matrix& tempP); // NOLINT

void calculateGmn(matrix &Gmn, double dx, double dy); // NOLINT

void transformGmnP(int Nx, int Ny, matrix& Gmn, cMatrix& Gmn_tild, // NOLINT
                  matrix& p, cMatrix& p_tild); // NOLINT

void multiplyTransformed(cMatrix& Gmn_tild, cMatrix& Umn_tild, // NOLINT
                        cMatrix& p_tild); // NOLINT

void transformToReal(cMatrix& Umn_tild, matrix& Umn, int Nx, int Ny); // NOLINT

void writeToResultArray(const matrix& Umn, matrix& Umn_res, // NOLINT
                        int Nx, int Ny); // NOLINT

matrix BoussinesqFFT(double size, int grid); 

#endif  // BOUSSINESQFFT_H_
