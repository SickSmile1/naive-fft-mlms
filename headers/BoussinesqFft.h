/* Copyright 2023 <Ilia Fedotov @ Uni Freiburg> */
#ifndef BOUSSINESQFFT_H_
#define BOUSSINESQFFT_H_

#include "Boussinesq.h"

void copyPressureArray(matrix& p, const matrix& tempP); // NOLINT

void copyPressureArray(matrix& p, const matrix& tempP); // NOLINT

void calculateGmn(matrix &Gmn, double dx, double dy); // NOLINT

void transformGmnP(matrix& Gmn, cMatrix& Gmn_tild, // NOLINT
                  matrix& p, cMatrix& p_tild); // NOLINT

void multiplyTransformed(cMatrix& Gmn_tild, cMatrix& Umn_tild, // NOLINT
                        cMatrix& p_tild); // NOLINT

void transformToReal(cMatrix& Umn_tild, matrix& Umn); // NOLINT

void writeToResultArray(const matrix& Umn, matrix& Umn_res); // NOLINT

/* initialize with centered pressure patch*/
matrix BoussinesqFFT(double size, int grid);
/* calculate with ready surface and topography*/
matrix BoussinesqFFT(double size, matrix surf, matrix topo);

#endif  // BOUSSINESQFFT_H_
