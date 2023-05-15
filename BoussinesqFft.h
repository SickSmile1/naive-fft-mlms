/* Copyright 2023 <Ilia Fedotov @ Uni Freiburg> */
#ifndef BOUSSINESQFFT_H_
#define BOUSSINESQFFT_H_

#include "Boussinesq.h"

void copyPressureArray(matrix& p, const matrix& tempP);

void copyPressureArray(matrix& p, const matrix& tempP);

void calculateGmn(matrix &Gmn, double dx, double dy);

void transformGmnP(int Nx, int Ny, matrix& Gmn, cMatrix& Gmn_tild,
                  matrix& p, cMatrix& p_tild);

void multiplyTransformed(cMatrix& Gmn_tild, cMatrix& Umn_tild,
                        cMatrix& p_tild);

void transformToReal(cMatrix& Umn_tild, matrix& Umn, int Nx, int Ny);

void writeToResultArray(const matrix& Umn, matrix& Umn_res, 
                        int Nx, int Ny);

#endif  // BOUSSINESQFFT_H_