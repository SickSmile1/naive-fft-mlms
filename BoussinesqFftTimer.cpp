/* Copyright 2023 <Ilia Fedotov @ Uni Freiburg> */
#include "BoussinesqFftTimer.h"
#include "Boussinesq.h"
#include <fftw3.h>
#include <cmath>
#include <complex>
#include <iostream>

void BoussinesqFFT() {
  double Lx = 2., Ly = 2.;
  int Nx = 16, Ny = 16;
  double pSize = 1;
  double dx = (Lx/Nx);
  double dy = (Ly/Ny);

  int lb = Nx/2-(pSize/dx)/2;
  int ub = Ny/2+(pSize/dy)/2;

  matrix tempP({Nx, Ny});
  matrix p({2*Nx, 2*Ny});
  cMatrix p_tild({2*Nx, 2*Ny});
  matrix Gmn({2*Nx, 2*Ny});
  cMatrix Gmn_tild({2*Nx, 2*Ny});
  matrix Umn({Nx*2, Ny*2});
  cMatrix Umn_tild({2*Nx, 2*Ny});
  matrix Umn_res({Nx, Ny});

  initializePressureArray(tempP, lb, ub, 1.);
  initializeDisplacementArray(p);

  copyPressureArray(p, tempP);

  calculateGmn(Gmn, dx/2, dy/2);

  transformGmnP(Nx, Ny, Gmn, Gmn_tild, p, p_tild);

  multiplyTransformed(Gmn_tild, Umn_tild, p_tild);

  transformToReal(Umn_tild, Umn, Nx, Ny);

  writeToResultArray(Umn, Umn_res, Nx, Ny);

  writeToFile(Umn_res, "fastestGrid");
}
