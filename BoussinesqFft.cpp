/* Copyright 2023 <Ilia Fedotov @ Uni Freiburg> */
#include "Boussinesq.h"
#include "BoussinesqFft.h"
#include <fftw3.h>

// __________________________________________________________________
void copyPressureArray(matrix& p, const matrix& tempP) {// NOLINT 
  for (int i = 0; i < tempP.shape[0]; i++) {
    for (int j = 0; j < tempP.shape[1]; j++) {
      p(i, j) = tempP(i, j);
    }
  }
}

// __________________________________________________________________
void calculateGmn(matrix &Gmn, double dx, double dy) { // NOLINT
  int shape = (Gmn.shape[0])/2;
  for (int i = 0; i <= shape; i++) {
    double res = calcBoussinesq(i, 0, dx, dy, dx, dy);
    Gmn(i, 0) = res;
    Gmn(0, i) = res;
    if (i>0) {Gmn((2*shape+1-i), 0) = res;
    Gmn(0, (2*shape+1)-i) =  res;}
  }
  for (int i = 1; i <= shape; i++) {
    for (int j = 1; j <= shape; j++) {
      double res = calcBoussinesq(i, j, dx, dy, dx, dy);
      Gmn(i, (2*shape+1)-j) =  res;
      Gmn(i, j) =  res;
      Gmn((2*shape+1)-i, j) = res;
      Gmn((2*shape+1)-i, (2*shape+1)-j) = res;
    }
  }
}

// __________________________________________________________________
void transformGmnP(int Nx, int Ny, matrix& Gmn, cMatrix& Gmn_tild, // NOLINT
                  matrix& p, cMatrix& p_tild) { // NOLINT
  fftw_plan p1;
  p1 = fftw_plan_dft_r2c_2d(Gmn.shape[0], Gmn.shape[1], Gmn.data.data(),
      reinterpret_cast<fftw_complex*>(Gmn_tild.data.data()), FFTW_ESTIMATE);
  // output array needs to be 2*nx / (ny*2/2)-1

  fftw_plan p2;
  p2 = fftw_plan_dft_r2c_2d(p.shape[0], p.shape[1], p.data.data(),
      reinterpret_cast<fftw_complex*>(p_tild.data.data()), FFTW_ESTIMATE);

  fftw_execute(p1);
  fftw_execute(p2);

  fftw_destroy_plan(p1);
  fftw_destroy_plan(p2);
}

// __________________________________________________________________
void transformToReal(cMatrix& Umn_tild, matrix& Umn, int Nx, int Ny) { // NOLINT
  fftw_plan p3;
  p3 = fftw_plan_dft_c2r_2d(Umn.shape[0], Umn.shape[1],
                            reinterpret_cast<fftw_complex*>
                            (Umn_tild.data.data()),
                            Umn.data.data(), FFTW_ESTIMATE);
  fftw_execute(p3);
  fftw_destroy_plan(p3);
}

// __________________________________________________________________
void writeToResultArray(const matrix& Umn, matrix& Umn_res, // NOLINT
                        int Nx, int Ny) { //NOLINT
  for (int i = 0; i < Umn_res.shape[0]; i++) {
    for (int j = 0; j < Umn_res.shape[1]; j++) {
      // devide each result by Nx*Ny 
      Umn_res(i, j) = Umn(i, j)/(Umn.shape[0]*Umn.shape[1]);
    }
  }
  // writeToFile(Umn_res, "whatthe");
}

// __________________________________________________________________
void multiplyTransformed(cMatrix& Gmn_tild, cMatrix& Umn_tild, // NOLINT
                        cMatrix& p_tild) { // NOLINT
  for (int i = 0; i < Gmn_tild.shape[0]; i++) {
    for (int j = 0; j < Gmn_tild.shape[1]; j++) {
      Umn_tild(i, j) = Gmn_tild(i, j)*p_tild(i, j);
    }
  }
}


// __________________________________________________________________
matrix BoussinesqFFT(double size, int grid) {
  double Lx = size, Ly = size;
  int Nx = grid, Ny = grid;
  double pSize = size/2.;
  double dx = (Lx/Nx);
  double dy = (Ly/Ny);

  int lb = Nx/2-(pSize/dx)/2;
  int ub = Ny/2+(pSize/dy)/2;

  matrix Gmn({(2*Nx)-1, (2*Ny)-1});
  cMatrix Gmn_tild({Gmn.shape[0], Gmn.shape[1]/2+1});
  matrix p({Gmn.shape[0], Gmn.shape[1]});
  cMatrix p_tild({Gmn.shape[0], Gmn.shape[1]/2+1});
  matrix tempP({Nx, Ny});
  matrix Umn({Gmn.shape[0], Gmn.shape[1]});
  cMatrix Umn_tild({Gmn.shape[0], Gmn.shape[1]/2+1});
  matrix Umn_res({Nx, Ny});

  initializePressureArray(tempP, lb, ub, 1.);
  initializeDisplacementArray(p);

  copyPressureArray(p, tempP);

  calculateGmn(Gmn, dx, dy);

  transformGmnP(Nx, Ny, Gmn, Gmn_tild, p, p_tild);

  multiplyTransformed(Gmn_tild, Umn_tild, p_tild);

  transformToReal(Umn_tild, Umn, Nx, Ny);

  // writeToResultArray(Umn, Umn_res, Nx, Ny);
  
  return Umn_res;
}
