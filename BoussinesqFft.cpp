/* Copyright 2023 <Ilia Fedotov @ Uni Freiburg> */
#include "Boussinesq.h"
#include "BoussinesqFft.h"
#include <fftw3.h>

// __________________________________________________________________
void copyPressureArray(matrix& p, const matrix& tempP) {
  for (int i = 0; i < tempP.shape[0]; i++) {
    for (int j = 0; j < tempP.shape[1]; j++) {
      p(i, j) = tempP(i, j);
    }
  }  
}

// __________________________________________________________________
void calculateGmn(matrix &Gmn, double dx, double dy) { // NOLINT
  for (int i = 0; i <= (Gmn.shape[0])/2; i++) {
    for (int j = 0; j <= (Gmn.shape[0])/2; j++) {
      Gmn(i, j) = calcBoussinesq(i, j, dx, dy, dx, dy);
      Gmn(i, (Gmn.shape[0]-j)) = calcBoussinesq(i, -j, dx, dy, dx, dy);
      Gmn((Gmn.shape[0]-i), j) = calcBoussinesq(-i, j, dx, dy, dx, dy);
      Gmn((Gmn.shape[0]-i), (Gmn.shape[0]-j)) =
        calcBoussinesq(-i, -j, dx, dy, dx, dy);
      // std::cout << Gmn(i, j) << std::endl;
    }
  }
  printarray(Gmn);
}

// __________________________________________________________________
void transformGmnP(int Nx, int Ny, matrix& Gmn, cMatrix& Gmn_tild,
                  matrix& p, cMatrix& p_tild) {
  fftw_plan p1;
  p1 = fftw_plan_dft_r2c_2d(Nx, Ny*2, Gmn.data.data(),
      reinterpret_cast<fftw_complex*>(Gmn_tild.data.data()), FFTW_ESTIMATE);
  // output array needs to be 2*nx / (ny*2/2)-1

  fftw_plan p2;
  p2 = fftw_plan_dft_r2c_2d(Nx, Ny*2, p.data.data(),
      reinterpret_cast<fftw_complex*>(p_tild.data.data()), FFTW_ESTIMATE);

  fftw_execute(p1);
  fftw_execute(p2);

  fftw_destroy_plan(p1);
  fftw_destroy_plan(p2);
}

// __________________________________________________________________
void transformToReal(cMatrix& Umn_tild, matrix& Umn, int Nx, int Ny) {
  fftw_plan p3;
  p3 = fftw_plan_dft_c2r_2d(Nx, Ny*2,
                            reinterpret_cast<fftw_complex*>
                            (Umn_tild.data.data()),
                            Umn.data.data(), FFTW_ESTIMATE);
  fftw_execute(p3);
  fftw_destroy_plan(p3);
}

// __________________________________________________________________
void writeToResultArray(const matrix& Umn, matrix& Umn_res, 
                        int Nx, int Ny) {
  int N = Nx*Ny;
  #pragma omp parallel for simd
  for (int i = 1; i < Umn_res.shape[0]; i++) {
    for (int j = 0; j < Umn_res.shape[1]; j++) {
    Umn_res(i, j) = Umn(i, j)/N;
    }
  }
  writeToFile(Umn_res, "whatthe");
}

// __________________________________________________________________
void multiplyTransformed(cMatrix& Gmn_tild, cMatrix& Umn_tild,
                        cMatrix& p_tild) {
  for (int i = 0; i < Gmn_tild.shape[0]; i++) {
    for (int j = 0; j < Gmn_tild.shape[1]; j++) {
      Umn_tild(i, j) = Gmn_tild(i, j)*p_tild(i, j);
    }
  }
}
