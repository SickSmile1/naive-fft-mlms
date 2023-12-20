/* Copyright 2023 <Ilia Fedotov @ Uni Freiburg> */
#include "Boussinesq.h"
#include "BoussinesqFft.h"
#include <fftw3.h>
#include <cmath>

/*static fftw_plan r2c = NULL;
static fftw_plan c2r = NULL;*/

// __________________________________________________________________
void copyPressureArray(matrix& p, const matrix& tempP) {// NOLINT 
  for (int i = 0; i < tempP.rows(); i++) {
    for (int j = 0; j < tempP.cols(); j++) {
      p(i, j) = tempP(i, j);
    }
  }
}

// __________________________________________________________________
void calculateGmn(matrix &Gmn, double dx, double dy) { // NOLINT
  int n = Gmn.rows();
  int oShape = n/2+1;
  int shape = (n/2)-1;
  Gmn.setZero();
  for (int i = 0; i < oShape; i++) {
    for (int j = 0; j < oShape; j++) {
      double res = calcBoussinesq(i, j, dx, dy, dx, dy);
      Gmn(i, j) = res;
    }
  }
  /* what actually is in the paper
  Gmn.block(0,oShape,oShape,shape+1) = -1*Gmn.block(0,1,oShape,shape+1).rowwise().reverse();
  Gmn.block(oShape,0,shape+1,oShape) = -1*Gmn.block(1,0,shape+1,oShape).colwise().reverse();
  Gmn.block(oShape,oShape,shape+1,shape+1) = -1*Gmn.block(1,oShape,shape+1,shape+1).colwise().reverse();
  what works: */
  Gmn.block(0,oShape,oShape,shape+1) = Gmn.block(0,1,oShape,shape+1).rowwise().reverse();
  Gmn.block(oShape,0,shape+1,oShape+shape+1) = Gmn.block(1,0,shape+1,oShape+shape+1).colwise().reverse();
  // Gmn.array() = Gmn.array()/Gmn.array();
  // writeToFile(Gmn, "gmn"+std::to_string(dx));
}

/*void createPlan_r2c(matrix& Gmn, cMatrix& Gmn_tild, fftw_plan& r2c,
                int mode) {
  auto plan = FFTW_ESTIMATE; // use fftw_wisdom with save/load?
  if (mode==1) plan = FFTW_PATIENT;
  r2c = fftw_plan_dft_r2c_2d(Gmn.rows(), Gmn.cols(), Gmn.data(),
      reinterpret_cast<fftw_complex*>(Gmn_tild.data()), plan);
}

void createPlan_c2r(matrix& Umn, cMatrix& Umn_tild, fftw_plan& c2r,
                int mode) {
  auto plan = FFTW_ESTIMATE;
  if (mode==1) plan = FFTW_PATIENT;
  c2r = fftw_plan_dft_c2r_2d(Umn.rows(), Umn.cols(),
                            reinterpret_cast<fftw_complex*>
                            (Umn_tild.data()),
                            Umn.data(), plan);
}*/
// __________________________________________________________________
void transformGmnP(matrix& Gmn, cMatrix& Gmn_tild, // NOLINT
                  matrix& p, cMatrix& p_tild) { // NOLINT
  fftw_plan p1;
  p1 = fftw_plan_dft_r2c_2d(Gmn.rows(), Gmn.cols(), Gmn.data(),
      reinterpret_cast<fftw_complex*>(Gmn_tild.data()), FFTW_ESTIMATE);
  // output array needs to be 2*nx / (ny*2/2)-1

  fftw_plan p2;
  p2 = fftw_plan_dft_r2c_2d(p.rows(), p.cols(), p.data(),
      reinterpret_cast<fftw_complex*>(p_tild.data()), FFTW_ESTIMATE);

  fftw_execute(p1);
  fftw_execute(p2);


  fftw_destroy_plan(p1);
  fftw_destroy_plan(p2);

  // https://www.fftw.org/doc/New_002darray-Execute-Functions.html#New_002darray-Execute-Functions
  // example: https://cplusplus.com/forum/general/285603/ */
}

// __________________________________________________________________
void transformToReal(cMatrix& Umn_tild, matrix& Umn) { // NOLINT
  fftw_plan p3; // TODO: use in place array
  p3 = fftw_plan_dft_c2r_2d(Umn.rows(), Umn.cols(),
                            reinterpret_cast<fftw_complex*>
                            (Umn_tild.data()),
                            Umn.data(), FFTW_ESTIMATE);
  fftw_execute(p3);
  // Umn = Umn_tild.real();
  fftw_destroy_plan(p3);
  /*fftw_execute_dft_c2r(c2r, reinterpret_cast<fftw_complex*>(Umn_tild.data()),
                       Umn.data());*/
}

// __________________________________________________________________
void writeToResultArray(matrix& Umn, matrix& Umn_res) { //NOLINT
  Umn_res = Umn.block(0,0, Umn_res.rows(), Umn_res.cols());
  Umn_res.array() = Umn_res.array()/(Umn.rows()* Umn.cols());
}

// __________________________________________________________________
void multiplyTransformed(cMatrix& Gmn_tild, cMatrix& Umn_tild, // NOLINT
                        cMatrix& p_tild) { // NOLINT
  Umn_tild.array() = Gmn_tild.array() * p_tild.array();
}


// __________________________________________________________________
matrix BoussinesqFFT(const double size, const int grid) {
  double Lx = size, Ly = size;
  int Nx = grid, Ny = grid;
  double pSize = size/2.;
  double dx = (Lx/Nx);
  double dy = (Ly/Ny);

  int lb = Nx/2-(pSize/dx)/2;
  int ub = Ny/2+(pSize/dy)/2;
  matrix tempP({Nx, Ny});
  initializePressureArray(tempP, lb, ub, 1.);
  matrix Gmn({(2*Nx)-1, (2*Ny)-1});
  return BoussinesqFFT(size, Gmn, tempP);
}

matrix BoussinesqFFT(const double size, matrix& surf, const matrix& topo) {
  matrix Gmn = Eigen::MatrixXd(surf);
  matrix tempP = Eigen::MatrixXd(topo);

  double Lx = size, Ly = size;
  int Nx = tempP.rows(), Ny = tempP.cols();
  double dx = (Lx/Nx);
  double dy = (Ly/Ny);
  // std::cout << "physical size: " << (Lx*Ly)/(Ny*Nx) << std::endl;
  cMatrix Gmn_tild({Gmn.rows(), Gmn.cols()/2+1});
  matrix p({Gmn.rows(), Gmn.cols()});
  cMatrix p_tild({Gmn.rows(), Gmn.cols()/2+1});
  matrix Umn({Gmn.rows(), Gmn.cols()});
  cMatrix Umn_tild({Gmn.rows(), Gmn.cols()/2+1});
  matrix Umn_res({tempP.rows(), tempP.cols()});

  copyPressureArray(p, tempP);

  calculateGmn(Gmn, dx, dy);

  transformGmnP(Gmn, Gmn_tild, p, p_tild);

  multiplyTransformed(Gmn_tild, Umn_tild, p_tild);

  transformToReal(Umn_tild, Umn);

  writeToResultArray(Umn, Umn_res);

  return Umn_res;
}

/*void BoussinesqFFT(const double size, matrix& surf, const matrix& topo,
                    fftw_plan& r2c, fftw_plan& c2r) {
  auto Gmn = surf;
  auto tempP = topo;

  double Lx = size, Ly = size;
  int Nx = tempP.rows(), Ny = tempP.cols();
  double dx = (Lx/Nx);
  double dy = (Ly/Ny);
  // std::cout << "physical size: " << (Lx*Ly)/(Ny*Nx) << std::endl;
  cMatrix Gmn_tild({Gmn.rows(), Gmn.cols()/2+1});
  matrix p({Gmn.rows(), Gmn.cols()});
  cMatrix p_tild({Gmn.rows(), Gmn.cols()/2+1});
  matrix Umn({Gmn.rows(), Gmn.cols()});
  cMatrix Umn_tild({Gmn.rows(), Gmn.cols()/2+1});
  matrix Umn_res({tempP.rows(), tempP.cols()});

  copyPressureArray(p, tempP);

  calculateGmn(Gmn, dx, dy);

  if(r2c == nullptr) {
    createPlan_r2c(Gmn, Gmn_tild, r2c, 0);
    fftw_execute_dft_r2c(r2c, Gmn.data(), reinterpret_cast<fftw_complex*>(Gmn_tild.data()));
    fftw_execute_dft_r2c(r2c, p.data(), reinterpret_cast<fftw_complex*>(p_tild.data()));
  }

  multiplyTransformed(Gmn_tild, Umn_tild, p_tild);

  if(c2r == nullptr) {
    createPlan_c2r(Umn, Umn_tild, c2r, 0);
  }
  transformToReal(Umn_tild, Umn);

  writeToResultArray(Umn, Umn_res, std::pow(size,2));

  Umn_res;
}*/