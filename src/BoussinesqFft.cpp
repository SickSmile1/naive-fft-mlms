/* Copyright 2023 <Ilia Fedotov @ Uni Freiburg> */
#include "Boussinesq.h"
#include "BoussinesqFft.h"
#include <fftw3.h>

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
  int n = Gmn.rows()+1;
  int oShape = (n-1)/2 + 1;
  int shape = (n/2)-1;
  std::cout << "gmn shapes "<< oShape << " : " << shape << "\n" << std::endl;
  Gmn.setZero();
  for (int i = 0; i <= oShape; i++) {
    for (int j = 0; j <= oShape; j++) {
      double res = calcBoussinesq(i, j, dx, dy, dx, dy);
      Gmn(j, i) = res;
      Gmn(i, j) = res;
    }
  }
  Gmn.block(0,oShape,oShape,oShape) = Gmn.block(0,0,oShape,shape).rowwise().reverse();
  Gmn.block(oShape,0,shape,shape+oShape) = Gmn.block(0,0,shape,shape+oShape).reverse();
  writeToFile(Gmn, "gmn");
}

// BACKUP
/*void calculateGmn(matrix &Gmn, double dx, double dy) { // NOLINT
  int shape = (Gmn.rows())/2;
  Gmn.setZero();
  for (int i = 0; i <= shape; i++) {
    double res = calcBoussinesq(i, 0, dx, dy, dx, dy);
    Gmn(i, 0) = res;
    Gmn(0, i) = res;
    if (i > 0) {
      Gmn((2*shape+1-i), 0) = res;
      Gmn(0, (2*shape+1)-i) = res;
    }
  }
  for (int i = 1; i <= shape; i++) {
    for (int j = 1; j <= shape; j++) {
      double res = calcBoussinesq(i, j, dx, dy, dx, dy);
      Gmn(i, (2*shape+1)-j) = res;
      Gmn(i, j) = res;
      Gmn((2*shape+1)-i, j) = res;
      Gmn((2*shape+1)-i, (2*shape+1)-j) = res;
    }
  }
}*/

// __________________________________________________________________
void transformGmnP(matrix& Gmn, cMatrix& Gmn_tild, // NOLINT
                  matrix& p, cMatrix& p_tild) { // NOLINT
  fftw_plan p1; // TODO: reuse plan? use inplace array
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
}

// __________________________________________________________________
void transformToReal(cMatrix& Umn_tild, matrix& Umn) { // NOLINT
  fftw_plan p3; // TODO: use in place array
  p3 = fftw_plan_dft_c2r_2d(Umn.rows(), Umn.cols(),
                            reinterpret_cast<fftw_complex*>
                            (Umn_tild.data()),
                            Umn.data(), FFTW_ESTIMATE);
  fftw_execute(p3);
  fftw_destroy_plan(p3);
}

// __________________________________________________________________
void writeToResultArray(const matrix& Umn, matrix& Umn_res) { //NOLINT
  for (int i = 0; i < Umn_res.rows(); i++) {
    for (int j = 0; j < Umn_res.cols(); j++) {
      // devide each result by Nx*Ny
      Umn_res(i, j) = Umn(i, j)/(Umn.rows()*Umn.cols()); 
      // #TODO replace by simplified eigen operatiom
    }
  }
  // writeToFile(Umn_res, "whatthe");
}

// __________________________________________________________________
void multiplyTransformed(cMatrix& Gmn_tild, cMatrix& Umn_tild, // NOLINT
                        cMatrix& p_tild) { // NOLINT
  /*for (int i = 0; i < Gmn_tild.rows(); i++) {
    for (int j = 0; j < Gmn_tild.cols(); j++) {
      Umn_tild(i, j) = Gmn_tild(i, j)*p_tild(i, j);
    }
  }*/
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
  // matrix Gmn({(2*Nx)-1, (2*Ny)-1});
  matrix Gmn({ ((Nx-1)/2)+((Nx/2)-1) , ((Ny-1)/2)+((Ny/2)-1)});
  std::cout << "nx shape " << Nx << " : " << Ny << "\n" << std::endl;
  std::cout << "constr shape " << Gmn.rows() << " : " << Gmn.cols() << "\n" << std::endl;
  return BoussinesqFFT(size, Gmn, tempP);
}

matrix BoussinesqFFT(const double size, matrix& surf, const matrix& topo) {
  
  auto Gmn = surf;
  auto tempP = topo;

  double Lx = size, Ly = size;
  int Nx = tempP.rows(), Ny = tempP.cols();
  double dx = (Lx/Nx);
  double dy = (Ly/Ny);

  cMatrix Gmn_tild({Gmn.rows(), Gmn.cols()/2+1});
  matrix p({Gmn.rows(), Gmn.cols()});
  cMatrix p_tild({Gmn.rows(), Gmn.cols()/2+1});
  matrix Umn({Gmn.rows(), Gmn.cols()});
  cMatrix Umn_tild({Gmn.rows(), Gmn.cols()/2+1});
  matrix Umn_res({tempP.rows(), tempP.cols()});


  initializeDisplacementArray(p);

  copyPressureArray(p, tempP);

  calculateGmn(Gmn, dx, dy);

  transformGmnP(Gmn, Gmn_tild, p, p_tild);

  multiplyTransformed(Gmn_tild, Umn_tild, p_tild);

  transformToReal(Umn_tild, Umn);

  writeToResultArray(Umn, Umn_res);

  return Umn_res;
}
