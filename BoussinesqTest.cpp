/* Copyright 2023 <Ilia Fedotov @ Uni Freiburg> */

#include <gtest/gtest.h>
#include <iostream>
#include <fstream>
#include <thread> //NOLINT
#include <string>
#include "./Boussinesq.h"
#include "./BoussinesqTimer.h"

/*TEST(filename, function){
  function values
  ASSERT_EQ(res, function call);
}*/

matrix res({20, 20});

TEST(BoussinesqNaive, calculate) {
  // test for symmetry of resulting matrix
  double size = 2;
  double size_p = 1;
  int grid = 20;
  double cell_size = size/grid;

  matrix Ic({grid, grid});
  matrix Pa({grid, grid});

  double pressure = 1.;
  double lower_b = (size/cell_size)/2 - (size_p/cell_size)/2;
  double upper_b = (size/cell_size)/2 + (size_p/cell_size)/2;

  initializePressureArray(Pa, lower_b, upper_b, pressure);
  initializeDisplacementArray(Ic);
  naiveCalculation(Ic, Pa, cell_size);

  bool equal = false;
  double eps = 0.001;
  for (int i = 0; i < Ic.shape[0]; i++) {
    for (int j = 0; j < Ic.shape[1]; j++) {
      // std::cout << Ic(i, j) << " : " << Ic(j, i) << std::endl;
      equal = std::abs(Ic(i, j) - Ic(j, i)) < eps;
      res(i, j) = Ic(i, j);
      if (equal == false) {
        std::cout << Ic(i, j) << " : " << Ic(j, i) << std::endl;
        break;
      }
    }
  }

  EXPECT_TRUE(equal);
}

TEST(BoussinesqTimer, timeTest) {
  MeasureTime mt = MeasureTime();
  std::this_thread::sleep_for(std::chrono::milliseconds(1000));
  auto stopped =  mt.stopTime();
  EXPECT_GE(stopped, 1);
}

TEST(WriteToFile, writeToFile) {
  // write/read file and check if values match
  matrix test({3, 3});
  int counter = 0;
  for (int i = 0; i < test.shape[0]; i++) {
    for (int j = 0; j < test.shape[1]; j++) {
      test(i, j) = counter++;
    }
  }

  writeToFile(test, "test.txt");
  std::ifstream inp("test.txt");
  std::stringstream buff;
  buff << inp.rdbuf();
  double res;
  buff >> res;
  ASSERT_EQ(test(0, 0), res);
  buff >> res;
  ASSERT_EQ(test(0, 1), res);
  buff >> res;
  ASSERT_EQ(test(0, 2), res);
  buff >> res;
  ASSERT_EQ(test(1, 0), res);
  buff >> res;
  ASSERT_EQ(test(1, 1), res);
  buff >> res;
  ASSERT_EQ(test(1, 2), res);
  buff >> res;
  ASSERT_EQ(test(2, 0), res);
  buff >> res;
  ASSERT_EQ(test(2, 1), res);
  buff >> res;
  ASSERT_EQ(test(2, 2), res);
}

TEST(Boussinesq, boundaryCheck) {
  matrix m({2, 2});
  EXPECT_TRUE(boundaryCheck(m, 1, 1));
  EXPECT_TRUE(boundaryCheck(m, 1, 0));
  EXPECT_TRUE(boundaryCheck(m, 0, 1));
  EXPECT_TRUE(boundaryCheck(m, 0, 0));
  EXPECT_FALSE(boundaryCheck(m, -1, 0));
  EXPECT_FALSE(boundaryCheck(m, 1, -1));
  EXPECT_FALSE(boundaryCheck(m, 3, 0));
  EXPECT_FALSE(boundaryCheck(m, 0, 3));
}

TEST(BoussinesqMlms, calculate) {
  const double size = 2;
  const double size_p = 1;
  const double pressure = 1.;

  // initial grid size for initialization
  int grid1 = 20;
  int grid2 = 20;

  double fineSizeA = size / grid1;
  double fineSizeB = size / grid2;

  matrix kM({grid1, grid2});
  matrix Ip({grid1, grid2});

  double lower_b = (grid1)/2. - (size_p/fineSizeA)/2.;
  double upper_b = (grid2)/2. + (size_p/fineSizeB)/2.;

  initializePressureArray(Ip, lower_b, upper_b, pressure);
  // material moduli and v

  double beta = 0.84;
  double min_g = std::min(grid1, grid2);
  // int t = beta*log(min_g);
  int t = 5;
  int mc = 0.7* pow(min_g, 1./t)-1;
  if (mc < 2*t) {
    mc = 2*t;
  }
  // std::cout << mc << std::endl;
  matrix st({2*t, 1});

  initializeStylusArray(st, t);

  std::vector<matrix> pfVec;
  std::vector<matrix> cDVec;

  std::vector<int> qs;

  double qLevel = std::log2(grid1 * grid2)/2-1;
  int q = grid1 / 2 + 2*t - 1;
  qs.push_back(q);

  for (int length = 1; length < qLevel-1; length++) {
    q = q / 2 + 2*t - 1;
    qs.push_back(q);
  }

  double d = qs.size();
  pfVec.reserve(d);
  cDVec.reserve(d);
  pfVec.push_back(Ip);
  cDVec.push_back(kM);
  double coarseSize = fineSizeA*pow(2, d);
  std::vector<matrix> cCVec;
  cCVec.reserve(3);
  createCorrectionArrays(cCVec, st, coarseSize, fineSizeA,
                          fineSizeB, mc);

  calcCoarsePressure(qs, pfVec, cDVec, t, st);

  calc_displacement(pfVec[d], coarseSize, fineSizeA, cDVec[d]);

  for (int i = 0; i < qs.size(); i++) {
    double hS = fineSizeA*pow(2, d-i-1);
    int temp_mc = (mc*2)+1;
    matrix cC({temp_mc, temp_mc});
    correctionSteps(cC, st, mc, t, fineSizeA, fineSizeB, hS);
    applyCorrection(cDVec[d-i], cC, pfVec[d-i-1], t);
    interpolateGrid(cDVec[d-i-1], cDVec[d-i], st);
    secondCorrectionStep(mc, st, fineSizeA, fineSizeB, hS,
                          pfVec[d-i-1], cDVec[d-i-1], cCVec);
  }

  bool equal = false;
  bool equal1 = false;
  double eps = 0.005;
  for (int i = 0; i < cDVec[0].shape[0]; i++) {
    for (int j = 0; j < cDVec[0].shape[1]; j++) {
      equal = std::abs(cDVec[0](i, j) - cDVec[0](j, i)) < eps;
      equal1 = std::abs(res(i, j) - cDVec[0](j, i)) < eps;
      if (equal == false) {
        std::cout << cDVec[0](i, j) << " : " << cDVec[0](j, i) << std::endl;
        break;
      }
    }
  }
  EXPECT_TRUE(equal);
  EXPECT_TRUE(equal1);
}

TEST(BoussinesqFFT, calculate) {
  double Lx = 2., Ly = 2.;
  int Nx = 20, Ny = 20;
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

  transformToReam(Umn_tild, Umn, Nx, Ny);

  writeToResultArray(Umn, Umn_res, Nx, Ny);

  bool equal = false;
  bool equal1 = false;
  double eps = 0.001;
  for (int i = 0; i < Umn_res.shape[0]; i++) {
    for (int j = 0; j < Umn_res.shape[1]; j++) {
      equal = std::abs(Umn_res(i, j) - Umn_res(j, i)) < eps;
      equal1 = std::abs(res(i, j) - Umn_res(j, i)) < eps;
      if (equal1 == false) {
        std::cout << res(i, j) << 
                        " : " << Umn_res(i, j) <<
                        " : " << i << " : " << j <<
                        std::endl;
        break;
      }
    }
  }
  EXPECT_TRUE(equal);
  EXPECT_TRUE(equal1);
}