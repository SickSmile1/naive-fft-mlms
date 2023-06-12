/* Copyright 2023 <Ilia Fedotov @ Uni Freiburg> */

#include <gtest/gtest.h>
#include <iostream>
#include <fstream>
#include <thread> //NOLINT
#include <string>
#include "Boussinesq.h"
#include "BoussinesqMlms.h"
#include "BoussinesqFft.h"
#include "BoussinesqTimer.h"

/*TEST(filename, function){
  function values
  ASSERT_EQ(res, function call);
}*/

int grids = 64;
matrix res1({grids, grids});
matrix resMlms_6({grids, grids});
matrix resFft({grids, grids});

TEST(BoussinesqNaive, calculate) {
  // GTEST_SKIP();
  // test for symmetry of resulting matrix
  double size = 2;
  double size_p = 1;
  int grid = grids;
  double cell_size = size/grid;

  matrix Ic({grid, grid});
  matrix Pa({grid, grid});

  double pressure = 1.;
  double lower_b = (size/cell_size)/2 - (size_p/cell_size)/2;
  double upper_b = (size/cell_size)/2 + (size_p/cell_size)/2;

  initializePressureArray(Pa, lower_b, upper_b, pressure);
  initializeDisplacementArray(Ic);
  // naiveCalculation(Ic, Pa, cell_size/2, cell_size);
  calc_displacement(Pa, cell_size, cell_size, Ic);

  bool equal = false;
  double eps = 0.001;
  for (int i = 0; i < Ic.shape[0]; i++) {
    for (int j = 0; j < Ic.shape[1]; j++) {
      // std::cout << Ic(i, j) << " : " << Ic(j, i) << std::endl;
      equal = std::abs(Ic(i, j) - Ic(j, i)) < eps;
      res1(i, j) = Ic(i, j);
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
  double res2;
  buff >> res2;
  ASSERT_EQ(test(0, 0), res2);
  buff >> res2;
  ASSERT_EQ(test(0, 1), res2);
  buff >> res2;
  ASSERT_EQ(test(0, 2), res2);
  buff >> res2;
  ASSERT_EQ(test(1, 0), res2);
  buff >> res2;
  ASSERT_EQ(test(1, 1), res2);
  buff >> res2;
  ASSERT_EQ(test(1, 2), res2);
  buff >> res2;
  ASSERT_EQ(test(2, 0), res2);
  buff >> res2;
  ASSERT_EQ(test(2, 1), res2);
  buff >> res2;
  ASSERT_EQ(test(2, 2), res2);
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

TEST(BoussinesqMlms, initializeStylusArray) {
  matrix st = initializeStylusArray(1);
  ASSERT_EQ(st(0, 0), (1./2.));
  ASSERT_EQ(st(1, 0), (1./2.));

  matrix st1 = initializeStylusArray(2);
  ASSERT_EQ(st1(0, 0), (-1./16.) );
  ASSERT_EQ(st1(1, 0), (1./16.)*9);
  ASSERT_EQ(st1(2, 0), (1./16.)*9);
  ASSERT_EQ(st1(3, 0), (-1/16.) );

  matrix st2 = initializeStylusArray(3);
  ASSERT_EQ(st2(0, 0), (1./256)*3);
  ASSERT_EQ(st2(1, 0), -(1./256)*25);
  ASSERT_EQ(st2(2, 0), (1./256)*150);
  ASSERT_EQ(st2(3, 0), (1./256)*150);
  ASSERT_EQ(st2(4, 0), -(1./256)*25);
  ASSERT_EQ(st2(5, 0), (1./256)*3);
}

class FixtureMLms : public::testing::Test {
 protected:
  int t, q, grid1, grid2, mc;
  double size, size_p, pressure, fineSize,
            coarseSize;
  std::vector<matrix> pfVec, cDVec, cCVec;
  std::vector<matrix> matrices;
  std::vector<matrix> pfVec1;
  void SetUp() override {
    size = 2;
    size_p = 1;
    pressure = 1.;

      // initial grid size for initialization
    grid1 = grids;
    grid2 = grids;

    fineSize = size / grids;

    matrix kM({grids, grids});
    matrix Ip({grids, grids});

    double lower_b = (grids)/2. - (size_p/fineSize)/2.;
    double upper_b = (grids)/2. + (size_p/fineSize)/2.;

    initializePressureArray(Ip, lower_b, upper_b, pressure);

    t = 5;
    mc = 2*t;

    matrix st = initializeStylusArray(t);

    initializeStack(st, t, Ip, kM, pfVec, cDVec);

    pfVec1 = pfVec;
    coarseSize = fineSize*pow(2, pfVec.size()-1);
    cCVec.reserve(3);
    matrices.push_back(st);
    matrices.push_back(Ip);
    matrices.push_back(kM);
  }
  void TearDown() override {}
};

TEST_F(FixtureMLms, calcCoarsePressure){
  GTEST_SKIP();
  matrix st = matrices[0];
  matrix Ip = matrices[1];
  matrix kM = matrices[2];
  createCorrectionArrays(cCVec, st, coarseSize, fineSize);

  calcCoarsePressure(pfVec1, st);
  writeToFile(pfVec1[pfVec.size()-1], "tests/new_press");
  std::cout << "\n there u go \n";
  old_calcCoarsePressure(pfVec, st);
  writeToFile(pfVec[pfVec.size()-1], "tests/old_press");
  bool equal = false;
  double eps = 0.00001;
  for (int i = 0; i < pfVec[pfVec.size()].shape[0]; i++) {
    for (int j = 0; j < pfVec[pfVec.size()].shape[1]; j++) {
      equal = std::abs(pfVec1[pfVec.size()](i, j) - pfVec[pfVec.size()](i, j)) < eps;
      /*if (equal == false) {
        std::cout << pfVec1[pfVec.size()](i, j) << " : " << pfVec[pfVec.size()](i, j) << 
          " , " << i << " : " << j <<std::endl;
        break;
      }*/
    }
  }
  EXPECT_TRUE(equal);
}

TEST_F(FixtureMLms, calculate) {
  // GTEST_SKIP();
  matrix st = matrices[0];
  matrix Ip = matrices[1];
  matrix kM = matrices[2];
  createCorrectionArrays(cCVec, st, coarseSize, fineSize);

  calcCoarsePressure(pfVec, st);

  int d = pfVec.size()-1;
  
  calc_displacement(pfVec[d], coarseSize, fineSize, cDVec[d]);
  // return;
  for (int i = 0; i < pfVec.size()-1; i++) {
    double hS = fineSize*pow(2, d-i-1);
    int temp_mc = (mc*2)+1;
    matrix cC({temp_mc, temp_mc});
    correctionSteps(cC, st, mc, t, fineSize, hS);
    applyCorrection(cDVec[d-i], cC, pfVec[d-i-1], t);
    interpolateGrid(cDVec[d-i-1], cDVec[d-i], st);
    secondCorrectionStep(st, hS,
                          pfVec[d-i-1], cDVec[d-i-1], cCVec);
  }
  resMlms_6 = cDVec[0];
  bool equal = false;
  bool equal1 = false;
  double eps = 0.05;
  for (int i = 0; i < cDVec[0].shape[0]; i++) {
    for (int j = 0; j < cDVec[0].shape[1]; j++) {
      equal = std::abs(cDVec[0](i, j) - cDVec[0](j, i)) < eps;
      equal1 = std::abs(res1(i, j) - cDVec[0](j, i)) < eps;
      /*if (equal == false) {
        std::cout << cDVec[0](i, j) << " : " << cDVec[0](j, i) << std::endl;
        break;
      } else if (equal1 == false) {
        std::cout << res1(i, j) << " : " << cDVec[0](j, i) << std::endl;
        break;
      }*/
    }
  }

  EXPECT_TRUE(equal);
  EXPECT_TRUE(equal1);
}

TEST(BoussinesqFFT2, calculate) {
  double Lx = 2., Ly = 2.;
  int Nx = grids, Ny = grids;
  double pSize = 1;
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
  // printarray(Gmn); return;
  transformGmnP(Nx, Ny, Gmn, Gmn_tild, p, p_tild);

  multiplyTransformed(Gmn_tild, Umn_tild, p_tild);

  transformToReal(Umn_tild, Umn, Nx, Ny);

  writeToResultArray(Umn, Umn_res, Nx, Ny);
  resFft = Umn_res;
  bool equal = false;
  bool equal1 = false;
  double eps = 0.001;
  for (int i = 0; i < Umn_res.shape[0]; i++) {
    for (int j = 0; j < Umn_res.shape[1]; j++) {
      equal = std::abs(Umn_res(i, j) - Umn_res(j, i)) < eps;
      equal1 = std::abs(res1(i, j) - Umn_res(i, j)) < eps;
      /*if (equal1 == false) {
        std::cout << res1(i, j) <<
                        " : " << Umn_res(i, j) <<
                        " : " << i << " : " << j <<
                        std::endl;
      }*/
    }
  }
  EXPECT_TRUE(equal);
  EXPECT_TRUE(equal1);
}

TEST(tsquare_norm, mlms_fft) {
  double size = 2;
  double size_p = 1;
  double pressure = 1.;
  double fineSize = size / grids;
  matrix kM({grids, grids});
  matrix Ip({grids, grids});
  double lower_b = (grids)/2. - (size_p/fineSize)/2.;
  double upper_b = (grids)/2. + (size_p/fineSize)/2.;
  initializePressureArray(Ip, lower_b, upper_b, pressure);
  int t = 2;
  int mc = 2*t;
  std::vector<matrix> pfVec, cDVec, cCVec;
  matrix st = initializeStylusArray(t);
  initializeStack(st, t, Ip, kM, pfVec, cDVec);
  double coarseSize = fineSize*pow(2, pfVec.size()-1);
  cCVec.reserve(3);
  createCorrectionArrays(cCVec, st, coarseSize, fineSize);
  calcCoarsePressure(pfVec, st);
  int d = pfVec.size()-1;
  calc_displacement(pfVec[d], coarseSize, fineSize, cDVec[d]);
  for (int i = 0; i < pfVec.size()-1; i++) {
    double hS = fineSize*pow(2, d-i-1);
    int temp_mc = (mc*2)+1;
    matrix cC({temp_mc, temp_mc});
    correctionSteps(cC, st, mc, t, fineSize, hS);
    applyCorrection(cDVec[d-i], cC, pfVec[d-i-1], t);
    interpolateGrid(cDVec[d-i-1], cDVec[d-i], st);
    secondCorrectionStep(st, hS,
                          pfVec[d-i-1], cDVec[d-i-1], cCVec);
  }

  double res_t2 = 0;
  for (int i = 0; i < cDVec[0].shape[0]; i++) {
    for (int j = 0; j < cDVec[0].shape[0]; j++) {
      res_t2 += std::abs(cDVec[0](i, j)-resFft(i, j));
    }
  }

  t = 3;
  mc = 2*t;
  std::vector<matrix> pfVec2, cDVec2, cCVec2;
  st = initializeStylusArray(t);
  initializeStack(st, t, Ip, kM, pfVec2, cDVec2);
  coarseSize = fineSize*pow(2, pfVec2.size()-1);
  cCVec2.reserve(3);
  createCorrectionArrays(cCVec2, st, coarseSize, fineSize);
  calcCoarsePressure(pfVec2, st);
  d = pfVec2.size()-1;
  calc_displacement(pfVec2[d], coarseSize, fineSize, cDVec2[d]);
  for (int i = 0; i < pfVec2.size()-1; i++) {
    double hS = fineSize*pow(2, d-i-1);
    int temp_mc = (mc*2)+1;
    matrix cC({temp_mc, temp_mc});
    correctionSteps(cC, st, mc, t, fineSize, hS);
    applyCorrection(cDVec2[d-i], cC, pfVec2[d-i-1], t);
    interpolateGrid(cDVec2[d-i-1], cDVec2[d-i], st);
    secondCorrectionStep(st, hS,
                          pfVec2[d-i-1], cDVec2[d-i-1], cCVec2);
  }

  double res_t3 = 0;
  for (int i = 0; i < cDVec2[0].shape[0]; i++) {
    for (int j = 0; j < cDVec2[0].shape[0]; j++) {
      res_t3 += std::abs(cDVec2[0](i, j)-resFft(i, j));
    }
  }

  
  double res_t4 = 0;
  for (int i = 0; i < cDVec2[0].shape[0]; i++) {
    for (int j = 0; j < cDVec2[0].shape[0]; j++) {
      res_t4 += std::abs(resMlms_6(i, j)-resFft(i, j));
    }
  }
  EXPECT_LT(res_t4, res_t3);
  EXPECT_LT(res_t3, res_t2);
  std::cout << res_t2 << " : " << res_t3 <<  " : " << res_t4 << std::endl;
}
