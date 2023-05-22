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

double grids = 8;
matrix res1({grids, grids});

TEST(BoussinesqNaive, calculate) {
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
  naiveCalculation(Ic, Pa, cell_size);

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
  int tl = 2;
  matrix st({2*tl, 1});
  initializeStylusArray(st, tl);
  // printarray(st);

  int t2 = 4;
  matrix st1({2*t2, 1});
  // std::cout << "\n\n" << std::endl;

  initializeStylusArray(st1, t2);
  // printarray(st1);
  // std::cout << "\n\n" << std::endl;
  /*ASSERT_EQ(st1(0,0), (-1./16.) );
  std::cout << st1.shape[0] << std::endl;
  // ASSERT_EQ(st1(1,0), 0 );
  ASSERT_EQ(st1(2,0), (1./16.)*9 );
  ASSERT_EQ(st1(3,0), (1/16.)*16 );
  ASSERT_EQ(st1(4,0), (1/16.)*16 );
  ASSERT_EQ(st1(5,0), (1/16.)*9);
  ASSERT_EQ(st1(6,0), 0);
  ASSERT_EQ(st1(7,1), -(1/16.)*1);*/

  tl = 6;
  matrix st2({2*tl, 1});
  initializeStylusArray(st2, tl);
  // printarray(st2);
  // std::cout << "\n\n" << std::endl;
  /*ASSERT_EQ(st2(0,0), -(1/256)*3);
  ASSERT_EQ(st2(1,0), 0 );
  ASSERT_EQ(st2(2,0), -(1/256)*25);
  ASSERT_EQ(st2(3,0), 0 );
  ASSERT_EQ(st2(4,0), (1/256)*150);
  ASSERT_EQ(st2(5,0), (1/256)*256);
  ASSERT_EQ(st2(6,0), (1/256)*256);
  ASSERT_EQ(st2(7,0), (1/256)*150);
  ASSERT_EQ(st2(8,0), 0);
  ASSERT_EQ(st2(9,0), -(1/256)*25);
  ASSERT_EQ(st2(10,0), 0);
  ASSERT_EQ(st2(11,0), (1/256)*3);*/


  // printarray(st2);
}

/*TEST(BoussinesqMlms, calcCoarsePressure) {
  std::vector<int>& qs;
  double qLevel = std::log2(grid * grid)/2-1;
  int q = grid / 2 + 2*t - 1;
  qs.push_back(q);

  for (int length = 1; length < qLevel-1; length++) {
    q = q / 2 + 2*t - 1;
    qs.push_back(q);
  }
  std::vector<matrix>& pFVec;
  std::vector<matrix>& cDVec;
  int t = 2; 
  const matrix& st(t,1);
}*/

/*TEST(BoussinesqMlms, correctionSteps) {}
TEST(BoussinesqMlms, applyCorrection) {}
TEST(BoussinesqMlms, correctionHelper) {}
TEST(BoussinesqMlms, interpolateGrid) {}

TEST(BoussinesqMlms, secondCorrectionStep) {}
TEST(BoussinesqMlms, createCorrectionArrays) {}
TEST(BoussinesqMlms, interpolateGrid) {}
TEST(BoussinesqMlms, interpolateGrid) {}*/

class FixtureMLms : public::testing::Test {
  protected:
    void SetUp() override {
      size = 2;
      size_p = 1;
      pressure = 1.;

      // initial grid size for initialization
      grid1 = grids;
      grid2 = grids;

      fineSizeA = size / grid1;
      fineSizeB = size / grid2;

      matrix kM({grid1, grid2});
      matrix Ip({grid1, grid2});

      double lower_b = (grid1)/2. - (size_p/fineSizeA)/2.;
      double upper_b = (grid2)/2. + (size_p/fineSizeB)/2.;

      initializePressureArray(Ip, lower_b, upper_b, pressure);
      // material moduli and v

      double beta = 0.84;
      double min_g = std::min(grid1, grid2);
      // int t = beta*log(min_g);
      t = 5;
      mc = 0.7* pow(min_g, 1./t)-1;
      if (mc < 2*t) {
        mc = 2*t;
      }
      // std::cout << "\”\n" << std::endl;
      matrix st({2*t, 1});

      initializeStylusArray(st, 1);

      double qLevel = std::log2(grid1 * grid2)/2-1;
      q = grid1 / 2 + 2*t - 1;
      qs.push_back(q);

      for (int length = 1; length < qLevel-1; length++) {
        q = q / 2 +  2*t - 1;
        qs.push_back(q);
      }

      d = qs.size();
      pfVec.reserve(d);
      cDVec.reserve(d);
      pfVec.push_back(Ip);
      cDVec.push_back(kM);
      coarseSize = fineSizeA*pow(2, d);
      cCVec.reserve(3);
    }
    using matrix = matrixTemplate<double>;
    virtual void TearDown() {}
    int t, q, grid1, grid2, mc;
    double size, size_p, pressure, fineSizeA,
            fineSizeB, coarseSize, d;
    std::vector<matrix> pfVec, cDVec, cCVec;
    matrix st, kM, Ip;
    std::vector<int> qs;
};

/*TEST_F(FixtureMLms, calculate) {
  /*const double size = 2;
  const double size_p = 1;
  const double pressure = 1.;

  // initial grid size for initialization
  int grid1 = grids;
  int grid2 = grids;

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
  // std::cout << "\”\n" << std::endl;
  matrix st({2*t, 1});

  initializeStylusArray(st, 1);

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
  cCVec.reserve(3);*//*
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
  double eps = 0.01;
  for (int i = 0; i < cDVec[0].shape[0]; i++) {
    for (int j = 0; j < cDVec[0].shape[1]; j++) {
      equal = std::abs(cDVec[0](i, j) - cDVec[0](j, i)) < eps;
      equal1 = std::abs(res1(i, j) - cDVec[0](j, i)) < eps;
      if (equal == false) {
        std::cout << cDVec[0](i, j) << " : " << cDVec[0](j, i) << std::endl;
        break;
      } else if (equal1 == false) {
        std::cout << res1(i, j) << " : " << cDVec[0](j, i) << std::endl;
        break;
      }
    }
  }
  EXPECT_TRUE(equal);
  EXPECT_TRUE(equal1);
}*/

TEST(BoussinesqFFT, calculateGmn) {
  double Lx = 2., Ly = 2.;
  int Nx = grids, Ny = grids;
  double pSize = 1;
  double dx = (Lx/Nx);
  double dy = (Ly/Ny);

  matrix Gmn({(2*Nx)-1, (2*Ny)-1});
  calculateGmn(Gmn, dx, dy);
  printarray(Gmn);
}

TEST(BoussinesqFFT2, calculate) {
  double Lx = 2., Ly = 2.;
  int Nx = grids, Ny = grids;
  double pSize = 1;
  double dx = (Lx/Nx);
  double dy = (Ly/Ny);

  int lb = Nx/2-(pSize/dx)/2;
  int ub = Ny/2+(pSize/dy)/2;

  matrix tempP({Nx, Ny});
  matrix p({2*Nx, 2*Ny});
  cMatrix p_tild({2*Nx, (2*Ny)/2-1});
  matrix Gmn({2*Nx, 2*Ny});
  cMatrix Gmn_tild({2*Nx, (2*Ny)/2-1});
  matrix Umn({Nx*2, Ny*2});
  cMatrix Umn_tild({2*Nx, (2*Ny)/2-1});
  matrix Umn_res({Nx, Ny});

  initializePressureArray(tempP, lb, ub, 1.);
  initializeDisplacementArray(p);

  copyPressureArray(p, tempP);
  printarray(p);

  calculateGmn(Gmn, dx, dy);

  transformGmnP(Nx, Ny, Gmn, Gmn_tild, p, p_tild);

  multiplyTransformed(Gmn_tild, Umn_tild, p_tild);

  transformToReal(Umn_tild, Umn, Nx, Ny);

  writeToResultArray(Umn, Umn_res, Nx, Ny);

  bool equal = false;
  bool equal1 = false;
  double eps = 0.001;
  for (int i = 0; i < Umn_res.shape[0]; i++) {
    for (int j = 0; j < Umn_res.shape[1]; j++) {
      equal = std::abs(Umn_res(i, j) - Umn_res(j, i)) < eps;
      equal1 = std::abs(res1(i, j) - Umn_res(i, j)) < eps;
      if (equal1 == false) {
        std::cout << res1(i, j) << 
                        " : " << Umn_res(i, j) <<
                        " : " << i << " : " << j <<
                        std::endl;
        // break;
      }
    }
  }
  EXPECT_TRUE(equal);
  EXPECT_TRUE(equal1);
}
