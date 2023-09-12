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
 TEMPLATE FOR TESTS
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
  double size_p = size/2;
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

TEST(Boussinesq, bCheck) {
  matrix ms({2, 2});
  int m = ms.shape[0];
  EXPECT_TRUE(bCheck(m, 1, 1));
  EXPECT_TRUE(bCheck(m, 1, 0));
  EXPECT_TRUE(bCheck(m, 0, 1));
  EXPECT_TRUE(bCheck(m, 0, 0));
  EXPECT_FALSE(bCheck(m, -1, 0));
  EXPECT_FALSE(bCheck(m, 1, -1));
  EXPECT_FALSE(bCheck(m, 3, 0));
  EXPECT_FALSE(bCheck(m, 0, 3));
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



TEST(TestMlms, calculate) {
  // GTEST_SKIP();

  resMlms_6 = BoussinesqMlms(2., grids, 3);
  bool equal = false;
  bool equal1 = false;
  double eps = 0.05;
  for (int i = 0; i < resMlms_6.shape[0]; i++) {
    for (int j = 0; j < resMlms_6.shape[1]; j++) {
      equal = std::abs(resMlms_6(i, j) - resMlms_6(j, i)) < eps;
      equal1 = std::abs(res1(i, j) - resMlms_6(j, i)) < eps;
      if (equal == false) {
        std::cout << resMlms_6(i, j) << " : " << resMlms_6(i, j) << std::endl;
        break;
      } else if (equal1 == false) {
        std::cout << res1(i, j) << " : " << resMlms_6(i, j) << std::endl;
        break;
      }
    }
  }
  writeToFile(resMlms_6, "results/BTest64");
  writeToFile(res1, "results/BTest64_naive");

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
      if (equal1 == false) {
        std::cout << res1(i, j) <<
                        " : " << Umn_res(i, j) <<
                        " : " << i << " : " << j <<
                        std::endl;
      }
    }
  }
  EXPECT_TRUE(equal);
  EXPECT_TRUE(equal1);
}

TEST(tsquare_norm, mlms_fft) {
  matrix rest1 = BoussinesqMlms(2., grids, 1);
  double res_t1 = 0;
  for (int i = 0; i < rest1.shape[0]; i++) {
    for (int j = 0; j < rest1.shape[0]; j++) {
      res_t1 += std::abs(rest1(i, j)-resFft(i, j));
    }
  }

  matrix rest2 = BoussinesqMlms(2., grids, 2);
  double res_t2 = 0;
  for (int i = 0; i < rest2.shape[0]; i++) {
    for (int j = 0; j < rest2.shape[0]; j++) {
      res_t2 += std::abs(rest2(i, j)-resFft(i, j));
    }
  }

  double res_t3 = 0;
  for (int i = 0; i < rest2.shape[0]; i++) {
    for (int j = 0; j < rest2.shape[0]; j++) {
      res_t3 += std::abs(resMlms_6(i, j)-resFft(i, j));
    }
  }
  EXPECT_LT(res_t3, res_t2);
  EXPECT_LT(res_t2, res_t1);
  std::cout << res_t1 << " : " << res_t2 <<  " : " << res_t3 << std::endl;
}
