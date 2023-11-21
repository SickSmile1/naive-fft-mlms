/* Copyright 2023 <Ilia Fedotov @ Uni Freiburg> */

#include <Eigen/Core>
#include <gtest/gtest.h>
#include <iostream>
#include <fstream>
#include <thread> //NOLINT
#include <string>
#include "Boussinesq.h"
#include "BoussinesqMlms.h"
#include "BoussinesqFft.h"
#pragma GCC optimize ("no-fast-math")
// Template
/*TEST(filename, function){
 TEMPLATE FOR TESTS
  function values
  ASSERT_EQ(res, function call);
}*/

int grids = 31;
matrix res1({grids, grids});
matrix resMlms_6({grids, grids});
matrix resFft({grids, grids});

TEST(BoussinesqNaive, calculate) {
  res1.setZero();
  // GTEST_SKIP();
  // test for symmetry of resulting matrix
  double size = 2.;
  double size_p = size/2;
  double cell_size = size/grids;

  matrix Ic({grids, grids});
  matrix Pa({grids, grids});

  double pressure = 1.;
  double lower_b = (size/cell_size)/2 - (size_p/cell_size)/2;
  double upper_b = (size/cell_size)/2 + (size_p/cell_size)/2;

  initializePressureArray(Pa, lower_b, upper_b, pressure);
  initializeDisplacementArray(Ic);
  calc_displacement(Pa, cell_size, cell_size, Ic);

  bool equal = false;
  double eps = 0.001;
  for (int i = 0; i < Ic.rows(); i++) {
    for (int j = 0; j < Ic.cols(); j++) {
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

TEST(KernelNotNan, BoussinesqTest) {
  GTEST_SKIP();
  matrix test({31,31});
  test.setZero();
  matrix press =  Eigen::MatrixXd::Ones(31,31);
  press.leftCols(2) = test.leftCols(2);
  press.rightCols(2) = test.rightCols(2);
  press.topRows(2) = test.topRows(2);
  press.bottomRows(2) = test.bottomRows(2);
  auto dxy = Eigen::ArrayXd::LinSpaced(100,-.02,.02);
  // dxy -= Eigen::ArrayXd::Ones(100);
  for(auto dx : dxy ) {
    calc_displacement(press, dx, dx, test);
    double s = test.sum();
    ASSERT_GT(s,0);
    bool nan, nan2;
    for(auto v : test.reshaped()) {
      nan = std::isnan(v);
      nan2 = (v!=v);
    }
    ASSERT_FALSE(nan);
    ASSERT_FALSE(nan2);
    test.setZero();
  }
}

TEST(WriteToFile, writeToFile) {
  // write/read file and check if values match,
  // --> works, not needed for algorithm changes
  GTEST_SKIP();
  matrix test({3, 3});
  int counter = 0;
  for (int i = 0; i < test.rows(); i++) {
    for (int j = 0; j < test.cols(); j++) {
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
  // test out of bounds function
  GTEST_SKIP();
  matrix ms({2, 2});
  int m = ms.rows();
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
  // Test stylus array as described in Lubrecht & Brandt 20a, 20b, 20c
  GTEST_SKIP();
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
  GTEST_SKIP();

  resMlms_6 = BoussinesqMlms(2., grids, 3);
  bool equal = false;
  bool equal1 = false;
  double eps = 0.05;
  for (int i = 0; i < resMlms_6.rows(); i++) {
    for (int j = 0; j < resMlms_6.cols(); j++) {
      equal = std::abs(resMlms_6(i, j) - resMlms_6(j, i)) < eps;
      equal1 = std::abs(res1(i, j) - resMlms_6(i, j)) < eps;
      if (equal == false) {
        std::cout << resMlms_6(i, j) << " : " << resMlms_6(i, j) << std::endl;
        break;
      } else if (equal1 == false) {
        std::cout << res1(i, j) << " : " << resMlms_6(i, j) << std::endl;
        break;
      }
    }
  }
  // writeToFile(resMlms_6, "results/BTest64");
  // writeToFile(res1, "results/BTest64_naive");

  EXPECT_TRUE(equal);
  EXPECT_TRUE(equal1);
}

TEST(BoussinesqFFT2, calculate) {
  GTEST_SKIP();
  double Lx = 2., Ly = 2.;
  int Nx = grids, Ny = grids;
  double pSize = 1;
  double dx = (Lx/Nx);
  double dy = (Ly/Ny);

  int lb = Nx/2-(pSize/dx)/2;
  int ub = Ny/2+(pSize/dy)/2;
  std::cout << "dis is dx: " << dx << std::endl;
  // std::cout << "dis is dx: " << dx << std::endl;

  matrix Gmn({(2*Nx)-1, (2*Ny)-1});
  cMatrix Gmn_tild({Gmn.rows(), Gmn.cols()/2+1});
  matrix p({Gmn.rows(), Gmn.cols()});
  cMatrix p_tild({Gmn.rows(), Gmn.cols()/2+1});
  matrix tempP({Nx, Ny});
  matrix Umn({Gmn.rows(), Gmn.cols()});
  cMatrix Umn_tild({Gmn.rows(), Gmn.cols()/2+1});
  matrix Umn_res({Nx, Ny});

  initializePressureArray(tempP, lb, ub, 1.);
  initializeDisplacementArray(p);
  copyPressureArray(p, tempP);
  calculateGmn(Gmn, dx, dy);

  transformGmnP(Gmn, Gmn_tild, p, p_tild);

  multiplyTransformed(Gmn_tild, Umn_tild, p_tild);

  transformToReal(Umn_tild, Umn);

  writeToResultArray(Umn, Umn_res);
  resFft.setZero();
  resFft = Umn_res;
  bool equal = false;
  bool equal1 = false;
  double eps = 0.01;
  for (int i = 0; i < Umn_res.rows(); i++) {
    for (int j = 0; j < Umn_res.cols(); j++) {
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

TEST(BoussinesqFft3, fftprime) {
  std::vector<int> n{7,15,29,31,35,41,49,53,57};
  // for (int i = 31; i < 260; i*=2) {
  for (std::size_t k = 0; k < n.size(); k++) {
    // replace 0.84*log(i) for constant t
    // matrix res_t2 = BoussinesqMlms(2., i, 0.84*log(i));
    matrix res_t2 = BoussinesqFFT(2., n[k]);
    double eps = 0.000001;
    bool equal = false;
    for (long int i = 0; i < res_t2.rows(); i++) {
      for (long int j = 0; j < res_t2.cols(); j++) {
        equal = std::abs(res_t2(i, j) - res_t2(j, i)) < eps;
        /*if (equal == false) {
          std::cout << res1(i, j) <<
                          " : " << res_t2(i, j) <<
                          " : " << i << " : " << j <<
                          std::endl;
        }*/
      }
    }
    EXPECT_TRUE(equal);
  }
}

TEST(tsquare_norm, mlms_fft) {
  double fftSum = resFft.sum();
  matrix rest1 = BoussinesqMlms(2., grids, 1);
  double res_t1 = std::abs(rest1.sum()-fftSum);
  /* for (int i = 0; i < rest1.rows(); i++) {
    for (int j = 0; j < rest1.cols(); j++) {
      res_t1 += std::abs(rest1(i, j)-resFft(i, j));
    }
  } */
  writeToFile(rest1, "mlms1");

  matrix rest2 = BoussinesqMlms(2., grids, 2);
  double res_t2 = std::abs(rest2.sum()-fftSum);
/*   for (int i = 0; i < rest2.rows(); i++) {
    for (int j = 0; j < rest2.cols(); j++) {
      res_t2 += std::abs(rest2(i, j)-resFft(i, j));
    }
  } */
  writeToFile(rest2, "mlms2");

  resMlms_6.setZero();
  resMlms_6 = BoussinesqMlms(2.,grids,3);
  double res_t3 = std::abs(resMlms_6.sum()-fftSum);
/*   for (int i = 0; i < rest2.rows(); i++) {
    for (int j = 0; j < rest2.cols(); j++) {
      res_t3 += std::abs(resMlms_6(i, j)-resFft(i, j));
    }
  } */
  writeToFile(resMlms_6, "mlms3");
  EXPECT_LT(res_t3, res_t2);
  EXPECT_LT(res_t2, res_t1);
  std::cout << res_t1 << " : " << res_t2 <<  " : " << res_t3 << std::endl;
}
