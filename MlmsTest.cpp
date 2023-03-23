/* Copyright 2023 <Ilia Fedotov @ Uni Freiburg> */

#include <gtest/gtest.h>
#include <iostream>
#include <fstream>
#include <thread> //NOLINT
#include <string>
#include "./Mlms.h"
#include "./MlmsWriteToFile.h"
#include "./MlmsTimer.h"

/*TEST(filename, function){
  function values
  ASSERT_EQ(res, function call);
}*/

TEST(MLMS, calculate) {
  // test for symmetry of resulting matrix
  double size = 200;  // default: 200

  std::size_t size_p = 100;  // default: 100

  std::size_t grid = 20;  // default: 50

  double cell_size = size/grid;

  double v = 0;  // default: 0
  double E = 1;  // default: 1

  matrix Ic({grid, grid});
  matrix Pa({grid, grid});

  double pressure = 1.;
  double lower_b = (size/cell_size)/2 - (size_p/cell_size)/2;
  double upper_b = (size/cell_size)/2 + (size_p/cell_size)/2;

  initializePressureArray(Pa, lower_b, upper_b, pressure);
  initializeDisplacementArray(Ic);
  calculation_loop(Ic, Pa, cell_size, v, E);

  bool equal = false;
  double eps = 0.001;
  for (std::size_t i = 0; i < Ic.shape[0]; i++) {
    for (std::size_t j = 0; j < Ic.shape[1]; j++) {
      // std::cout << Ic(i, j) << " : " << Ic(j, i) << std::endl;
      equal = std::abs(Ic(i, j) - Ic(j, i)) < eps;
      if (equal == false) break;
    }
  }

  EXPECT_TRUE(equal);
}

TEST(MlmsTimer, timeTest) {
  MeasureTime mt = MeasureTime();
  std::this_thread::sleep_for(std::chrono::milliseconds(1000));
  auto stopped =  mt.stopTime();
  EXPECT_GE(stopped, 1);
}

TEST(MlmsWriteToFile, writeToFile) {
  // write/read file and check if values match
  matrix test({3, 3});
  int counter = 0;
  for (std::size_t i = 0; i < test.shape[0]; i++) {
    for (std::size_t j = 0; j < test.shape[1]; j++) {
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
