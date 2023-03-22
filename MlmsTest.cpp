/* Copyright 2023 <Ilia Fedotov @ Uni Freiburg> */

#include <gtest/gtest.h>
#include "./Mlms.h"
#include "./MlmsWriteToFile.h"

/*TEST(filename, function){
  function values
  ASSERT_EQ(res, function call);
}*/

TEST(MLMS, precalculation) {
  const double PI = 3.141592653589793238463;
  double res;
  double res1 = 129.27107860874096;
  int grid_size = 20;
  double cell_size = 10;
  for (int i = 0; i < grid_size*grid_size; i++) {
    res += precalculation(5, 5, (i%grid_size)*cell_size,
                          (i/grid_size)*cell_size);
  }
  res = res *(1/(PI));
  ASSERT_DOUBLE_EQ(res, res1);
}

TEST(MlmsWriteToFile, writeTimeToFile) {
  ASSERT_DOUBLE_EQ(0, 0);
}

TEST(MlmsWriteToFile, writeVecToFile) {
  ASSERT_DOUBLE_EQ(0, 0);
}
