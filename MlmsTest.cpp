/* Copyright 2023 <Ilia Fedotov @ Uni Freiburg> */

#include <gtest/gtest.h>
#include "./Mlms.h"

/*TEST(filename, function){
  function values
  ASSERT_EQ(res, function call);
}*/

TEST(MLMS, precalculation) {
  double res = 16.731751079260256;
  ASSERT_DOUBLE_EQ(res, precalculation(7.5, 7.5, 1., 1.));
}
