/* Copyright 2023 <Ilia Fedotov @ Uni Freiburg> */

#include <gtest/gtest.h>
#include "./Mlms.h"

/*TEST(filename, function){
  function values
  ASSERT_EQ(res, function call);
}*/

TEST(MLMS, precalculation) {
  const double PI = 3.141592653589793238463;
  ASSERT_EQ(0, precalculation(5., 5., 5., 5., 1., 1., 1.));
}
