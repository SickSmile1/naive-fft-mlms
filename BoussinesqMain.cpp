/* Copyright 2023 <Ilia Fedotov @ Uni Freiburg> */

#include "./Boussinesq.h"
#include "./BoussinesqMlms.h"
#include "./BoussinesqFft.h"
#include "./Boussinesq.h"

int main() {

  matrix res({9, 3});
  int ct = 0;
  for (int i = 2048; i < 2200; i*=2) {
    // replace 0.84*log(i) for constant t
    matrix res_t2 = BoussinesqMlms(2., i, 0.84*log(i));
    ct++;
  }

  // writeToFile(res, "results/res_l2_const-t-8");
  return 0;
}

