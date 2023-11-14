/* Copyright 2023 <Ilia Fedotov @ Uni Freiburg> */

#include "./Boussinesq.h"
#include "./BoussinesqMlms.h"
#include "./BoussinesqFft.h"
#include "./Boussinesq.h"

int main() {
  int ct = 0;
  for (int i = 16; i < 270; i*=2) {
    // replace 0.84*log(i) for constant t
    matrix res_t2 = BoussinesqMlms(2., i, 0.84*log(i));
    writeToFile(res_t2, "mlms"+std::to_string(i));
    ct++;
  }


  return 0;
}

