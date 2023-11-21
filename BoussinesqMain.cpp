/* Copyright 2023 <Ilia Fedotov @ Uni Freiburg> */

#include "./Boussinesq.h"
#include "./BoussinesqMlms.h"
#include "./BoussinesqFft.h"
#include "./Boussinesq.h"

int main() {
  int ct = 0;
  std::vector<int> n{7,15,29,31,35,41,49,53,57};
  // for (int i = 31; i < 260; i*=2) {
  for (int i = 0; i < n.size(); i++) {
    // replace 0.84*log(i) for constant t
    // matrix res_t2 = BoussinesqMlms(2., i, 0.84*log(i));
    matrix res_t2 = BoussinesqFFT(2., n[i]);
    writeToFile(res_t2, "mlms"+std::to_string(n[i]));
    ct++;
  }


  return 0;
}