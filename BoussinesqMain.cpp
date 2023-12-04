/* Copyright 2023 <Ilia Fedotov @ Uni Freiburg> */

#include "./Boussinesq.h"
#include "./BoussinesqMlms.h"
#include "./BoussinesqFft.h"
#include "./Boussinesq.h"

int main() {
  std::vector<int> n{15,29,35,41,53,57};
  // std::vector<int> n{16,20,24,32,38,48,64};
  // for (int i = 31; i < 260; i*=2) {
  for (std::size_t i = 0; i < n.size(); i++) {
    // replace 0.84*log(i) for constant t
    // matrix res_t2 = BoussinesqMlms(2., i, 0.84*log(i));
    matrix res_t2 = BoussinesqFFT(2., n[i]);
    writeToFile(res_t2, "mlms"+std::to_string(n[i]));
    res_t2.setZero();
  }


  return 0;
}