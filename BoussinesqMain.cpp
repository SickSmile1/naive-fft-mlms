/* Copyright 2023 <Ilia Fedotov @ Uni Freiburg> */

#include "./Boussinesq.h"
// #include "./MlmsTimer.h"
#include "./BoussinesqMlms.h"
#include "./BoussinesqFft.h"
#include "./Boussinesq.h"

int main() {
  // initial size and pressure values
  // runFastTimerLoops();
  // runTimerLoops();
  // BoussinesqFFT();

  auto means = [](double &min, double &max, double &mean, const matrix v1, const matrix v2) {
    double res1 = 0;
    for (int i = 0; i < v1.shape[0]; i++) {
      for (int j = 0; j < v1.shape[0]; j++) {
        double temp = std::pow(v1(i, j)-v2(i, j), 2);
        if (temp < min) min = temp;
        if (temp > max) max = temp;
        res1 += temp;
      }
    }
    double n = ((v1.shape[0]-1),2);
    mean = std::sqrt(res1/n);
  };

  matrix res({9, 3});
  int ct = 0;
  for (int i = 8; i < 2100; i*=2) {
    matrix res_t2 = BoussinesqMlms(2., i, 0.84*log(i));
    //matrix res_t3 = BoussinesqMlms(2., i, 3);
    // matrix res_t4 = BoussinesqMlms(2., i, 4);
    matrix res_fft = BoussinesqFFT(2., i+1);
    double min_0 = std::numeric_limits<double>::min();
    double max_0 = std::numeric_limits<double>::max();
    double min = max_0, max = min_0, mean = 0;
    means(min, max, mean, res_t2, res_fft);
    res(ct, 0) = min;
    res(ct, 1) = max;
    res(ct, 2) = mean;
    /*min = max_0, max = min_0, mean = 0;
    means(min, max, mean, res_t3, res_fft);
    res(ct, 3) = min;
    res(ct, 4) = max;
    res(ct, 5) = mean;
    min = max_0, max = min_0, mean = 0;
    means(min, max, mean, res_t4, res_fft);
    res(ct, 6) = min;
    res(ct, 7) = max;
    res(ct, 8) = mean;*/
    ct++;
  }

  // matrix res_t2 = BoussinesqMlms(2., 512, 2);
  // writeToFile(res_t2, "results/res_mlms2048");
  // writeToFile(BoussinesqFFT(2.,2048), "results/res_fft3");
  writeToFile(res, "results/res_l2_corr");
  return 0;
}

