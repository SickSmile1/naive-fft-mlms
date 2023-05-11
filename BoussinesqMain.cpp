/* Copyright 2023 <Ilia Fedotov @ Uni Freiburg> */

#include "./Boussinesq.h"
// #include "./MlmsTimer.h"
// #include "./MlmsFastTimer.h"
#include "./BoussinesqFftTimer.h"

int main() {
  // initial size and pressure values
  // runFastTimerLoops();
  // runTimerLoops();
  BoussinesqFFT();
  return 0;
}

