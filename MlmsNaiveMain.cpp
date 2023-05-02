/* Copyright 2023 <Ilia Fedotov @ Uni Freiburg> */
#include <vector>
#include <iostream>
#include <cmath>
#include "./MlmsWriteToFile.h"
#include "./Mlms.h"
#include "./MlmsTimer.h"
#include "./MlmsFastTimer.h"

int main() {
  // initial size and pressure values
  runFastTimerLoops();
  runTimerLoops();
  return 0;
}

