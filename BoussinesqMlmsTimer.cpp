/* Copyright 2023 <Ilia Fedotov @ Uni Freiburg> */
#include "BoussinesqMlmsTimer.h"
#include "BoussinesqMlms.h"
#include "Boussinesq.h"
#include "BoussinesqTimer.h"
#include <algorithm>
#include <vector>
#include <cmath>

void runFastTimerLoops() {
  const double size = 2;
  const double size_p = 1;
  const double pressure = 1.;

  // initial grid size for initialization
  int grid1 = 2048;

  matrix result({100, 2});
  for (int iteration = 0; iteration < 2; iteration++) {
    double fineSizeA = size / grid1;

    matrix kM({grid1, grid1});
    matrix Ip({grid1, grid1});

    double lower_b = (grid1)/2. - (size_p/fineSizeA)/2.;
    double upper_b = (grid1)/2. + (size_p/fineSizeA)/2.;

    initializePressureArray(Ip, lower_b, upper_b, pressure);
    // material moduli and v

    int t = 4;
    int mc = 2*t;
    
    matrix st = initializeStylusArray(t);

    std::vector<matrix> pfVec;
    std::vector<matrix> cDVec;

    std::vector<int> qs;

    initializeStack(st, t, Ip, kM, grid1, pfVec, cDVec);

    double d = qs.size();

    double coarseSize = fineSizeA*pow(2, d);
    std::vector<matrix> cCVec;
    cCVec.reserve(3);
    createCorrectionArrays(cCVec, st, coarseSize, fineSizeA);

    MeasureTime mt = MeasureTime();
    calcCoarsePressure(pfVec, st);

    // calc_displacement(pfVec[d], coarseSize, fineSizeA, cDVec[d]);
    naiveCalculation(cDVec[d], pfVec[d], coarseSize);
    
    for (int i = 0; i < qs.size(); i++) {
      double hS = fineSizeA*pow(2, d-i-1);
      int temp_mc = (mc*2)+1;
      matrix cC({temp_mc, temp_mc});
      correctionSteps(cC, st, mc, t, fineSizeA, hS);
      applyCorrection(cDVec[d-i], cC, pfVec[d-i-1], t);
      interpolateGrid(cDVec[d-i-1], cDVec[d-i], st);
      secondCorrectionStep(st, hS,
                           pfVec[d-i-1], cDVec[d-i-1], cCVec);
    }
    auto stopped = mt.stopTime();

    result(iteration, 0) = grid1;
    result(iteration, 1) = stopped;
    writeToFile(cDVec[0], "./grid_faster_"+std::to_string(grid1));
    grid1 += 20;
  }
  writeToFile(result, "timings_mlms");
}
