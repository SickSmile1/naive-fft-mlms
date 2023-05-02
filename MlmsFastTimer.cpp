/* Copyright 2023 <Ilia Fedotov @ Uni Freiburg> */
#include "MlmsTimer.h"
#include "cmath"
#include "Mlms.h"
#include "MlmsWriteToFile.h"

void runFastTimerLoops() {
  const double size = 2;
  const double size_p = 1;
  const double pressure = 1.;

  // initial grid size for initialization
  int grid1 = 60;
  int grid2 = 60;
  matrix result({100, 2});
  for (int iteration = 0; iteration < 200; iteration++) {
    double fineSizeA = size / grid1;
    double fineSizeB = size / grid2;

    matrix kM({grid1, grid2});
    matrix Ip({grid1, grid2});

    double lower_b = (grid1)/2. - (size_p/fineSizeA)/2.;
    double upper_b = (grid2)/2. + (size_p/fineSizeB)/2.;

    initializePressureArray(Ip, lower_b, upper_b, pressure);
    // material moduli and v
    double v = 0.;
    double E = 1.;

    double beta = 0.84;
    double min_g = std::min(grid1, grid2);
    // int t = beta*log(min_g);
    int t = 4;
    int mc = 0.7* pow(min_g, 1./t)-1;
    if (mc < 2*t) {
      mc = 2*t;
    }
    // std::cout << mc << std::endl;
    matrix st({2*t, 1});

    initializeStylusArray(st, t);

    std::vector<matrix> pfVec;
    std::vector<matrix> cDVec;

    std::vector<int> qs;

    double qLevel = std::log2(grid1 * grid2)/2-1;
    int q = grid1 / 2 + 2*t - 1;
    qs.push_back(q);

    for (int length = 1; length < qLevel-1; length++) {
      q = q / 2 + 2*t - 1;
      qs.push_back(q);
    }

    double d = qs.size();
    pfVec.reserve(d);
    cDVec.reserve(d);
    pfVec.push_back(Ip);
    cDVec.push_back(kM);
    double coarseSize = fineSizeA*pow(2, d);
    std::vector<matrix> cCVec;
    cCVec.reserve(3);
    createCorrectionArrays(cCVec, st, coarseSize, fineSizeA,
                           fineSizeB, mc);

    MeasureTime mt = MeasureTime();
    calcCoarsePressure(qs, pfVec, cDVec, t, st);
    
    calc_displacement(pfVec[d], coarseSize, fineSizeA, cDVec[d]);

    for (int i = 0; i < qs.size(); i++) {
      double hS = fineSizeA*pow(2, d-i-1);
      int temp_mc = (mc*2)+1;
      matrix cC({temp_mc, temp_mc});
      correctionSteps(cC, st, mc, t, fineSizeA, fineSizeB, hS);
      // return 0;
      applyCorrection(cDVec[d-i], cC, pfVec[d-i-1], t);
      interpolateGrid(cDVec[d-i-1], cDVec[d-i], st);
      secondCorrectionStep(mc, st, fineSizeA, fineSizeB, hS,
                           pfVec[d-i-1], cDVec[d-i-1], cCVec);
      // return 0;
    }
    auto stopped = mt.stopTime();

    result(iteration, 0) = grid1;
    result(iteration, 1) = stopped;
    // writeToFile(cDVec[0], "./test_mlms/grid_faster_"+std::to_string(grid1));
    grid1 += 20;
    grid2 += 20;
  }
  writeToFile(result, "timings_mlms");
}
