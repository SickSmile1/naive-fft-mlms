/* Copyright 2023 <Ilia Fedotov @ Uni Freiburg> */
#include <vector>
#include <iostream>
#include <cmath>
#include "./MlmsWriteToFile.h"
#include "./Mlms.h"

int main() {
  const double PI = 3.141592653589793238463;

  // initial size and pressure values
  const double size = 2;
  const double size_p = 1;
  double pressure = 1.;

  // initial grid size for initialization
  int grid1 = 2048;
  int grid2 = 2048;

  double fineSizeA = size / grid1;
  double fineSizeB = size / grid2;

  matrix kM({grid1, grid2});
  matrix Ip({grid1, grid2});

  double lower_b = (grid1)/2. - (size_p/fineSizeA)/2.;
  double upper_b = (grid2)/2. + (size_p/fineSizeB)/2.;

  initializePressureArray(Ip, lower_b, upper_b, pressure);
  writeToFile(Ip, "finePressure");
  // material moduli and v
  double v = 0.;
  double E = 1.;

  double beta = 0.84;
  double min_g = std::min(grid1, grid2);
  //int t = beta*log(min_g);
  int t = 4;
  // double mc = 0.7* pow(min_g, 1./t)-1;

  //if (mc < 2*t) {
    int mc = 2*t;
  //}

  matrix st({2*t, 1});
  
  initializeStylusArray(st, t);

  // for (auto &a: st) std::cout << a << std::endl;
  std::vector<matrix> pfVec;
  std::vector<matrix> cDVec;
  
  std::vector<double> qs;
  
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
  // std::cout << "halfSize: "<< coarseSize << std::endl;
  // std::cout << "d:" << d << " gridSize[5]" << gridSize[d-1] << std::endl;
  // std::cout << fineSizeA << " finesizea\n"; 
  
  calcCoarsePressure(qs, pfVec, cDVec, t, st);
  // for (auto& e: qs) { std::cout << e << std::endl; }
  std::cout << d << std::endl;
  int shape = pfVec[d].shape[0];
  matrix coarseDisplacement({shape, shape});
  for (int i = 0; i < shape; ++i) {
    for (int j = 0; j < shape; ++j) {
      for (int k = 0; k < shape; ++k) {
        for (int l = 0; l < shape; ++l) {
          cDVec[d](i, j) += calculate(i-k, j-l,  coarseSize, coarseSize,
                          fineSizeA,
                          fineSizeB)*pfVec[d](k, l); //* ((1-pow(v, 2))/(PI*E)) * pF(i, j);
        }
      }
    }
  }
  writeToFile(cDVec[d], "./tests/ds");
  // return 0;
  // printarray(coarseDisplacement);
  
    std::cout << d << " <- thats the d\n";
  for (int i = 0; i < qs.size(); i++) {
    myBreakpoint();
    std::cout << cDVec[d-i-1].shape[0]-1 << " thats the size of displacement array\n";
    double hS = fineSizeA*pow(2, d-i-1);
    std::cout << hS << " : halfsize \n";
    int temp_mc = (mc*2)+1;
    matrix cC({temp_mc, temp_mc});
    correctionSteps(cC, st, mc, t, fineSizeA, fineSizeB, hS);
    writeToFile(cC, "./tests/correctionArray_"+std::to_string(qs[d-i-1]));
    applyCorrection(cDVec[d-i], cC, pfVec[d-i-1], t);
    interpolateGrid(cDVec[d-i-1], cDVec[d-i], st);
    secondCorrectionStep(mc, st, fineSizeA, fineSizeB, hS, pfVec[d-i-1], cDVec[d-i-1]);
    writeToFile(cDVec[d-i-1], "./tests/cC_"+std::to_string(qs[d-i-1]));
    writeToFile(cDVec[d-i-1], "./tests/cC_"+std::to_string(qs[d-i-1]));
    // writeToFile(cD, "./tests/cC_"+std::to_string(qs[d-i]));
  }

  return 0;
}

