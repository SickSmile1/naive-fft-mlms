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
  const std::size_t size_p = 1;
  double pressure = 1.;

  // initial grid size for initialization
  int grid1 = 260;
  int grid2 = 260;

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
  double delta = 0;
  double g_old = 0;
  double min_g = std::min(grid1, grid2);
  int t = beta*log(min_g);
  double mc = 0.7* pow(min_g, 1./t)-1;
  if (mc < 2*t) {
    mc = 2*t;
  }
  double eps = pow(min_g, (-3./2.));

  double warnings = eps + size + size_p + pressure + delta + g_old + mc;
  warnings += warnings;

  matrix st({2*t, 1});
  initializeMultiplicationArray(st);
  // calculate transfer stylus (8)
  for (int i = 1; i <= 2*t; i++) {
    for (int j = 1; j <= 2*t; j++) {
      if (j != i) {
        // std::cout << st(i, 0) << " : ";
        double divider = (2.0*t-2.0*j)+1;
        double divisor = (2.0*i-2.0*j);
        st(i-1, 0) *= (divider/ divisor);
        // std::cout << st(i, 0) << std::endl;
      }
    }
  }

  int gridNum = grid1 / 2 + 2 * t - 1;
  matrix &pF = Ip;
  std::vector<int> gridSize;
  // gridSize.push_back(grid1);
  gridSize.push_back(gridNum);

  while (gridNum*gridNum >= grid1) {
    matrix pC({gridNum, gridNum});
    for (int m = 0; m < pC.shape[0]; m++) {
      for (int n = 0; n < pC.shape[1]; n++) {
        // 1 <= k/l <= 2*t
        int i = 2*(m-t+1);
        int j = 2*(n-t+1);
        double res = 0;
        bool notInFineGrid = boundaryCheck(pF, i, j);
        if (!notInFineGrid) {
          continue;
        } else {
          res += pF(i, j);
          for (int k = 1; k <= 2*t; k++) {
            int pi = i+2*(k-t)-1;
            res += boundaryCheck(pF, pi, j) ? st(k-1, 0)*pF(pi, j): 0;
          }
          for (int l = 1; l <= 2*t; l++) {
            int pj = j+2*(l-t)-1;
            res += boundaryCheck(pF, i, pj) ? st(l-1, 0)*pF(i, pj): 0;
          }

          for (int k = 1; k <= 2*t; k++) {
            for (int l = 1; l < 2*t; l++) {
              int pi = i+2*(k-t)-1;
              int pj = j+2*(l-t)-1;
              res += boundaryCheck(pF, pi, pj) ?
                      st(k-1, 0)*st(l-1, 0)*pF(pi, pj) : 0;
            }
          }
          pC(m, n) = res;
        }
      }
    }
    writeToFile(pC, "./tests/coarsened_"+std::to_string(gridNum));
    gridNum = gridNum / 2 + 2 * t - 1;
    if (gridNum*gridNum >= grid1) {
      gridSize.push_back(gridNum);
    }
    pF = pC;
  }

  // for (auto& e: gridSize) std::cout << e << std::endl;
  // std::cout << fineSizeA << " finesizeA\n"; 
  double d = gridSize.size()-1.;
  double halfSize = fineSizeA*pow(2, d-1);
  std::cout << "halfSize: "<<halfSize << std::endl;
  std::cout << "d:" << d << " gridSize[5]" << gridSize[5]<< std::endl;


  matrix coarseDisplacement({pF.shape[0], pF.shape[1]});
  for (int i = 0; i < pF.shape[0]; i++) {
    for (int j = 0; j < pF.shape[1]; j++) {
      for (int k = 0; k < pF.shape[0]; k++) {
        for (int l = 0; l < pF.shape[1]; l++) {
          double xj = (i*halfSize*2)-(k*halfSize*2);
          double yi = (j*halfSize*2)-(l*halfSize*2);
          coarseDisplacement(i, j) += calculate(fineSizeA, fineSizeB, xj, yi) *
                            ((1-v)/(PI*E)) * pF(i, j);
        }
      }
    }
  }
    writeToFile(coarseDisplacement, "./tests/ds");
  // printarray(coarseDisplacement);
  for (int i = 0; i < gridSize.size(); i++) {
    std::cout << d-i << " <- thats the d\n";
    double hS = fineSizeA*pow(2, d-i);
    std::size_t temp_mc = (mc*2)+1;
    matrix cC({temp_mc, temp_mc});
    correctionSteps(cC, st, mc, t, fineSizeA, fineSizeB, hS);
    applyCorrection(coarseDisplacement, cC, Ip, t);

    matrix nextDisplacement({gridSize[d-i] , gridSize[d-i]});
    interpolateGrid(nextDisplacement, coarseDisplacement, st);
    writeToFile(nextDisplacement, "./tests/cC_"+std::to_string(gridSize[d-i]));
  }

  return 0;
}

  /*
  // std::cout <<  << " is t" << std::endl;
  // printarray(st);

  // condition for coarse grid and pressure calculation
  std::vector<std::size_t> gridLen1;
  std::vector<std::size_t> gridLen2;
  std::vector<std::size_t> coarseGrid;


  // (9) in the paper
  prepareCoarseSizes(gridLen1, gridLen2, coarseGrid, kM.shape[0], t);
  myBreakpoint();
  matrix correctionCoefficients({mc*2, mc*2});
  calcCorrMatrix(correctionCoefficients, mc, st, t, fineSizeA, fineSizeB);
  writeToFile(Ip, "pressureGrid");
  // double q = coarseGrid.size();
  std::vector<matrix> coarseVector;

  for (int i = 0; i < coarseGrid.size(); i++) {
    // (10) in the paper
    // start iterative calculation_loop
    std::size_t gridLen = gridLen1[i];
    matrix coarsePressure({gridLen, gridLen});
    coarseVector.push_back(coarsePressure);
    double cellSize = size/ gridLen;
    myBreakpoint();

    if (i == 0){
      calc_coarse_pressure(Ip, st, coarseVector[0], t);
    } else {
      calc_coarse_pressure(coarseVector[i-1], st, coarseVector[i], t);
    }
    writeToFile(coarseVector[i], "./tests/coarse_grid"+std::to_string(gridLen));
  }

  // direct summation from Johnson, equal to (11) in paper
  std::size_t gridLen = gridLen1[coarseVector.size()-1];
  matrix coarsestMatrix({gridLen, gridLen});
  double cellSize = size/ gridLen;
  calculation_loop(coarsestMatrix, coarseVector[coarseVector.size()-1],
                    cellSize, v, E);

  for (int i = coarseVector.size()-1; i > 0; --i ) {
    std::size_t gridLen = gridLen1[i-1];
    matrix coarseMatrix({gridLen, gridLen});
    double cellSize = size/ gridLen;

    // writeToFile(coarseMatrix, "./tests/coarse_matrix"+std::to_string(gridLen));
    // (13) and (14) in paper
    deflectionCorrection(coarsestMatrix, mc, correctionCoefficients, Ip, t);
    writeToFile(coarsestMatrix, "./tests/coarse_corr_matrix"+
                std::to_string(gridLen));
    // (15) in paper
    coarseToFineGrid(coarsestMatrix, coarseMatrix, st, t);
    coarsestMatrix = coarseMatrix;
  }

  // (16) and (17) in paper
  fineGridCorrection(kM, st, mc, t, Ip);
  writeToFile(kM, "grid_faster");
  writeToFile(Ip, "grid_faster_pressure");
*/

