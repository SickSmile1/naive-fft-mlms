/* Copyright 2023 <Ilia Fedotov @ Uni Freiburg> */
#include <vector>
#include <iostream>
#include <cmath>
#include "./MlmsWriteToFile.h"
#include "./Mlms.h"

int main() {
  const double PI = 3.141592653589793238463;

  // initial size and pressure values
  const double size = 500;
  const std::size_t size_p = 250;
  double pressure = 1.;

  // initial grid size for initialization
  std::size_t grid1 = 260;
  std::size_t grid2 = 260;

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
  std::size_t t = beta*log(min_g);
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
  for (std::size_t i = 1; i <= 2*t; i++) {
    for (std::size_t j = 1; j <= 2*t; j++) {
      if (j != i) {
        // std::cout << st(i, 0) << " : ";
        double divider = (2.0*t-2.0*j)+1;
        double divisor = (2.0*i-2.0*j);
        st(i-1, 0) *= (divider/ divisor);
        // std::cout << st(i, 0) << std::endl;
      }
    }
  }

  std::size_t gridNum = grid1 / 2 + 2 * t - 1;
  matrix &pF = Ip;
  std::vector<std::size_t> gridSize;
  gridSize.push_back(grid1);
  gridSize.push_back(gridNum);
  // printarray(pF);
  std::cout << gridNum << std::endl;
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
              res += boundaryCheck(pF, pi, pj) ? st(k-1, 0)*st(l-1, 0)*pF(pi, pj): 0;
            }
          }
          pC(m, n) = res;
        }
      }
    }
    writeToFile(pC, "./tests/coarsened_"+std::to_string(gridNum));
    gridNum = gridNum / 2 + 2 * t - 1;
    gridSize.push_back(gridNum);
    pF = pC;
  }

  double cell = size / pF.shape[0];
  double halfSize = cell/2;
  matrix coarseDisplacement({pF.shape[0], pF.shape[1]});
  for (int i = 0; i < pF.shape[0]; i++) {
    for (int j = 0; j < pF.shape[1]; j++) {
      for (int k = 0; k < pF.shape[0]; k++) {
        for (int l = 0; l < pF.shape[1]; l++) {
          double xj = (i*cell)-(k*cell);
          double yi = (j*cell)-(l*cell);
          coarseDisplacement(i, j) += calculate(halfSize, halfSize, xj, yi) *
                            ((1-v)/(PI*E)) * pF(i, j);
        }
      }
    }
  }
  printarray(coarseDisplacement);
  double d = gridSize.size()-1.;
  std::cout << halfSize << std::endl;
  std::size_t temp_mc = mc;
  matrix cC({temp_mc*2+1, temp_mc*2+1});
  for (int i = -mc; i <= mc; i++) {
    for (int j = -mc; j <= mc; j++) {
      bool iEven = (i%2==0);
      bool jEven = (j%2==0);
      double powi = (i*pow(2,d));
      double powj = (j*pow(2,d));
      double res1, res2, res3 = 0;
      double K = calculate(halfSize, halfSize, powj-(i*cell), powi-(i*cell))*
                  ((1-v)/(PI*E)) * pressure;

      for (int k = 1; k <= 2*t; k++) {
        res1 += st(k-1, 0) * calculate(halfSize, halfSize, powi-2*(k-t)+1,
                              powj)* ((1-v)/(PI*E)) * pressure;
      }
      for (int l = 1; l <= 2*t; l++) {
        res2 += st(l-1, 0) * calculate(halfSize, halfSize, powi,
                              powj-2*(l-t)+1) * ((1-v)/(PI*E)) * pressure;
      }
      for (int k = 1; k <= 2*t; k++) {
        for (int l = 1; l <= 2*t; l++) {
          res3 += st(k-1,0)*st(l-1, 0) * calculate(halfSize, halfSize, powi-2*(k-t)+1,
                              powj-2*(l-t)+1) * ((1-v)/(PI*E)) * pressure;
        }
      }
      cC(i+mc, j+mc) = (!iEven*jEven)*(K-res1);
      cC(i+mc, j+mc) = (iEven*!jEven)*(K-res2);
      cC(i+mc, j+mc) = (!iEven*!jEven)*(K-res3);
      cC(i+mc, j+mc) = (iEven*jEven)*0;
    }
    std::cout << (i*pow(2,d))-2*(-8-t)+1 << " thats my pow\n";
  }
  writeToFile(cC, "./tests/cC");
  // printarray(cC);



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
  return 0;
}

