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

  // initial grid size forinitialization
  std::size_t grid1 = 200;
  std::size_t grid2 = 200;

  double fineSizeA = size / grid1;
  double fineSizeB = size / grid2;

  matrix kM({grid1, grid2});
  matrix Ip({grid1, grid2});

  double lower_b = (grid1)/2 - (size_p/fineSizeA)/2;
  double upper_b = (grid2)/2 + (size_p/fineSizeB)/2;

  initializePressureArray(Ip, lower_b, upper_b, pressure);

  // material moduli and v
  double v = 0.;
  double E = 1.;

  double beta = 0.84;
  double delta = 0;
  double g_old = 0;
  double min_g = std::min(grid1, grid2);
  std::size_t t = beta*log(min_g);
  double mc = 0.7* pow(min_g, 1/t)-1;
  double eps = pow(min_g, (-3/2));

  double warnings = eps + size + size_p + pressure + delta + g_old + mc;
  warnings += warnings;

  matrix st({2*t, 1});
  initializeMultiplicationArray(st);
  // calculate transfer stylus (8)
  for (std::size_t i = 1; i <= 2*t; i++) {
    for (std::size_t j = 1; j <= 2*t; j++) {
      if (j != i) {
        // std::cout << st(i, 0) << " : ";
        double divider = (2*t-2*j)+1;
        double divisor = (2*i-2*j);
        st(i-1, 0) *= (divider/ divisor);
        // std::cout << st(i, 0) << std::endl;
      }
    }
  }

  // std::cout << t << " is t" << std::endl;
  // printarray(st);

  // condition for coarse grid and pressure calculation
  std::vector<std::size_t> gridLen1;
  std::vector<std::size_t> gridLen2;
  std::vector<std::size_t> coarseGrid;

  prepareCoarseSizes(gridLen1, gridLen2, coarseGrid, kM.shape[0], t);
  matrix correctionCoefficients({mc*2, mc*2});
  transferCoarseGrid(correctionCoefficients, mc, st, t, fineSizeA, fineSizeB);

  // double q = coarseGrid.size();

  for (auto qM = coarseGrid.rbegin(); qM != coarseGrid.rend(); ++qM) {
    // std::cout << *qM << std::endl;
    // start iterative calculation_loop
    std::size_t gridLen = sqrt(*qM);
    matrix coarseMatrix({gridLen, gridLen});
    matrix coarsePressure({gridLen, gridLen});
    double cellSize = size/ gridLen;

    std::cout << "startin press calc\n";
    calc_coarse_pressure(Ip, st, coarsePressure, t);
    writeToFile(coarsePressure, "./tests/coarse_grid"+std::to_string(gridLen));
    std::cout << "startin loop calcs\n";
    calculation_loop(coarseMatrix, coarsePressure, cellSize, v, E);
    writeToFile(coarseMatrix, "./tests/coarse_matrix"+std::to_string(gridLen));
    std::cout << "startin defl correction\n";
    deflectionCorrection(coarseMatrix, mc, correctionCoefficients, Ip, t);
    writeToFile(coarseMatrix, "./tests/coarse_corr_matrix"+
                std::to_string(gridLen));
    coarseToFineGrid(coarseMatrix, kM, st, t);
  }

  fineGridCorrection(kM, st, mc, t, Ip);
  writeToFile(kM, "grid_faster");
  writeToFile(Ip, "grid_faster2");

  return 0;
}

