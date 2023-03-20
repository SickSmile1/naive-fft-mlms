/* Copyright 2023 <Ilia Fedotov @ Uni Freiburg> */

#include <cstdio>
#include <omp.h>
#include "./Mlms.h"

int main() {
  std::vector<std::vector<double>> Ic;
  std::vector<std::vector<double>> Pa;

  std::vector<double> help;
  int grid_a = 20;
  int grid_b = 20;
  double a = 1.5;
  double b = 1.5;
  double v = 1;
  double E = 1;

  help.clear();
  Ic.clear();
  Pa.clear();
  for (int i = 0; i < grid_b; i++) {
    for (int j = 0; j < grid_a; j++) help.push_back(0.);
    Ic.push_back(help);
    help.clear();
  }

  // initialize pressure array
  for (int i = 0; i < grid_b; i++) {
    for (int j = 0; j < grid_a; j++) {
      if (j >= 5 && i >= 5 && j < 15 && i < 15) {
        help.push_back(1.0);
      } else {
        help.push_back(0.);
      }
    }
    Pa.push_back(help);
    help.clear();
  }

  // printarray(Ic, grid_a, grid_b);

  // printarray(Pa, grid_a, grid_b);
  // run calculation loops
  #pragma omp parallel for
  for (int i = 0; i < grid_b; i++) {
    for (int j = 0; j < grid_a; j++) {
      Ic[i][j] = run_grid(Ic, Pa, i, j, grid_a, grid_b, a, b, v, E);
    }
  }

  printarray(Ic, grid_a, grid_b);
}
