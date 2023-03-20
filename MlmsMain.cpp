/* Copyright 2023 <Ilia Fedotov @ Uni Freiburg> */

#include <cstdio>
#include <omp.h>
#include "./Mlms.h"

int main() {
  std::vector<std::vector<double>> Ic;
  std::vector<std::vector<double>> Pa;

  std::vector<double> help;
  double size_x = 100;
  double size_y = 100;

  double size_p_x = 50;
  double size_p_y = 50;

  int grid_x = 100;
  int grid_y = 100;

  double cell_size_x = size_x/grid_x;
  double cell_size_y = size_y/grid_y;

  double a = size_p_x/2;
  double b = size_p_y/2;
  double v = 1;
  double E = 1;

  help.clear();
  Ic.clear();
  Pa.clear();
  for (int i = 0; i < grid_y; i++) {
    for (int j = 0; j < grid_x; j++) help.push_back(0.);
    Ic.push_back(help);
    help.clear();
  }

  // initialize pressure array
  for (int i = 0; i < grid_y; i++) {
    for (int j = 0; j < grid_x; j++) {
      if (j >= 25 && i >= 25 && j < 75 && i < 75) {
        help.push_back(1.0);
      } else {
        help.push_back(0.);
      }
    }
    Pa.push_back(help);
    help.clear();
  }

  // printarray(Ic, grid_x, grid_y);

  // printarray(Pa, grid_x, grid_y);
  // run calculation loops
  #pragma omp parallel for
  for (int i = 0; i < grid_y; i++) {
    for (int j = 0; j < grid_x; j++) {
      Ic[i][j] = run_grid(Ic, Pa, i, j, grid_x, grid_y, a, b, v, E,
                          cell_size_x, cell_size_y);
    }
  }

  printarray(Ic, grid_x, grid_y);
}
