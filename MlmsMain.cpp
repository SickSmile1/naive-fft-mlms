/* Copyright 2023 <Ilia Fedotov @ Uni Freiburg> */
/* Using lava/matplotlib-cpp for Matplotlib */

#include <cstdio>
// #include "matplotlibcpp.h"
#include <omp.h>
#include "./Mlms.h"

int main() {
  std::vector<std::vector<double>> Ic;
  std::vector<std::vector<double>> Pa;

  std::vector<double> help;
  double size = 200;

  double size_p = 100;

  int grid = 20;

  double cell_size = size/grid;

  double v = 0;
  double E = 1;

  help.clear();
  Ic.clear();
  Pa.clear();
  for (int i = 0; i < grid; i++) {
    for (int j = 0; j < grid; j++) help.push_back(0.);
    Ic.push_back(help);
    help.clear();
  }
  double lower_b = ((size/grid)*(size_p/cell_size));
  double upper_b = ((size/grid)+(size_p/cell_size));
  printf("%.1f lower, %.1f uppwe \n", lower_b, upper_b);
  // initialize pressure array
  for (int i = 0; i < grid; i++) {
    for (int j = 0; j < grid; j++) {
      if (j >= lower_b && i >= lower_b && j < upper_b && i < upper_b) {
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
  for (int i = 0; i < grid; i++) {
    for (int j = 0; j < grid; j++) {
      Ic[i][j] = run_grid(Pa, i, j, grid, cell_size/2, cell_size/2, v, E,
                          cell_size);
    }
  }

  printarray(Ic);
}
