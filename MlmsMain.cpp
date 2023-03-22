/* Copyright 2023 <Ilia Fedotov @ Uni Freiburg> */
/* Using lava/matplotlib-cpp for Matplotlib */

#include <omp.h>
#include <cstdio>
#include <cstring>
#include <sstream>
#include <iostream>
#include "./Mlms.h"

int main(int argc, char* argv[]) {
  std::stringstream strs1;

  if (argc == 1) {
    strs1 << "200" << std::endl;
    strs1 << "100" << std::endl;
    strs1 << "50" << std::endl;
    strs1 << "0" << std::endl;
    strs1 << "1" << std::endl;
    std::cout << "running test with parameters:" << std::endl;
    std::cout << "size:" << 200 << " size_p:" << 100
      << " grid:" << 50 << " v:" << 0 << " E:" << 1 << std::endl;

  } else if (argc > 1 && argc <5) {
    std::cout << "not enough parameters. 5 numbers needed!";
    std::cout <<"\n<double>Size\n<double>pressure patch size\n";
    std::cout << "<int>grid size\n<double>v\n<double>E" << std::endl;
  } else {
    for (int i = 1; i < argc; i++) {
      strs1 << argv[i] << std::endl;
    }
  }

  double size;
  strs1 >> size;  // default: 200

  std::size_t size_p;
  strs1 >> size_p;  // default: 100

  std::size_t grid;
  strs1 >> grid;  // default: 50

  double cell_size = size/grid;

  double v;
  strs1 >> v;  // default: 0
  double E;
  strs1 >> E;  // default: 1

  // return 0;

  /*std::vector<std::vector<double>> Ic;
  std::vector<std::vector<double>> Pa;

  std::vector<double> help;

  help.clear();
  Ic.clear();
  Pa.clear();*/
  matrix Ic({grid, grid});
  matrix Pa({grid, grid});

  for (std::size_t i = 0; i < Ic.shape[0]; i++) {
    for (std::size_t j = 0; j < Ic.shape[1]; j++) {
      Ic(i, j) = 0;
    }
  }
  double lower_b = (size/cell_size)/2 - (size_p/cell_size)/2;
  double upper_b = (size/cell_size)/2 + (size_p/cell_size)/2;
  printf("%.1f lower, %.1f uppwe \n", lower_b, upper_b);
  // initialize pressure array
  for (std::size_t i = 0; i < Pa.shape[0]; i++) {
    for (std::size_t j = 0; j < Pa.shape[1]; j++) {
      if (j >= lower_b && i >= lower_b && j < upper_b && i < upper_b) {
        Pa(i, j) = 1;
      } else {
        Pa(i, j) = 0;
      }
    }
  }

  // printarray(Ic, grid_x, grid_y);

  // printarray(Pa, grid_x, grid_y);
  // run calculation loops
  // #pragma omp parallel for
  for (std::size_t i = 0; i < grid; i++) {
    for (std::size_t j = 0; j < grid; j++) {
      Ic(i, j) = run_grid(Pa, i, j, grid, cell_size/2, cell_size/2, v, E,
                          cell_size);
    }
  }

  printarray(Ic);
}
