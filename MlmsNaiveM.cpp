/* Copyright 2023 <Ilia Fedotov @ Uni Freiburg> */
/* Using lava/matplotlib-cpp for Matplotlib */

#include <omp.h>
#include <cstdio>
#include <cstring>
#include <sstream>
#include <iostream>
#include "./Mlms.h"
#include "./MlmsTimer.h"

void thisMain(int argc, char* argv[]) {
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

 /*  double size;
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

  matrix Ic({grid, grid});
  matrix Pa({grid, grid});

  double pressure = 1;
  double lower_b = (size/cell_size)/2 - (size_p/cell_size)/2;
  double upper_b = (size/cell_size)/2 + (size_p/cell_size)/2;

  initializePressureArray(Pa, lower_b, upper_b, pressure);
  initializeDisplacementArray(Ic);

  // run calculation loops
  calculation_loop(Ic, Pa, cell_size, grid, v, E);

  printarray(Ic); */
  runTimerLoops();
}
