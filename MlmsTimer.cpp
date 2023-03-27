/* Copyright 2023 <Ilia Fedotov @ Uni Freiburg> */

#include <iostream>
// #include <fstream>
// #include <sstream>
#include <vector>
#include <string>
#include "MlmsTimer.h"
#include "Mlms.h"
#include "MlmsWriteToFile.h"

void runTimerLoops() {
  const double size = 1000;

  const std::size_t size_p = 500;
  std::size_t grid = 60;

  const double v = 0;
  const double E = 1;
  double pressure = 1.;

  matrix result({2, 10});
  std::vector<matrix> matrix_results;
  int iteration;
  for (iteration = 0; iteration < 5; iteration++) {
    double cell_size = size/grid;

    matrix Ic({grid, grid});
    matrix Pa({10, 2});

    double lower_b = (size/cell_size)/2 - (size_p/cell_size)/2;
    double upper_b = (size/cell_size)/2 + (size_p/cell_size)/2;

    initializePressureArray(Pa, lower_b, upper_b, pressure);
    initializeDisplacementArray(Ic);

    MeasureTime mt = MeasureTime();
    calculation_loop(Ic, Pa, cell_size, v, E);
    auto stopped =  mt.stopTime();

    result(0, iteration) = grid;
    result(1, iteration) = stopped;

    std::string name = std::to_string(grid);
    name += "_grid";

    writeToFile(Ic, name);
    grid += 20;
  }
  writeToFile(result, "timings.txt");
}

