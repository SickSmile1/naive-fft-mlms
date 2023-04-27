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
  const double size = 2;

  const double size_p = 1;
  int grid = 64;

  const double v = 0;
  const double E = 1;
  double pressure = 1.;

  matrix result({20, 2});
  int iteration;
  for (iteration = 0; iteration < 1; iteration++) {
    double cell_size = size/grid;

    matrix Ic({grid, grid});
    matrix Pa({grid, grid});

    double lower_b = (size/cell_size)/2 - (size_p/cell_size)/2;
    double upper_b = (size/cell_size)/2 + (size_p/cell_size)/2;

    initializePressureArray(Pa, lower_b, upper_b, pressure);
    initializeDisplacementArray(Ic);

    MeasureTime mt = MeasureTime();
    calculation_loop(Ic, Pa, cell_size, v, E);
    auto stopped =  mt.stopTime();

    result(iteration, 0) = grid;
    result(iteration, 1) = stopped;

    std::string name = std::to_string(grid);
    name += "_grid";

    writeToFile(Ic, name);
    grid += 20;
  }
  writeToFile(result, "timings.txt");
}

