/* Copyright 2023 <Ilia Fedotov @ Uni Freiburg> */

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include "MlmsTimer.h"
#include "Mlms.h"
#include "MlmsWriteToFile.h"

void runTimerLoops() {
  const double size = 1000;

  const std::size_t size_p = 500;
  std::size_t grid = 240;

  const double v = 0;
  const double E = 1;

  matrix result({10, 2});
  std::vector<matrix> matrix_results;
  int iteration;
  for (iteration = 0; iteration < 5; iteration++) {
    double cell_size = size/grid;

    matrix Ic({grid, grid});
    matrix Pa({grid, grid});

    double pressure = 1.;
    double lower_b = (size/cell_size)/2 - (size_p/cell_size)/2;
    double upper_b = (size/cell_size)/2 + (size_p/cell_size)/2;
    // std::cout << lower_b << "->lower_b " << upper_b << "->upper_b" << std::endl;

    initializePressureArray(Pa, lower_b, upper_b, pressure);
    initializeDisplacementArray(Ic);

    MeasureTime mt = MeasureTime();
    calculation_loop(Ic, Pa, cell_size, v, E);
    auto stopped =  mt.stopTime();

    std::string name = std::to_string(grid);
    name += "_grid";
    writeToFile(Ic, name);
    writeToFile(result, "timings.txt");

    result(iteration, 0) = grid;
    result(iteration, 1) = stopped;
    grid += 20;
    std::cout << "time:" << stopped << " grid:" << grid << std::endl;
  }
}

