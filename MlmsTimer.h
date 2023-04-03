/* Copyright 2023 <Ilia Fedotov @ Uni Freiburg> */
#ifndef MLMSTIMER_H_
#define MLMSTIMER_H_

#include <chrono> //NOLINT

struct MeasureTime{
  // time measurement struct for calculations
  // accurate until ,000x order
  std::chrono::_V2::steady_clock::time_point start;
  std::chrono::_V2::steady_clock::time_point stop;
  double duration;
  MeasureTime() {
    start = std::chrono::steady_clock::now();
  }

  double stopTime() {
    stop = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed_seconds = stop-start;
    return elapsed_seconds.count();
  }
};

// function calling the calculation with different grid sizes
// for the purpose of time measurement, saving the resulting matrix
// as 
//  xx_grid.txt -> xx being the grid size
//  timing.txt  -> gridsize \t time
void runTimerLoops();

#endif  // MLMSTIMER_H_
