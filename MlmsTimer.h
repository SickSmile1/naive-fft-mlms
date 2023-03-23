/* Copyright 2023 <Ilia Fedotov @ Uni Freiburg> */
#ifndef MLMSTIMER_H_
#define MLMSTIMER_H_

#include <chrono> //NOLINT

struct MeasureTime{
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

void runTimerLoops();

#endif  // MLMSTIMER_H_
