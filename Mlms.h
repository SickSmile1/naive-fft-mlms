/* Copyright 2023 <Ilia Fedotov @ Uni Freiburg> */
#ifndef MLMS_H_
#define MLMS_H_

// #include <gtest/gtest.h>
#include <vector>

void printarray(std::vector<std::vector<double>> array,
                 double grid_a, double grid_b);

double run_grid(std::vector<std::vector<double>> array,
                std::vector<std::vector<double>> pressure,
                double y, double x,
                int grid_a, int grid_b, double a,
                double b, double v, double E);

double precalculation(double a, double b, double x, double y);


#endif  // MLMS_H_
