/* Copyright 2023 <Ilia Fedotov @ Uni Freiburg> */
#include <gtest/gtest.h>
#include <vector>

#ifndef MLMS_H_
#define MLMS_H_

std::vector<std::vector<double>> Ic;
std::vector<std::vector<double>> Pa;

std::vector<double> help;

double grid_a;
double grid_b;
double a;
double b;
double v;
double E;

double ln(double x);

void printarray(std::vector<std::vector<double>> array,
                 size_t grid_a, size_t grid_b);

double run_grid(std::vector<std::vector<double>> array, size_t y, size_t x,
                double pressure, double grid_a, double grid_b, double a,
                double b, double v, double E);

double precalculation(double a, double b, double x, double y, double v,
                      double E, double pressure);


#endif  // MLMS_H_
