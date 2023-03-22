/* Copyright 2023 <Ilia Fedotov @ Uni Freiburg> */
#include <fstream>
#include <cmath>
#include <vector>
#include <array>
#include "./Mlms.h"

const double PI = 3.141592653589793238463;

// __________________________________________________________________
void printarray(const matrix &array) {
  for (int i = 0; i < array.shape[0]; i++) {
    for (int j = 0; j < array.shape[1]; j++) {
      printf("%.1f\t", array(i, j));
    }
    printf("\n");
  }
  return;
}


// __________________________________________________________________
double precalculation(double a, double b, double x, double y) {
  double first, second, third, fourth;

  first = ((y + b) + sqrt(pow((y + b), 2) + pow((x + a), 2))) /
          ((y - b) + sqrt(pow((y - b), 2) + pow((x + a), 2)));
  first = (x + a) * log(first);
  // printf("%.2f \n", first);

  second = ((x + a) + sqrt(pow((y + b), 2) + pow((x + a), 2))) /
          ((x - a) + sqrt(pow((y + b), 2) + pow((x - a), 2)));
  second = (y + b) * log(second);
  // printf("%.2f \n", second);

  third = ((y - b) + sqrt(pow((y - b), 2) + pow((x - a), 2))) /
          ((y + b) + sqrt(pow((y + b), 2) + pow((x - a), 2)));
  third = (x - a) * log(third);
  // printf("%.2f \n", third);

  fourth = ((x - a) + sqrt(pow((y - b), 2) + pow((x - a), 2))) /
            ((x + a) + sqrt(pow((y - b), 2) + pow((x + a), 2)));
  fourth = (y - b) * log(fourth);
  // printf("%.2f \n", fourth);
  // printf("%.2f , %.2f, %.2f, %.2f",first,second,third,fourth );
  // if (std::isnan(first)) first = 0;
  // if (std::isnan(second)) second = 0;
  // if (std::isnan(third)) third = 0;
  // if (std::isnan(fourth)) fourth = 0;
  double res = (first+second+third+fourth);
  return res;
}

// __________________________________________________________________
double run_grid(const matrix &pressure,
                double y, double x, int grid,
                double a, double b , double v,
                double E, double cell) {
  double res = 0;

  for (int i = 0; i < grid; i+=1) {
    for (int j = 0; j < grid; j+=1) {
      res += precalculation(a, b, (x-j)*cell,
                            (y-i)*cell) *
                            ((1-v)/(PI*E)) * pressure(i, j);
      // if (j == 10) printf("%2.f \n", res);
    }
  }
  return res;
}
