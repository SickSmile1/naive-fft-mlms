/* Copyright 2023 <Ilia Fedotov @ Uni Freiburg> */
#include <cmath>
#include <vector>
#include <fstream>
// #include "./Mlms.h"

const double PI = 3.141592653589793238463;

// __________________________________________________________________
/*void test_precalculation(double a, double b, double v,
                          double E, double pressure) {
  // future test for the precalculate function
  for (double i = 0.; i <= 10; i++) {
    for (double j = 0.; j <= 10; j++) {
      double res = precalculation(a, b, i, j, v, E, pressure);
      printf("res = %.8f i= %f, j=%f \n", res, i, j);
    }
  }
}*/

// __________________________________________________________________
void printarray(std::vector<std::vector<double>> array,
                 size_t grid_a, size_t grid_b) {
  for (size_t i = 0; i < grid_b; i++) {
    for (size_t j = 0; j < grid_a; j++) {
      printf("%.2f\t", array[i][j]);
    }
    printf("\n");
  }
  return;
}


// __________________________________________________________________
double precalculation(double a, double b, double x, double y,
                       double v, double E, double pressure) {
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
            ((x + a) + sqrt(pow((y - b), 2) + pow((x - a), 2)));
  fourth = (y - b) * log(fourth);
  // printf("%.2f \n", fourth);
  //  printf("%.2f , %.2f, %.2f, %.2f",first,second,third,fourth );
  double res = (first+second+third+fourth)*(1 /  // (1-pow(v, 2) ) /
                (PI*E)) * pressure;
  if (std::isnan(res)) res = 0;
  return res;
}

// __________________________________________________________________
double run_grid(std::vector<std::vector<double>> array, size_t y,
                size_t x, double pressure, double grid_a,
                double grid_b, double a, double b , double v,
                double E) {
  double res = 0;
  for (size_t i = 0; i < grid_b; i++) {
    for (size_t j = 0; j < grid_a; j++) {
      res += precalculation(a, b, x-j, y-i, v, E, pressure);
      // printf("%2.f", res);
    }
  }
  return res;
}
