
/* Copyright 2023 <Ilia Fedotov @ Uni Freiburg> */
#include <cstdio>
#include <cmath>
#include <numbers>
#include <fstream>
#include <string>
#include <iostream>

const size_t grid_a = 20;
const size_t grid_b = 20;
double Ic[grid_a][grid_b];
double Pa[grid_a][grid_b];
double a = 5.;
double b = 5.;
double v = 1;
double E = 1;
const double PI = 3.141592653589793238463;
void printarray( double( &array)[grid_a][grid_b], size_t grid_a, size_t grid_b);  //NOLINT
double run_grid(double(&array)[grid_a][grid_b], size_t y, size_t x, double pressure);  //NOLINT
double precalculation(double a, double b, double x, double y, double v, double E, double pressure);  //NOLINT

int main() {
  /*std::ifstream file("filetext");
  std::string str;
  while(std::getline(file,str))
  {
    std::cout << str << " -> is the read line \n";
    return 0;
  }*/
  double * g_point;
  // initialize testarray, 20x20, 10x10 pressure
  for (size_t i = 0; i < grid_b; i++) {
    for (size_t j = 0; j < grid_a; j++) {
        g_point = &Ic[i][j];
        *g_point = 0;
    }
  }

  // initialize pressure array
  for (size_t i = 0; i < grid_b; i++) {
    for (size_t j = 0; j < grid_a; j++) {
      if (j <= 10 && i <= 10) {
        g_point = &Pa[i][j];
        *g_point = 1;
      } else {
        g_point = &Pa[i][j];
        *g_point = 0;
      }
    }
  }
  // run calculation loops
  for (size_t i = 0; i < grid_b; i++) {
    for (size_t j = 0; j < grid_a; j++) {
      g_point = &Ic[i][j];
      // printf("%.2f \n",Ic[i][j]);
      *g_point = run_grid(Ic, i, j, Pa[i][j]);
      //  printf("%.2f \n",Ic[i][j]);
    }
  }

  printarray(Ic, grid_a, grid_b);
}

double run_grid(double(&array)[grid_a][grid_b], size_t y, size_t x,
                double pressure) {
  double res = 0;
  for (size_t i = 0; i < grid_b; i++) {
    for (size_t j = 0; j < grid_a; j++) {
      res += precalculation(a, b, x-j, y-i, v, E, pressure);
      // printf("%2.f", res);
    }
  }
  return res;
}

void test_precalculation(double a, double b, double v,
                          double E, double pressure) {
  // future test for the precalculate function
  for (double i = 0.; i <= 10; i++) {
    for (double j = 0.; j <= 10; j++) {
      double res = precalculation(a, b, i, j, v, E, pressure);
      printf("res = %.8f i= %f, j=%f \n", res, i, j);
    }
  }
}

void printarray(double(&array)[grid_a][grid_b],
                 size_t grid_a, size_t grid_b) {
  size_t i;
  size_t j;
  for (i = 0; i < grid_b; i++) {
    for (j = 0; j < grid_a; j++) {
      printf("%.2f\t", array[i][j]);
    }
    printf("\n");
  }
  return;
}

double precalculation(double a, double b, double x, double y,
                       double v, double E, double pressure) {
  double first, second, third, fourth;
  if ((x+a) == 0) first = 0;  //NOLINT
  else {
    first = (x + a) * log(
            ((y + b) + pow((pow((y + b), 2) + pow((x + a), 2)), 0.5)) /
            ((y - b) + pow((pow((y - b), 2) + pow((x + a), 2)), 0.5)));
  }

  if ((y+b) == 0) second = 0;  //NOLINT
  else {
    second = (y + b) * log(
            ((x + a) + pow((pow((y + b), 2) + pow((x + a), 2)), 0.5)) /
            ((x - a) + pow((pow((y + b), 2) + pow((x - a), 2)), 0.5)));
  }

  if ((x-a) == 0) third = 0;  //NOLINT
  else {
    third = (x - a) * log(
            ((y - b) + pow((pow((y - b), 2) + pow((x - a), 2)), 0.5) ) /
            ((y + b) + pow((pow((y + b), 2) + pow((x - a), 2)), 0.5)));
  }

  if ((y-b) == 0) fourth = 0;  //NOLINT
  else {
    fourth = (y - b) * log(
              ((x - a) + pow((pow((y - b), 2) + pow((x - a), 2)), 0.5)) /
              ((x + a) + pow((pow((y - b), 2) + pow((x - a), 2)), 0.5)));
  }

  double res = (first+second+third+fourth)*(1 /  // (1-pow(v, 2) ) /
                (PI*E)) * pressure;
  return res;
}
