/* Copyright 2023 <Ilia Fedotov @ Uni Freiburg> */

// #include <cstdio>
// #include <cmath>
// #include <numbers>
// #include <fstream>
// #include <string>
// #include <iostream>
// #include <vector>
#include "./Mlms.h"

int main() {
  grid_a = 20;
  grid_b = 20;
  a = 5.;
  b = 5.;
  v = 1;
  E = 1;

  help.clear();
  Ic.clear();
  Pa.clear();
  for (size_t i = 0; i < grid_b; i++) {
    for (size_t j = 0; j < grid_a; j++) help.push_back(0.);
    Ic.push_back(help);
    help.clear();
  }

  // initialize pressure array
  for (size_t i = 0; i < grid_b; i++) {
    for (size_t j = 0; j < grid_a; j++) {
      if (j <= 10 && i <= 10) {
        help.push_back(1.);
      } else {
        help.push_back(0.);
      }
      Pa.push_back(help);
      help.clear();
    }
  }
  // run calculation loops
  for (size_t i = 0; i < grid_b; i++) {
    for (size_t j = 0; j < grid_a; j++) {
      // g_point = &Ic[i][j];
      // printf("%.2f \n",Ic[i][j]);
      // *g_point
      Ic[i][j] = run_grid(Ic, i, j, Pa[i][j], grid_a, grid_b, a, b, v, E);
      //  printf("%.2f \n",Ic[i][j]);
    }
  }

  printarray(Ic, grid_a, grid_b);
}
