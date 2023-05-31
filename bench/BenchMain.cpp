#include <benchmark/benchmark.h>
#include "Boussinesq.h"
#include "BoussinesqMlms.h"
#include <iostream>
#include <cmath>
#include <fftw3.h>

//#include "Boussinesq1.hh"
//#include "BoussinesqMatrix1.hh"
//#include "BoussinesqMlms1.hh"

//#include <iostream>
//#include <iomanip>

/*static void bench_alloc_and_null (benchmark::State &state) {
  std::size_t grid = 100;
  for (auto _ : state) {
  // naive formula implementation for calculating a pressure patch
    matrix s({grid,grid});
    initializeDisplacementArray(s);
  }
}

static void bench_alloc_1dvec (benchmark::State &state) {
  std::size_t grid = 100;
  for (auto _ : state) {
  // naive formula implementation for calculating a pressure patch
    matrix s({grid,grid});
  }
}

static void bench_init_pressure_array (benchmark::State &state) {
  std::size_t grid = 100;
  for (auto _ : state) {
  // naive formula implementation for calculating a pressure patch
    matrix s({grid,grid});
    initializePressureArray(s, 25, 75, 1);
  }
}

static void bench_rand_calculation_loop (benchmark::State &state) {
  const std::size_t grid = 40;
  const std::size_t size_p = 30;
  const std::size_t size = grid * 5;
  matrix s({grid,grid});
  matrix p({grid,grid});

  double cell_size = 5;

  double lower_b = (size/cell_size)/2 - (size_p/cell_size)/2;
  double upper_b = (size/cell_size)/2 + (size_p/cell_size)/2;

  initializeRandomPressureArray(p);

  initializeDisplacementArray(s);
  for (auto _ : state) {
  // naive formula implementation for calculating a pressure patch
    calculation_loop(s, p, cell_size, 0, 1);
    
  }
}

static void bench_calculation_loop (benchmark::State &state) {
  const std::size_t grid = 40;
  const std::size_t size_p = 30;
  const std::size_t size = grid * 5;
  matrix s({grid,grid});
  matrix p({grid,grid});

  double cell_size = 5;

  double lower_b = (size/cell_size)/2 - (size_p/cell_size)/2;
  double upper_b = (size/cell_size)/2 + (size_p/cell_size)/2;

  initializePressureArray(p, lower_b, upper_b, 1);

  initializeDisplacementArray(s);
  for (auto _ : state) {
  // naive formula implementation for calculating a pressure patch
    calculation_loop(s, p, cell_size, 0, 1);
  }
}

BENCHMARK(bench_alloc_1dvec);
BENCHMARK(bench_alloc_and_null);
BENCHMARK(bench_init_pressure_array);
BENCHMARK(bench_rand_calculation_loop)->Iterations(100);
BENCHMARK(bench_calculation_loop)->Iterations(100);
BENCHMARK(bench_calculation_loop)->Iterations(100);
BENCHMARK(bench_calculation_loop)->Iterations(100);
BENCHMARK(bench_calculation_loop)->Iterations(100);
*/

void MlmsLoop(size_t grid) {
  double size = 2;
  double size_p = 1;
  double pressure = 1.;
  int grids = grid;
  double fineSize = size / grids;
  matrix kM({grids, grids});
  matrix Ip({grids, grids});
  double lower_b = (grids)/2. - (size_p/fineSize)/2.;
  double upper_b = (grids)/2. + (size_p/fineSize)/2.;
  initializePressureArray(Ip, lower_b, upper_b, pressure);
  int t = 6;
  int mc = 2*t;
  std::vector<matrix> pfVec, cDVec, cCVec;
  matrix st = initializeStylusArray(t);
  initializeStack(st, t, Ip, kM, pfVec, cDVec);
  double coarseSize = fineSize*pow(2, pfVec.size()-1);
  cCVec.reserve(3);
  createCorrectionArrays(cCVec, st, coarseSize, fineSize);
  calcCoarsePressure(pfVec, st);
  int d = pfVec.size()-1;
  calc_displacement(pfVec[d], coarseSize, fineSize, cDVec[d]);
  for (int i = 0; i < pfVec.size()-1; i++) {
    double hS = fineSize*pow(2, d-i-1);
    int temp_mc = (mc*2)+1;
    matrix cC({temp_mc, temp_mc});
    correctionSteps(cC, st, mc, t, fineSize, hS);
    applyCorrection(cDVec[d-i], cC, pfVec[d-i-1], t);
    interpolateGrid(cDVec[d-i-1], cDVec[d-i], st);
    secondCorrectionStep(st, hS,
                          pfVec[d-i-1], cDVec[d-i-1], cCVec);
  }
}

static void Mlms(benchmark::State &state) {
  for (auto _ : state) {
      MlmsLoop(state.range(0));
  }

}

BENCHMARK(Mlms)->Range(8, 8<<10);

/*static void bench_Mlms_displacement(benchmark::State &state) {
  for (auto _:state) {
    calc_displacement(pfVec[d], coarseSize, fineSizeA, cDVec[d]);
  }
}

static void bench_Mlms_corrArr(benchmark::State &state) {
  for (auto _:state) {
    createCorrectionArrays(cCVec, st, coarseSize, fineSizeA,
                         fineSizeA, mc);
  }
}

static void bench_Mlms_corrStep(benchmark::State &state) {
  for (auto _:state) {
    correctionSteps(cC, st, mc, t, fineSizeA, fineSizeA, hS);
  }
}

static void bench_Mlms_corrStep1(benchmark::State &state) {
  for (auto _:state) {
    applyCorrection(cDVec[0], cC, pfVec[1], t);
  }
}

static void bench_Mlms_interpolate(benchmark::State &state) {
  for (auto _:state) {
    interpolateGrid(cDVec[1], cDVec[0], st);
  }
}
static void bench_Mlms_2ndCor(benchmark::State &state) {
  for (auto _:state) {
    secondCorrectionStep(mc, st, fineSizeA, fineSizeA, hS,
                       pfVec[0], cDVec[0], cCVec);
  }
}*/

// my implementation
/*BENCHMARK(bench_Mlms_displacement)->Iterations(25);
BENCHMARK(bench_Mlms_corrArr);
BENCHMARK(bench_Mlms_corrStep);
BENCHMARK(bench_Mlms_corrStep1);
BENCHMARK(bench_Mlms_interpolate);
BENCHMARK(bench_Mlms_2ndCor);*/

// Lukas implementation
/*
double Lx, Ly;
std::size_t Nx, Ny;

vec2d pixel;

matrix1 pressures({1, 1}), displacement({1, 1});
int t_l, levels;

grid_stack pressure_stack(pressures, {2 / 64, 2 / 64}, 4, 5),
  displacement_stack(displacement, {2 / 64, 2 / 64}, 4, 5);

std::vector<grid_stack> stack;

std::vector<matrix1> corrections;

std::vector<std::vector<matrix1>> fcorrections;
*/

/*
static void bench_Mlms_lukas(benchmark::State &state) {
  double Lx = 2, Ly = 2;
  std::size_t Nx = 256, Ny = 256;
  vec2d pixel = {Lx / Nx, Ly / Ny};
  matrix pressures({Nx, Ny}), displacement({Nx, Ny});
  // initialize pressure
  for (std::size_t i = 0; i < pressures.shape[0]; ++i) {
    for (std::size_t j = 0; j < pressures.shape[1]; ++j) {
      if (i >= Nx/4 and i < 3 * Nx / 4 and j >= Nx/4 and j < 3 * Ny / 4) {
        pressures(i, j) = 1.;
      } else {
        pressures(i, j) = 0.;
      }
    }
  }
  int t_l = 5, levels = optimal_levels(pressures.shape);
  
  grid_stack pressure_stack(pressures, pixel, t_l, levels),
    displacement_stack(displacement, pixel, t_l, levels);

  auto mc = correction_size(t, pressures.shape);
  pressure_stack.coarsen();
  
  for (auto _: state) {
    compute_coarse_displacements(pressure_stack, displacement_stack);
  }

  // stack.push_back(pressure_stack);
  // stack.push_back(displacement_stack);
  // auto corrections = coarse_corrections(pressure_stack, mc);
  // auto fcorrections = fine_corrections(displacement_stack, mc);

  // refine(pressure_stack, displacement_stack, corrections, fcorrections);
}*/


/*
static void bench_Mlms_displacement_l(benchmark::State &state) {
  for (auto _: state) {
    compute_coarse_displacements(stack[0], stack[1]);
  }
}

static void bench_Mlms_coarseCorr_l(benchmark::State &state) {
  for (auto _: state) {
    corrections = coarse_corrections(stack[0], mc);
  }
}

static void bench_Mlms_fineCor_l(benchmark::State &state) {
  for (auto _: state) {
   fcorrections = fine_corrections(stack[0], mc);
  }
}

static void bench_Mlms_refine_l(benchmark::State &state) {
  for (auto _: state) {
    refine(stack[0], stack[1], corrections, fcorrections);
  }
}

BENCHMARK(bench_Mlms_lukas);
BENCHMARK(bench_Mlms_displacement_l);
BENCHMARK(bench_Mlms_coarseCorr_l);
BENCHMARK(bench_Mlms_fineCor_l);
BENCHMARK(bench_Mlms_refine_l);*/

// BENCHMARK(bench_Mlms_lukas);

BENCHMARK_MAIN();
