#include <benchmark/benchmark.h>
// #include "Boussinesq.h"
// #include "BoussinesqFft.h"
// #include "BoussinesqMlms.h"
#include <iostream>
#include <cmath>
#include <fftw3.h>

#include "Boussinesq.hh"
#include "BoussinesqMatrix.hh"
#include "BoussinesqMlms.hh"

void MlmsLoop(size_t grid) {
  double Lx = 2, Ly = 2;
  std::size_t Nx = grid, Ny = grid;

  vec2d pixel = {Lx / Nx, Ly / Ny};

  matrix pressure({Nx, Ny}), displacement({Nx, Ny});

  // initialize pressure
  for (std::size_t i = 0; i < pressure.shape[0]; ++i) {
    for (std::size_t j = 0; j < pressure.shape[1]; ++j) {
      if (i >= Nx/4 and i < 3 * Nx / 4 and j >= Nx/4 and j < 3 * Ny / 4) {
        pressure(i, j) = 1.;
      } else {
        pressure(i, j) = 0.;
      }
    }
  }

  // compute deflection
  // naive(pressure, displacement, pixel, pixel);
  // std::cout << displacement;
  // return 0;

  int t = 6, levels = optimal_levels(pressure.shape);

  grid_stack pressure_stack(pressure, pixel, t, levels),
    displacement_stack(displacement, pixel, t, levels);
  auto mc = correction_size(t, pressure.shape);

  pressure_stack.coarsen();

  compute_coarse_displacements(pressure_stack, displacement_stack);

  auto corrections = coarse_corrections(pressure_stack, mc);
  auto fcorrections = fine_corrections(displacement_stack, mc);


  refine(pressure_stack, displacement_stack, corrections, fcorrections);
}

static void Mlms(benchmark::State &state) {
  for (auto _ : state) {
      MlmsLoop(state.range(0));
  }
  state.SetComplexityN(state.range(0));
}

BENCHMARK(Mlms)->RangeMultiplier(2)->Range(8, 8<<9)->Unit(benchmark::kMillisecond)->Complexity();

BENCHMARK_MAIN();
