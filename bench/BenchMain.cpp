#include <benchmark/benchmark.h>
#include "../Mlms.h"
#include <cmath>

static void bench_alloc_and_null (benchmark::State &state) {
  std::size_t grid = 100;
  while (state.KeepRunning()) {
  // naive formula implementation for calculating a pressure patch
    matrix s({grid,grid});
    initializeDisplacementArray(s);
  }
}

static void bench_alloc_1dvec (benchmark::State &state) {
  std::size_t grid = 100;
  while (state.KeepRunning()) {
  // naive formula implementation for calculating a pressure patch
    matrix s({grid,grid});
  }
}

static void bench_init_pressure_array (benchmark::State &state) {
  std::size_t grid = 100;
  while (state.KeepRunning()) {
  // naive formula implementation for calculating a pressure patch
    matrix s({grid,grid});
    initializePressureArray(s, 25, 75, 1);
  }
}

static void bench_calculation_loop (benchmark::State &state) {
  std::size_t grid = 10;
  matrix s({grid,grid});
  matrix p({grid,grid});
  initializePressureArray(p, 2, 7, 1);
  initializeDisplacementArray(s);
  while (state.KeepRunning()) {
  // naive formula implementation for calculating a pressure patch
    calculation_loop(s, p, 5, 0, 1);
    
  }
}

BENCHMARK(bench_alloc_and_null);
BENCHMARK(bench_alloc_1dvec);
BENCHMARK(bench_init_pressure_array);
BENCHMARK(bench_calculation_loop);

BENCHMARK_MAIN();
