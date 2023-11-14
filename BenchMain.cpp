#include <benchmark/benchmark.h>
#include "Boussinesq.h"
#include "BoussinesqFft.h"
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
void Naive_c(size_t grid) {
  int grids = grid;
  double size = 2;
  double size_p = 1;
  double cell_size = size/grid;

  matrix Ic({grids, grids});
  matrix Pa({grids, grids});

  double pressure = 1.;
  double lower_b = (size/cell_size)/2 - (size_p/cell_size)/2;
  double upper_b = (size/cell_size)/2 + (size_p/cell_size)/2;

  initializePressureArray(Pa, lower_b, upper_b, pressure);
  initializeDisplacementArray(Ic);
  calc_displacement(Pa, cell_size, cell_size, Ic);
}

static void Naive(benchmark::State &state) {
  for (auto _ : state) {
      Naive_c(state.range(0));
  }
}

BENCHMARK(Naive)->RangeMultiplier(2)->Range(8, 8<<4)->Unit(benchmark::kMillisecond);

void MlmsLoop1(size_t grid) {
  double size = 2;
  double size_p = .5;
  double pressure = 1.;
  int grids = grid+1;
  double fineSize = size / grids;
  matrix kM({grids, grids});
  matrix Ip({grids, grids});
  double lower_b = (grids)/2. - (size_p/fineSize)/2.;
  double upper_b = (grids)/2. + (size_p/fineSize)/2.;
  initializePressureArray(Ip, lower_b, upper_b, pressure);
  int t = 0.84*log(grids);
  // int t = 2;
  int mc = std::max(0.7*t*std::pow(grids,1./t)-1,t*1.);
  std::vector<matrix> pfVec, cDVec, cCVec;
  matrix st = initializeStylusArray(t);
  initializeStack(t, Ip, kM, pfVec, cDVec);
  double coarseSize = fineSize*pow(2, pfVec.size()-1);
  cCVec.reserve(3);
  createCorrectionArrays(cCVec, st, coarseSize, fineSize, mc);
  calcCoarsePressure(pfVec, st);
  int d = pfVec.size()-1;
  calc_displacement(pfVec[d], coarseSize, fineSize, cDVec[d]);
  for (std::size_t i = 0; i < pfVec.size()-1; i++) {
    double hS = fineSize*pow(2, d-i-1);
    int temp_mc = (mc*2)+1;
    matrix cC({temp_mc, temp_mc});
    correctionSteps(cC, st, mc, t, fineSize, hS);
    applyCorrection(cDVec[d-i], cC, pfVec[d-i-1], t, mc);
    interpolateGrid(cDVec[d-i-1], cDVec[d-i], st);
    secondCorrectionStep(pfVec[d-i-1], cDVec[d-i-1], cCVec, mc);
  }
}

static void Mlms1(benchmark::State &state) {
  for (auto _ : state) {
      MlmsLoop1(state.range(0));
  }
}

BENCHMARK(Mlms1)->RangeMultiplier(2)->Range(8, 8<<9)->Unit(benchmark::kMillisecond);

void FftLoop(size_t grids) {
  double Lx = 2., Ly = 2.;
  int Nx = grids+1, Ny = grids+1;
  double pSize = .5;
  double dx = (Lx/Nx);
  double dy = (Ly/Ny);

  int lb = Nx/2-(pSize/dx)/2;
  int ub = Ny/2+(pSize/dy)/2;

  matrix Gmn({(2*Nx)-1, (2*Ny)-1});
  cMatrix Gmn_tild({Gmn.rows(), Gmn.cols()/2+1});
  matrix p({Gmn.rows(), Gmn.cols()});
  cMatrix p_tild({Gmn.rows(), Gmn.cols()/2+1});
  matrix tempP({Nx, Ny});
  matrix Umn({Gmn.rows(), Gmn.cols()});
  cMatrix Umn_tild({Gmn.rows(), Gmn.cols()/2+1});
  matrix Umn_res({Nx, Ny});

  initializePressureArray(tempP, lb, ub, 1.);
  initializeDisplacementArray(p);

  copyPressureArray(p, tempP);

  calculateGmn(Gmn, dx, dy);

  transformGmnP(Gmn, Gmn_tild, p, p_tild);

  multiplyTransformed(Gmn_tild, Umn_tild, p_tild);

  transformToReal(Umn_tild, Umn);

  writeToResultArray(Umn, Umn_res);
  
  // return Umn_res;
}

static void FFT(benchmark::State &state) {
  for (auto _ : state) {
      FftLoop(state.range(0));
  }
}

BENCHMARK(FFT)->RangeMultiplier(2)->Range(8, 8<<9)->Unit(benchmark::kMillisecond);
// BENCHMARK(FFT)->RangeMultiplier(2)->DenseRange(300, 400,1)->Unit(benchmark::kMillisecond);

BENCHMARK_MAIN();
