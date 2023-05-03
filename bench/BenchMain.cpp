#include <benchmark/benchmark.h>
#include "Mlms.h"
#include <cmath>

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

static void bench_correction (benchmark::State &state) {
  const double size = 2;
  const double size_p = 1;
  const double pressure = 1.;

  // initial grid size for initialization
  int grid1 = 440;

  double fineSizeA = size / grid1;
  double fineSizeB = size / grid1;

  matrix kM({grid1, grid1});
  matrix Ip({grid1, grid1});

  double lower_b = (grid1)/2. - (size_p/fineSizeA)/2.;
  double upper_b = (grid1)/2. + (size_p/fineSizeB)/2.;

  initializePressureArray(Ip, lower_b, upper_b, pressure);
  // material moduli and v
  double v = 0.;
  double E = 1.;

  double beta = 0.84;
  double min_g = std::min(grid1, grid1);
  int t = beta*log(min_g);
  // int t = 4;
  int mc = 0.7* pow(min_g, 1./t)-1;
  if (mc < 2*t) {
    mc = 2*t;
  }
  matrix st({2*t, 1});

  initializeStylusArray(st, t);

  std::vector<matrix> pfVec;
  std::vector<matrix> cDVec;

  std::vector<int> qs;

  double qLevel = std::log2(grid1 * grid1)/2-1;
  int q = grid1 / 2 + 2*t - 1;
  qs.push_back(q);

  for (int length = 1; length < qLevel-1; length++) {
    q = q / 2 + 2*t - 1;
    qs.push_back(q);
  }

  double d = qs.size();
  pfVec.reserve(d);
  cDVec.reserve(d);
  pfVec.push_back(Ip);
  cDVec.push_back(kM);
  double coarseSize = fineSizeA*pow(2, d);
  std::vector<matrix> cCVec;
  cCVec.reserve(3);

  // create correction array for second correction -- duration appr.: 2412227 ns
  createCorrectionArrays(cCVec, st, coarseSize, fineSizeA,
                         fineSizeB, mc);
  
  // interpolate pressure to coarse grid -- duration appr.: 10580453 ns 
  calcCoarsePressure(qs, pfVec, cDVec, t, st);
  
  
  // calc displacement of coarsest grid with bousinesque kernel
  // appr.: 8528085 ns
    calc_displacement(pfVec[d], coarseSize, fineSizeA, cDVec[d]);
  
  
  // for (int i = 0; i < qs.size(); i++) {
    double hS = fineSizeA*pow(2, d-1) ; // i-1);
    int temp_mc = (mc*2)+1;
    matrix cC({temp_mc, temp_mc});
    // corr step appr. 578155 ns
    correctionSteps(cC, st, mc, t, fineSizeA, fineSizeB, hS);
    // return 0;
    // 178566 ns
    applyCorrection(cDVec[d], cC, pfVec[d-1], t);// i-1], t);
      
      // appr. 7694 ns
      interpolateGrid(cDVec[d-1], cDVec[d], st);// i-1], cDVec[d-i], st);
    
    for (auto _ : state) {
      secondCorrectionStep(mc, st, fineSizeA, fineSizeB, hS,
                         pfVec[d-1], cDVec[d-1], cCVec);//i-1], cDVec[d-i-1], cCVec);
    }
  // }
  

}

BENCHMARK(bench_correction)->Iterations(100);

BENCHMARK_MAIN();
