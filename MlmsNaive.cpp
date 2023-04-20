/* Copyright 2023 <Ilia Fedotov @ Uni Freiburg> */
#include <fstream>
#include <iostream>
#include <cmath>
#include <vector>
#include "./Mlms.h"

const double PI = 3.141592653589793238463;

// __________________________________________________________________
void printarray(const matrix &array) {
  for (int i = 0; i < array.shape[0]; i++) {
    for (int j = 0; j < array.shape[1]; j++) {
      printf("%.4f\t", array(i, j));
    }
    printf("\n");
  }
  return;
}

// __________________________________________________________________
void initializePressureArray(matrix &Pa, double lower_b, // NOLINT
                             double upper_b, double pressure) {
  initializeDisplacementArray(Pa);
  for (int i = lower_b; i < upper_b; i++) {
    for (int j = lower_b; j < upper_b; j++) {
      Pa(i, j) = pressure;
    }
  }
}

// __________________________________________________________________
void initializeDisplacementArray(matrix &Ic) { // NOLINT
  for (int i = 0; i < Ic.shape[0]; i++) {
    for (int j = 0; j < Ic.shape[1]; j++) {
      Ic(i, j) = 0;
    }
  }
}

// __________________________________________________________________
inline double calculate(double a, double b, double x, double y) {
  double first, second, third, fourth;
  // naive formula implementation for calculating a pressure patch
  first = ((y + b) + sqrt(pow((y + b), 2) + pow((x + a), 2))) /
          ((y - b) + sqrt(pow((y - b), 2) + pow((x + a), 2)));
  first = (x + a) * log(first);

  second = ((x + a) + sqrt(pow((y + b), 2) + pow((x + a), 2))) /
           ((x - a) + sqrt(pow((y + b), 2) + pow((x - a), 2)));
  second = (y + b) * log(second);

  third = ((y - b) + sqrt(pow((y - b), 2) + pow((x - a), 2))) /
          ((y + b) + sqrt(pow((y + b), 2) + pow((x - a), 2)));
  third = (x - a) * log(third);

  fourth = ((x - a) + sqrt(pow((y - b), 2) + pow((x - a), 2))) /
            ((x + a) + sqrt(pow((y - b), 2) + pow((x + a), 2)));
  fourth = (y - b) * log(fourth);
  double res = (first+second+third+fourth);
  return res;
}


// __________________________________________________________________
void calculation_loop(matrix &Ic, const matrix &Pa, // NOLINT
                double cell_size,
                double v, double E) {
  // outer loop over grid to call displacement calculation
  // for parallelisation uncomment, compiler option needed: -fopen
  #pragma omp parallel for simd
  for (int i = 0; i < Ic.shape[0]; i++) {
    for (int j = 0; j < Ic.shape[1]; j++) {
      Ic(i, j) = calc_displacement(Pa, Ic, i, j, cell_size/2,
                          cell_size/2, v, E, cell_size);
    }
  }
}

// __________________________________________________________________
double calc_displacement(const matrix &pressure,
              const matrix &Ic,
              double y, double x,
              double a, double b , double v,
              double E, double cell) {
  double res = 0;

  for (int i = 0; i < Ic.shape[0]; i+=1) {
    for (int j = 0; j < Ic.shape[1]; j+=1) {
      double xj = (x*cell)-(j*cell);
      double yi = (y*cell)-(i*cell);
      // std::cout << "xi: " << xj << " yi: " << yi << std::endl;
      res += calculate(a, b, xj, yi) *
                            ((1-v)/(PI*E)) * pressure(i, j);
    }
  }
  return res;
}



// __________________________________________________________________
void initializeMultiplicationArray(matrix &Ic) { // NOLINT
  for (int i = 0; i < Ic.shape[0]; i++) {
    for (int j = 0; j < Ic.shape[1]; j++) {
      Ic(i, j) = 1;
    }
  }
}

// __________________________________________________________________
bool boundaryCheck(const matrix &m, int i, int j) { // NOLINT
  bool res = (i < m.shape[0]) & (j < m.shape[1]) & (i >= 0) & (j >= 0);
  return res;
}

// __________________________________________________________________
bool boundaryCheck(matrix &m, int i, int j) { // NOLINT
  bool res = (i < m.shape[0]) & (j < m.shape[1]) & (i >= 0) & (j >= 0);
  return res;
}

// __________________________________________________________________
void myBreakpoint() {}

// __________________________________________________________________
void correctionSteps(matrix& cC, const matrix& st, int mc, int t, // NOLINT
    int fineSizeA, int fineSizeB, int halfSize) {
  for (int i = -mc; i <= mc; i++) {
    for (int j = -mc; j <= mc; j++) {
      bool iEven = (i%2 == 0);
      bool jEven = (j%2 == 0);
      double res1, res2, res3 = 0;
      int ts = t;
      double K = calculate(fineSizeA/2, fineSizeB/2, i*halfSize,
                  j*halfSize);

      for (int k = 1; k <= 2*t; k++) {
        res1 += st(k-1, 0) * calculate(fineSizeA/2, fineSizeB/2,
                              (i-2.*(k-ts)+1.)*halfSize,
                              j*halfSize);
      }
      for (int l = 1; l <= 2*t; l++) {
        res2 += st(l-1, 0) * calculate(fineSizeA/2, fineSizeB/2,
                              i*halfSize,
                              halfSize*(j-2.*(l-ts)+1.));
      }
      for (int k = 1; k <= 2*t; k++) {
        for (int l = 1; l <= 2*t; l++) {
          res3 += st(k-1, 0)*st(l-1, 0) * calculate(fineSizeA/2,
                              fineSizeB/2, (i-2.*(k-ts)+1.)*
                              halfSize, halfSize*(j-2.*(l-ts)+1.));
        }
      }
      cC(i+mc, j+mc) += ((!iEven)&&jEven)*(K-res1);
      cC(i+mc, j+mc) += (iEven&&(!jEven))*(K-res2);
      cC(i+mc, j+mc) += ((!iEven)&&(!jEven))*(K-res3);
      cC(i+mc, j+mc) += (iEven&&jEven)*0;
    }
  }
  // writeToFile(cC, "./tests/cC");
  // printarray(cC);
}

// __________________________________________________________________
void applyCorrection(matrix &cD, const matrix cC, const matrix Ip, // NOLINT
    int t) {
  for (int i = 0; i < cD.shape[0]; i++) {
    for (int j = 0; j < cD.shape[1]; j++) {
      cD(i, j) += correctionHelper(cC, Ip, t, i, j);
    }
  }
}

// __________________________________________________________________
double correctionHelper(const matrix& cC, const matrix& Ip, int t,
    int i, int j) {
  int mc = cC.shape[0]/2;
  double res = 0;
  for (int k = -mc; k <= mc; k++) {
    for (int l = -mc; l <= mc; l++) {
      int pi = 2*(i-t+1);
      int pj = 2*(j-t+1);
      if (boundaryCheck(Ip, pi, pj)) {
        res += cC(k+mc, l+mc) * Ip(pi, pj);
      }
    }
  }
  return res;
}

// __________________________________________________________________
void interpolateGrid(matrix &fD, const matrix cD, const matrix st) { // NOLINT
  int t = st.shape[0]/2;
  for (int i = 0; i < fD.shape[0]; i++) {
    for (int j = 0; j < fD.shape[1]; j++) {
      bool iEven = (i%2) == 0;
      bool jEven = (j%2) == 0;
      bool once = true;
      for (int k = 1; k <= 2*t; k++) {
        for (int l = 1; l <= 2*t; l++) {
          int m = (i + 1) / 2 + t - 1;
          int n = (j + 1) / 2 + t - 1;
          int pm = m+k-t-1;
          int pn = n+l-t-1;
          if ((!iEven)&jEven) {
            fD(i, j) += st(k-1, 0) * cD(pm, n);
          }
          if (iEven&(!jEven)) {
            fD(i, j) += (st(l-1, 0) * cD(m, pn)) * once;
          } else if ((!iEven)&(!jEven)) {
            fD(i, j) += st(k-1, 0) * st(l-1, 0) * cD(pm, pn);
          }
        }
        once = false;
      }
      // fD(i, j) += (iEven&jEven)*res;
      // fD(i, j) += (iEven&jEven)*res1;
      // fD(i, j) += (iEven&jEven)*res2;
      if (boundaryCheck(cD, (i+1)/2+t-1, (j+1)/2+t-1)) {
        fD(i, j) += (iEven&jEven)*cD((i+1)/2+t-1, (j+1)/2+t-1);
      }
    }
  }
}
/*
// __________________________________________________________________
void calc_coarse_pressure(const matrix &fP, const matrix &st, matrix &cP,
                          std::size_t ts) {
  int t = ts;
  for (std::size_t m = 0; m < cP.shape[0]; m++) {
    for (std::size_t n = 0; n < cP.shape[1]; n++) {
      int i = 2*(m-t+1);
      int j = 2*(n-t+1);
      // cP(m,n) = 0;

      double c_pressure;
      if (i < fP.shape[0] && j < fP.shape[1] && 
          i >= 0 && j >= 0) {
        c_pressure = fP(i, j);
      } else {
        c_pressure = 0;
        continue;
      }
      for (int k = 1; k <= 2*t; k++) {
        int pi = i+2*(k-t)-1;
        if (pi > 0 && pi < fP.shape[0] && j < fP.shape[1]) {
          c_pressure += st(k - 1, 0)*fP(pi, j);
        }
      }
      for (int l = 1; l <= 2*t; l++) {
        int pj = j+2*(l-t)-1;
        if (pj > 0 && pj < fP.shape[1] && i < fP.shape[0]) {
          c_pressure += st(l - 1, 0)*fP(i, pj);
        }
      }

      for (int k = 1; k <= 2*t; k++) {
        for (int l = 1; l <= 2* t; l++) {
          int si = i+2*(k-t)-1;
          int sj = j+2*(l-t)-1;
          if (si > 0 && sj > 0 && si < fP.shape[0] && sj < fP.shape[1]) {
            c_pressure += st(k - 1, 0) * st(l - 1, 0) * fP(si, sj);
          }
        }
      }

      cP(m, n) += c_pressure;
    }
  }
}

// __________________________________________________________________
void prepareCoarseSizes(std::vector<std::size_t> &gridLen1,
                        std::vector<std::size_t> &gridLen2,
                        std::vector<std::size_t> &coarseGrid,
                        const std::size_t shape,
                        const std::size_t ts) {
  int t = ts;
  double sizeN = shape*shape;
  double sqrt_sizeN = shape;

  std::size_t newGrid1 = shape/2 + 2*t - 1;
  std::size_t newGrid2 = shape/2 + 2*t - 1;

  coarseGrid.reserve(120);
  gridLen1.reserve(120);
  gridLen2.reserve(120);

  while (sizeN > sqrt_sizeN) {
    sizeN = newGrid1 * newGrid2;
    newGrid1 = newGrid1/2 + 2*t - 1;
    newGrid2 = newGrid2/2 + 2*t - 1;
    coarseGrid.push_back(sizeN);
    gridLen1.push_back(newGrid1);
    gridLen2.push_back(newGrid2);
    // std::cout << sizeN << std::endl;
  }
  coarseGrid.shrink_to_fit();
  gridLen1.shrink_to_fit();
  gridLen2.shrink_to_fit();
}

// __________________________________________________________________
void calcCorrMatrix(matrix &correctionCoefficients, const double mc, 
                        const matrix &st,
                        const std::size_t ts, const double fineSizeA,
                        const double fineSizeB) {
  int t = ts;
  for (int i = -mc; i <= mc; i++) {
    for (int j = -mc; j <= mc; j++) {
      bool iEven = (i % 2 == 0);
      bool jEven = (j % 2 == 0);

      double K;
      // odd i, even j
      double res = 0;
      for (int k = 1; k <= 2*t; k++) {
        int pi = i - 2 * (k - t) + 1;
        K = calculate(fineSizeA, fineSizeB, pi, j);
        res += st(k-1, 0) * K;
      }

      // even i, odd j
      double res1 = 0;
      for (int l = 1; l <= 2*t; l++) {
        int pj = j - 2 * (l - t) + 1;
        K = calculate(fineSizeA, fineSizeB, i, pj);
        res1 += st(l-1, 0) * K;
      }

      // odd i, j
      double res2 = 0;
      for (int l = 1; l <= 2*t; l++) {
        for (int k = 1; k <= 2*t; k++) {
          int pj = j - 2 * (l - t) + 1;
          int pi = i - 2 * (k - t) + 1;
          K = calculate(fineSizeA, fineSizeB, pi, pj);
          res2 += st(l-1, 0) * st(k-1, 0) * K;
        }
      }
      double result = 0;
      K = calculate(fineSizeA, fineSizeB, i, j);
      result += (jEven*!iEven) * (K - res);
      result += (!jEven*iEven) * (K - res1);
      result += (!jEven*!iEven) * (K - res2);
      if (jEven*iEven && i > 0 && j > 0) {
        correctionCoefficients(i, j) = 0;
      } else if (i > 0 && j > 0) {
        correctionCoefficients(i, j) = K - result*(jEven*iEven);
      }
    }
  }
}

void deflectionCorrection(matrix &kM, const double mc,
                          const matrix &cC, const matrix &cP,
                          const std::size_t ts) {
  // (13) (14)
  int t = ts;
  for (int i = 0; i < kM.shape[0]; i++) {
    for (int j = 0; j < kM.shape[1]; j++) {
      // summation over -mc <= k/l <= mc
      for (int k = -mc; k <= mc; k++) {
        for (int l = -mc; l <= mc; l++) {
          int pi = 2 * (i - t + 1) - k;
          int pj = 2 * (j - t + 1) - l;
          // std::cout << "startin correction\n";
          if (pi > 0 && pj > 0 && pi < cP.shape[0] && pj < cP.shape[1]&&
              (k+mc)>=0 && (k+mc)<cC.shape[0] && (l+mc)>=0 && (l+mc)<cC.shape[1]) {
            kM(i, j) += cC(k+mc, l+mc) * cP(pi, pj);
          }
        }
      }
    }
  }
}


void coarseToFineGrid(const matrix &cM, matrix &kM, const matrix &st,
                      const std::size_t ts) {
  // (15) and (16) in paper
  int t = ts;
  for (std::size_t i = 0; i < kM.shape[0]; i++) {
    for (std::size_t j = 0; j < kM.shape[1]; j++) {
      bool iEven = (i % 2 == 0);
      bool jEven = (j % 2 == 0);
      int m = (i + 1)/2 + t - 1;
      int n = (j + 1)/2 + t - 1;

      int res0 = 0;
      for (int k = 1; k <= 2*t; k++) {
        int pi = m + k - t - 1;
        if (pi >= 0 && pi< cM.shape[0] && n < cM.shape[1]) {
          res0 += st(k-1, 0) * cM(pi, n);
        }
      }
      int res1 = 0;
      for (int l = 1; l <= 2*t; l++) {
        int pj = n + l - t - 1;
        if (pj >= 0 && pj < cM.shape[1] && m < cM.shape[0]) {
          res1 += st(l-1, 0) * cM(m, pj);
        }
      }
      int res2 = 0;
      for (int k = 1; k <= 2*t; k++) {
        for (int l = 1; l <= 2*t; l++) {
          int pi = m + k - t - 1;
          int pj = n + l - t - 1;
          if (pi >= 0 && pj >= 0 && pi < cM.shape[0] && pj < cM.shape[1]) {
            res2 += st(k-1, 0) * st(l-1, 0) * cM(pi, pj);
          }
        }
      }
      double result = 0;
      result += (!iEven*jEven)*res0;
      result += (iEven*!jEven)*res1;
      result += (!iEven*jEven)*res2;
      result += (iEven*jEven)*kM(i, j);
      if (i > 0 && j > 0 && i < kM.shape[0] && j < kM.shape[1]) {
        kM(i, j) = result;
      }
    }
  }
}

void fineGridCorrection(matrix &kM, const matrix &st,
                  const double mc, const std::size_t ts,
                  const matrix fP) {
  // (16) and (17) in paper
  int t = ts;
  matrix cC1({2*mc, 2*mc});
  initializeDisplacementArray(cC1);
  matrix cC2({2*mc, 2*mc});
  matrix cC3({2*mc, 2*mc});
  matrix cC4({2*mc, 2*mc});

  for (int i = -mc; i < mc; i++) {
    for (int j = -mc; j < mc; j++) {
      bool iEven = (i % 2 == 0);
      bool jEven = (j % 2 == 0);

      int res0 = 0;
      for (int k = 1; k < 2*t; k++) {
        int pi = i + 2 * (k - t) - 1;
        if (pi >= 0 && pi < kM.shape[0] && j < kM.shape[1] && j>=0) {
          res0 += st(k-1, 0) * kM(pi, j);
        }
      }
      int res1 = 0;
      for (int l = 1; l < 2*t; l++) {
        int pj = j + 2 * (l - t) - 1;
        if (pj >= 0 && pj < kM.shape[1] && i < kM.shape[0] && i >=0) {
          res1 += st(l-1, 0) * kM(i, pj);
        }
      }
      int res2 = 0;
      for (int k = 1; k < 2*t; k++) {
        for (int l = 1; l < 2*t; l++) {
          int pi = i + 2 * (k - t) - 1;
          int pj = j + 2 * (l - t) - 1;
          if (pi >= 0 && pj >= 0 && pi < kM.shape[0] && pj < kM.shape[1]) {
            res2 += st(k-1, 0) * st(l-1, 0) * kM(pi, pj);
          }
        }
      }
      cC2(i+mc, j+mc) = kM(i+mc, j+mc) - ((!iEven*jEven) * res0);
      cC3(i+mc, j+mc) = kM(i+mc, j+mc) - ((iEven*!jEven) * res1);
      cC4(i+mc, j+mc) = kM(i+mc, j+mc) - ((!iEven*jEven) * res2);
    }
  }

  for (std::size_t i = 0; i < kM.shape[0]; i++) {
    for (std::size_t j = 0; j < kM.shape[1]; j++) {
      bool iEven = (i % 2 == 0);
      bool jEven = (j % 2 == 0);

      matrix &c = cC1;
      if (!iEven*jEven) c = cC2;
      else if (iEven*!jEven) c = cC3;
      else if (!iEven*jEven) c = cC4;
      double res = 0;
      for (int k = -mc; k < mc; k++) {
        for (int l = -mc; l < mc; l++) {
          int pi = i - k;
          int pj = j - l;
          if (pi > 0 && pj > 0 && pi < fP.shape[0] && pj < fP.shape[1]) {
            res += c(k+mc, l+mc)* fP(pi, pj);
          }
        }
      }
      kM(i, j) = kM(i, j) + res;
    }
  }
}
*/
