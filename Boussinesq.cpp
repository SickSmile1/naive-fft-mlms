/* Copyright 2023 <Ilia Fedotov @ Uni Freiburg> */
#include <fftw3.h>
#include <fstream>
#include <iostream>
#include <cmath>
#include <vector>
#include <complex>
#include "./Boussinesq.h"

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
double calcBoussinesq(double a, double b, double x, double y) {
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
  double res = (first+second+third+fourth)/M_PI;
  return res;
}

// __________________________________________________________________
double calcBoussinesq(int i, int j, double dxc, double dyc,
    double dxf, double dyf) {
  return calcBoussinesq(dxf/2, dyf/2, i*dxc, j*dyc);
}
// __________________________________________________________________
void naiveCalculation(matrix &Ic, const matrix &Pa, // NOLINT
                double cell_size) {
  // outer loop over grid to call displacement calculation
  // for parallelisation uncomment, compiler option needed: -fopen
  #pragma omp parallel for simd
  for (int i = 0; i < Ic.shape[0]; i++) {
    for (int j = 0; j < Ic.shape[1]; j++) {
      Ic(i, j) = calc_displacement(Pa, Ic, i, j, cell_size/2,
                          cell_size/2, cell_size);
    }
  }
}

// __________________________________________________________________
double calc_displacement(const matrix &pressure,
              const matrix &Ic,
              double y, double x,
              double a, double b, double cell) {
  double res = 0;

  for (int i = 0; i < Ic.shape[0]; i+=1) {
    for (int j = 0; j < Ic.shape[1]; j+=1) {
      double xj = (x*cell)-(j*cell);
      double yi = (y*cell)-(i*cell);
      res += calcBoussinesq(a, b, xj, yi) * pressure(i, j);
    }
  }
  return res;
}

void calc_displacement(const matrix &pF, double cS, double fS,
                         matrix &cD) { // NOLINT
  int shape = pF.shape[0];
  // #pragma omp parallel for simd
  for (int i = 0; i < shape; ++i) {
    for (int j = 0; j < shape; ++j) {
      // inner loop
      for (int k = 0; k < shape; ++k) {
        for (int l = 0; l < shape; ++l) {
          cD(i, j) += calcBoussinesq(i-k, j-l,  cS, cS,
                          fS, fS)*pF(k, l);
        }
      }
    }
  }
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
void initializeStylusArray(matrix &st, int t) { // NOLINT
  for (int i = 1; i <= 2*t; i++) {
    st(i-1, 0) = 1;
    for (int j = 1; j <= 2*t; j++) {
      if (j != i) {
        double divider = (2.0*t-2.0*j)+1;
        double divisor = (2.0*i-2.0*j);
        st(i-1, 0) *= (divider/ divisor);
      }
    }
  }
}

// __________________________________________________________________
void calcCoarsePressure(const std::vector<int>& qs,
                        std::vector<matrix>& pFVec, // NOLINT
                        std::vector<matrix>& cDVec, // NOLINT
                        int t, const matrix& st) {
  int coarse = 0;
  for (int level = 0; level <= qs.size()-1; level++) {
    // std::cout << qs[level] << "thats grid size\n";
    matrix pC({qs[level], qs[level]});
    matrix displacement({qs[level], qs[level]});
    for (int m = 0; m < pC.shape[0]; m++) {
      for (int n = 0; n < pC.shape[1]; n++) {
        // 1 <= k/l <= 2*t
        int i = 2*(m-t+1);
        int j = 2*(n-t+1);
        pC(m, n) = 0;
        displacement(m, n) = 0;
        double res = 0;
        if (!boundaryCheck(pFVec[coarse], i, j)) {
          res += 0;
        } else {
          res += pFVec[coarse](i, j);
        }
        bool once = true;
        for (int k = 1; k <= 2*t; k++) {
          int pi = i+2*(k-t)-1;
          int pj2 = j+2*(k-t)-1;
          for (int l = 1; l <= 2*t; l++) {
            int pj = j+2*(l-t)-1;
            res += boundaryCheck(pFVec[coarse], pi, pj) ?
                    st(k-1, 0)*st(l-1, 0)*pFVec[coarse](pi, pj) : 0;
          }
          if (boundaryCheck(pFVec[coarse], pi, j)) {
            res += st(k-1, 0)*pFVec[coarse](pi, j);
          }
          res += (boundaryCheck(pFVec[coarse], i, pj2)) ?
                  st(k-1, 0)*pFVec[coarse](i, pj2): 0;
        }
        pC(m, n) = res;
      }
    }
    pFVec.push_back(pC);
    cDVec.push_back(displacement);
    coarse++;
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
void correctionSteps(matrix& cC, const matrix& st, int mc, int t, // NOLINT
    double fineSizeA, double fineSizeB, double halfSize) {
  int ts = t;
  /*for (int i = -mc; i <= mc; i++) {
    for (int j = -mc; j <= mc; j++) {
      cC(i+mc, j+mc) = 0;
      bool iEven = (i%2) == 0;
      bool jEven = (j%2) == 0;
      double res1 = 0, res2 = 0, res3 = 0;
      double K = calcBoussinesq(i, j, halfSize, halfSize, fineSizeA, fineSizeB);
      for (int k = 1; k <= 2*t; k++) {
        res1 += st(k-1, 0) * calcBoussinesq(i-2*(k-ts)+1, j, halfSize, halfSize, fineSizeA, fineSizeB);
        res2 += st(k-1, 0) * calcBoussinesq(i, j-2*(k-ts)+1, halfSize, halfSize, fineSizeA, fineSizeB);
        for (int l = 1; l <= 2*t; l++) {
          res3 += st(k-1, 0)*st(l-1, 0) * 
            calcBoussinesq(i-2*(k-ts)+1, j-2*(l-ts)+1, halfSize, halfSize, fineSizeA, fineSizeB);
        }
      }
      cC(i+mc, j+mc) += ((!iEven)&jEven)*(K-res1);
      cC(i+mc, j+mc) += (iEven&(!jEven))*(K-res2);
      cC(i+mc, j+mc) += ((!iEven)&(!jEven))*(K-res3);
    }
  }*/
  for (int i = -mc; i < mc; i+=2) {
    for (int j = -mc; j < mc; j+=2) {
      cC(i+mc, j+mc) = 0;
      cC(i+mc, j+1+mc) = 0;
      cC((i+1)+mc, j+mc) = 0;
      cC((i+1)+mc, (j+1)+mc) = 0;
      double k_odd_j = calcBoussinesq(i, j+1, halfSize, halfSize,
                                  fineSizeA, fineSizeB);
      double k_odd_i = calcBoussinesq(i+1, j, halfSize, halfSize,
                                  fineSizeA, fineSizeB);
      double k_odd_i_j = calcBoussinesq(i+1, j+1, halfSize, halfSize,
                                    fineSizeA, fineSizeB);
      double odd_j = 0, odd_i = 0, odd_i_j = 0;
      for (int k = 1; k <= 2*t; k++) {
        odd_j += st(k-1, 0) * calcBoussinesq(i, j+1-2*(k-ts)+1,
                                        halfSize, halfSize,
                                        fineSizeA, fineSizeB);
        odd_i += st(k-1, 0) * calcBoussinesq(i+1-2*(k-ts)+1, j,
                                        halfSize, halfSize,
                                        fineSizeA, fineSizeB);
        for (int l = 1; l <= 2*t; l++) {
          odd_i_j += st(k-1, 0) * st(l-1, 0) *
            calcBoussinesq((i+1)-2*(k-ts)+1, (j+1)-2*(l-ts)+1,
                      halfSize, halfSize, fineSizeA, fineSizeB);
        }
      }
      cC(i+mc, (j+1)+mc) += k_odd_j - odd_j;
      cC((i+1)+mc, j+mc) += k_odd_i - odd_i;
      cC((i+1)+mc, (j+1)+mc) += k_odd_i_j - odd_i_j;
    }
  }
  for (int i = -mc; i < mc; i+=2) {
    double K1 = calcBoussinesq(mc, i+1, halfSize, halfSize,
                              fineSizeA, fineSizeB);
    cC(2*mc, i+1) = 0;
    double res1 = 0;
    for (int k = 1; k <= 2*t; k++) {
      res1 += st(k-1, 0) * calcBoussinesq(mc, (i+1)-2*(k-ts)+1,
                                      halfSize, halfSize,
                                      fineSizeA, fineSizeB);
    }
    cC(2*mc, i+1) += K1 - res1;
  }
  // writeToFile(cC, "./tests/cc2");
  // cC(2*mc, 2*mc) = 0;
}

// __________________________________________________________________
void applyCorrection(matrix &cD, const matrix cC, // NOLINT
                      const matrix Ip, // NOLINT
                      int t) {
  #pragma omp parallel for simd
  for (int i = 0; i < cD.shape[0]; i++) {
    for (int j = 0; j < cD.shape[1]; j++) {
      cD(i, j) += correctionHelper(cC, Ip, t, i, j);
    }
  }
}

// __________________________________________________________________
double correctionHelper(const matrix& cC, const matrix& Ip, int t,
    int i, int j) {
  int mc = (cC.shape[0]-1)/2;
  double res = 0;
  for (int k = -mc; k <= mc; k++) {
    for (int l = -mc; l <= mc; l++) {
      int pi = 2*(i-t+1)-k;
      int pj = 2*(j-t+1)-l;
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
  int shape = fD.shape[0]-1;
  int evenGrid = (shape%2) == 0;
  if (evenGrid) shape -= 1;
  for (int i = 0; i < shape; i+=2) {
    for (int j = 0; j < shape; j+=2) {
      int m = (i + 1) / 2 + t - 1;
      int n = (j + 1) / 2 + t - 1;
      int m_1 = (i + 1 + 1) / 2 + t - 1;
      int n_1 = (j + 1 + 1) / 2 + t - 1;

      double fDoddi = 0;
      double fDoddj = 0;
      double fDoddij = 0;
      for (int k = 1; k <= 2*t; k++) {
        int pm_1 = m_1+k-t-1;
        int pn_1 = n_1+k-t-1;

        // fD(i+1, j) #TOTEST
        fDoddi += st(k-1, 0) * cD(pm_1, n);
        // fD(i, j+1)
        fDoddj += (st(k-1, 0) * cD(m, pn_1));

        for (int l = 1; l <= 2*t; l++) {
          int pn_1 = n_1+l-t-1;
          // fD(i+1, j+1)
          fDoddij += st(k-1, 0) * st(l-1, 0) * cD(pm_1, pn_1);
        }
      }
      fD(i+1, j) += fDoddi;
      fD(i, j+1) += fDoddj;
      fD(i+1, j+1) += fDoddij;
      fD(i, j) = cD(m, n);
    }
  }
  if (evenGrid) {
    for (int i = 0; i < shape+1; i+=2) {
      for (int k = 1; k <= 2*t; k++) {
        int m = (i + 1 + 1) / 2 + t - 1;
        int n_1 = (shape + 1) / 2 + t - 1;
        int pn_1 = n_1+k-t-1;
        fD(shape, i+1) += (st(k-1, 0) * cD(m, pn_1));
      }
    }
  }
  /*for (int i = 0; i < fD.shape[0]; i++) {
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
          if ((!iEven)&jEven & (boundaryCheck(cD, pm ,n))) {
            fD(i, j) += st(k-1, 0) * cD(pm, n) *(l==1);
          }
          if (iEven&(!jEven) & (boundaryCheck(cD, m, pn))) {
            fD(i, j) += (st(l-1, 0) * cD(m, pn)) * once;
          } else if ((!iEven)&(!jEven) & (boundaryCheck(cD, pm, pn))) {
            fD(i, j) += st(k-1, 0) * st(l-1, 0) * cD(pm, pn);
          }
        }
        once = false;
      }
      if (boundaryCheck(cD, (i+1)/2+t-1, (j+1)/2+t-1)) {
        fD(i, j) += (iEven&jEven)*cD((i+1)/2+t-1, (j+1)/2+t-1);
      }
    }
  }*/
}

// __________________________________________________________________
void secondCorrectionStep(int mc, const matrix& st, double fineSizeA,
                          double fineSizeB, double hS,
                          const matrix& pF, matrix& cD, // NOLINT
                          const std::vector<matrix> &cCVec) {
  int shape = cD.shape[0];
  // std::cout << shape << std::endl;
  bool evenGrid = false;
  if ((shape%2) == 1) {
    shape -= 1;
    evenGrid = true;
  }
  int t = mc/2;
  for (int i = 0; i < shape; i+=2) {
    for (int j = 0; j < shape; j+=2) {
      double oddi = 0;
      double oddj = 0;
      double oddij = 0;
      for (int k = -mc; k <= mc; k++) {
        for (int l = -mc; l <= mc; l++) {
          int pi = i-k;
          int pj = j-l;
          if (boundaryCheck(pF, pi+1, pj+1)) {
            oddj += cCVec[1](k+mc, l+mc)*pF(pi, pj+1);
            oddi += cCVec[0](k+mc, l+mc)*pF(pi+1, pj);
            oddij += cCVec[2](k+mc, l+mc)*pF(pi+1, pj+1);
          }
        }
      }
      cD(i+1, j) += oddi;
      cD(i, j+1) += oddj;
      cD(i+1, j+1) += oddij;
    }
  }
  /*for (int i = 0; i <= mc; i+=2) {
    for (int j = 0; j <= mc; j+=2) {
      double oddi = 0;
      double oddj = 0;
      double oddij = 0;
      for (int k = i; k <= mc; k++) {
        for (int l = j; l <= mc; l++) {
          int pi = i-k;
          int pj = j-l;
          if (boundaryCheck(pF, pi+1, pj+1)){
            oddj += cCVec[1](k+mc, l+mc)*pF(pi, pj+1);
            oddi += cCVec[0](k+mc, l+mc)*pF(pi+1, pj);
            oddij += cCVec[2](k+mc, l+mc)*pF(pi+1, pj+1);
          }
        }
      }
      cD(i+1, j) += oddi;
      cD(i, j+1) += oddj;
      cD(i+1, j+1) += oddij;

    }
  }*/
  if (evenGrid) {
    for (int j = 0; j < shape; j+=2) {
      // std:: cout << "i:" << shape << " j:" << j << std::endl;
      for (int k = -mc; k <= mc; k++) {
        for (int l = -mc; l <= mc; l++) {
          int pi = shape-k;
          int pj = j-l;
          if (boundaryCheck(pF, pi, pj+1)) {
            cD(shape, j+1) += cCVec[1](k+mc, l+mc)*pF(pi, pj+1);
          }
        }
      }
    }
  }
  /*int shape = cD.shape[0];
  for (int i = 0; i < shape; i++) {
    for (int j = 0; j < shape; j++) {
      double correction = 0;
      bool iEven = (i%2) == 0;
      bool jEven = (j%2) == 0;
      bool bounds = true;
      for (int k = -mc; k <= mc; k++) {
        for (int l = -mc; l <= mc; l++) {
          int pi = i-k;
          int pj = j-l;
          if (boundaryCheck(pF, pi, pj)){
            correction += ((!iEven)&jEven) * cCVec[0](k+mc, l+mc)*pF(pi, pj);
            correction += (iEven&(!jEven)) *  cCVec[1](k+mc, l+mc)*pF(pi, pj);
            correction += ((!iEven)&(!jEven)) *  cCVec[2](k+mc, l+mc)*pF(pi, pj);
          }
        }
      }
      cD(i, j) += correction;
    }
  }*/
}


// __________________________________________________________________
void createCorrectionArrays(std::vector<matrix> &cCVec, // NOLINT
                            const matrix &st, double hS, // NOLINT
                            double fineSizeA, double fineSizeB,
                            int mc) {
  int tempMc = 2*mc+1;
  int t = mc/2;
  matrix cC2({tempMc, tempMc});
  matrix cC3({tempMc, tempMc});
  matrix cC4({tempMc, tempMc});
  for (int i = -mc; i <= mc; i++) {
    for (int j = -mc; j <= mc; j++) {
      double K = calcBoussinesq(i, j, hS, hS, fineSizeA, fineSizeB);
      double res = 0, res2 = 0, res3 = 0, res4 = 0;
      for (int k = 1; k <= t*2; k++) {
        for (int l = 1; l <= t*2; l++) {
          res4 +=  st(k-1, 0)*st(l-1, 0) *
                calcBoussinesq(i+2*(k-t)-1, j+2*(l-t)-1, hS, hS,
                          fineSizeA, fineSizeB);
        }
        res2 += st(k-1, 0) * calcBoussinesq(i+2*(k-t)-1, j, hS, hS,
                                        fineSizeA, fineSizeB);
        res3 += st(k-1, 0) * calcBoussinesq(i, j+2*(k-t)-1, hS, hS,
                                        fineSizeA, fineSizeB);
      }
      cC2(i+mc, j+mc) = K-res2;
      cC3(i+mc, j+mc) = K-res3;
      cC4(i+mc, j+mc) = K-res4;
    }
  }
  cCVec.push_back(cC2);
  cCVec.push_back(cC3);
  cCVec.push_back(cC4);
}

// __________________________________________________________________
void copyPressureArray(matrix& p, const matrix& tempP) {
  for (int i = 0; i < tempP.shape[0]; i++) {
    for (int j = 0; j < tempP.shape[1]; j++) {
      p(i, j) = tempP(i, j);
    }
  }  
}

// __________________________________________________________________
void calculateGmn(matrix &Gmn, double dx, double dy) { // NOLINT
  for (int i = 0; i < (Gmn.shape[0]-1)/2; i++) {
    for (int j = 0; j < (Gmn.shape[0]-1)/2; j++) {
      Gmn(i, j) = calcBoussinesq(i, j, dx, dy, dx, dy);
      Gmn(i, j+((Gmn.shape[0]-1)/2)) = calcBoussinesq(i, -j, dx, dy, dx, dy);
      Gmn(i+((Gmn.shape[0]-1)/2), j) = calcBoussinesq(-i, j, dx, dy, dx, dy);
      Gmn(i+((Gmn.shape[0]-1)/2), j+((Gmn.shape[0]-1)/2)) =
        calcBoussinesq(-i, -j, dx, dy, dx, dy);
      // std::cout << Gmn(i, j) << std::endl;
    }
  }
}

// __________________________________________________________________
void transformGmnP(int Nx, int Ny, matrix& Gmn, cMatrix& Gmn_tild,
                  matrix& p, cMatrix& p_tild) {
  fftw_plan p1;
  p1 = fftw_plan_dft_r2c_2d(Nx, Ny*2, Gmn.data.data(),
      reinterpret_cast<fftw_complex*>(Gmn_tild.data.data()), FFTW_ESTIMATE);
  // output array needs to be 2*nx / (ny*2/2)-1

  fftw_plan p2;
  p2 = fftw_plan_dft_r2c_2d(Nx, Ny*2, p.data.data(),
      reinterpret_cast<fftw_complex*>(p_tild.data.data()), FFTW_ESTIMATE);

  fftw_execute(p1);
  fftw_execute(p2);

  fftw_destroy_plan(p1);
  fftw_destroy_plan(p2);
}

// __________________________________________________________________
void multiplyTransformed(cMatrix& Gmn_tild, cMatrix& Umn_tild,
                        cMatrix& p_tild) {
  for (int i = 0; i < Gmn_tild.shape[0]; i++) {
    for (int j = 0; j < Gmn_tild.shape[1]; j++) {
      Umn_tild(i, j) = Gmn_tild(i, j)*p_tild(i, j);
    }
  }
}

// __________________________________________________________________
void transformToReal(cMatrix& Umn_tild, matrix& Umn, int Nx, int Ny) {
  fftw_plan p3;
  p3 = fftw_plan_dft_c2r_2d(Nx, Ny*2,
                            reinterpret_cast<fftw_complex*>
                            (Umn_tild.data.data()),
                            Umn.data.data(), FFTW_ESTIMATE);
  fftw_execute(p3);
  fftw_destroy_plan(p3);
}

// __________________________________________________________________
void writeToResultArray(const matrix& Umn, matrix& Umn_res, 
                        int Nx, int Ny) {
  int N = Nx*Ny;
  #pragma omp parallel for simd
  for (int i = 1; i < Umn_res.shape[0]; i++) {
    for (int j = 0; j < Umn_res.shape[1]; j++) {
    Umn_res(i, j) = Umn(i, j)/N;
    }
  }
  writeToFile(Umn_res, "whatthe");
}
