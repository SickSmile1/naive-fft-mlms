/* Copyright 2023 <Ilia Fedotov @ Uni Freiburg> */
#include "BoussinesqMlms.h"
#include "Boussinesq.h"
#include <algorithm>
#include <string>
#include <vector>

// __________________________________________________________________
matrix initializeStylusArray(int t) { // NOLINT
  matrix st({2*t, 1});
  for (int i = 1; i <= 2*t; i++) {
    st(i-1, 0) = 1;
    for (int j = 1; j <= 2*t; j++) {
      if (j != i) {
        st(i-1, 0) *= ((2.0*t-2.0*j)+1)/
                        (2.0*i-2.0*j);
      }
    }
  }
  return st;
}

// __________________________________________________________________
void initializeStack(matrix &st, const int t, const matrix Ip, // NOLINT
                    const matrix kM, // NOLINT
                    std::vector<matrix>& pfVec, // NOLINT
                    std::vector<matrix>& cDVec) { // NOLINT
  const int grid = kM.shape[0];
  pfVec.push_back(Ip);
  cDVec.push_back(kM);
  double qLevel = std::ceil(std::log2(grid * grid)/4);
  int q;
  for (int length = 0; length < qLevel-1; length++) {
    if (length == 0) {
      q = grid / 2 + 2*t - 1;
    } else {
      q = q / 2 +  2*t - 1;
    }
    matrix temp({q, q});
    pfVec.push_back(temp);
    cDVec.push_back(temp);
  }
}

// __________________________________________________________________
void calcCoarsePressure(std::vector<matrix>& pFVec, // NOLINT
                        const matrix& st) {
  int t = st.shape[0]/2;
  for (int level = 0; level <= pFVec.size()-2; level++) {
    matrix &pC = pFVec[level+1];
    int bound = pFVec[level].shape[0];
    for (int m = 0; m < pC.shape[0]; m++) {
      for (int n = 0; n < pC.shape[1]; n++) {
        int i = 2*(m-t+1);
        int j = 2*(n-t+1);
        pC(m, n) = 0;
        double res = 0;
        if (!bCheck(bound, i, j)) {
          res += 0;
        } else {
          res += pFVec[level](i, j);
        }
        for (int k = 1; k <= 2*t; k++) {
          int pi = i+2*(k-t)-1;
          int pj2 = j+2*(k-t)-1;
          for (int l = 1; l <= 2*t; l++) {
            int pj = j+2*(l-t)-1;
            res += (bCheck(bound, pi, pj) && pFVec[level](pi, pj)!=0) ?
                    st(k-1, 0)*st(l-1, 0)*pFVec[level](pi, pj) : 0;
          }
          if (bCheck(bound, pi, j) && pFVec[level](pi, j)!=0) {
            res += st(k-1, 0)*pFVec[level](pi, j);
          }
          res += (bCheck(bound, i, pj2) && pFVec[level](i, pj2) != 0) ?
                  st(k-1, 0)*pFVec[level](i, pj2): 0;
        }
        pC(m, n) = res;
      }
    }
  }
}


// __________________________________________________________________
bool bCheck(int b, int i, int j) { // NOLINT
  bool res = (i < b) & (j < b) & (i >= 0) & (j >= 0);
  return res;
}

// __________________________________________________________________
void correctionSteps(matrix& cC, const matrix& st, int mc, int t, // NOLINT
    double fineSize, double halfSize) {
  int ts = t;
  for (int i = -mc; i <= mc; i++) {
    for (int j = -mc; j <= mc; j++) {
      cC(i+mc, j+mc) = 0;
      bool iEven = (i%2) == 0;
      bool jEven = (j%2) == 0;
      double res1 = 0, res2 = 0, res3 = 0;
      double K = calcBoussinesq(i, j, halfSize, halfSize, fineSize, fineSize);
      if (iEven && jEven) {
        cC(i+mc, j+mc) = 0;
      } else if (!iEven && jEven) {
        for (int k = 1; k <= 2*t; k++) {
          res1 += st(k-1, 0) * calcBoussinesq(i-2*(k-ts)+1, j,
                    halfSize, halfSize, fineSize, fineSize);
        }
        cC(i+mc, j+mc) += K-res1;
      } else if (iEven && !jEven) {
        for (int k = 1; k <= 2*t; k++) {
          res2 += st(k-1, 0) * calcBoussinesq(i, j-2*(k-ts)+1,
                    halfSize, halfSize, fineSize, fineSize);
        }
        cC(i+mc, j+mc) += K-res2;
      } else if (!iEven && !jEven) {
        for (int k = 1; k <= 2*t; k++) {
          for (int l = 1; l <= 2*t; l++) {
            res3 += st(k-1, 0)*st(l-1, 0) *
              calcBoussinesq(i-2*(k-ts)+1, j-2*(l-ts)+1, halfSize,
              halfSize, fineSize, fineSize);
          }
        }
        cC(i+mc, j+mc) += K-res3;
      }
    }
  }
}


// __________________________________________________________________
void applyCorrection(matrix &cD, const matrix cC, // NOLINT
                      const matrix Ip, // NOLINT
                      int t, int mc) {
  #pragma omp parallel for simd
  for (int i = 0; i < cD.shape[0]; i++) {
    for (int j = 0; j < cD.shape[1]; j++) {
      cD(i, j) += correctionHelper(cC, Ip, t, i, j, mc);
    }
  }
}

// __________________________________________________________________
double correctionHelper(const matrix& cC, const matrix& Ip, int t,
    int i, int j, int mc) {
  // int mc = (cC.shape[0]-1)/2;
  double res = 0;
  int bound = Ip.shape[0];
  for (int k = -mc; k <= mc; k++) {
    for (int l = -mc; l <= mc; l++) {
      int pi = 2*(i-t+1)-k;
      int pj = 2*(j-t+1)-l;
      if (bCheck(bound, pi, pj)) {
        if(Ip(pi, pj)!=0) {
          res += cC(k+mc, l+mc) * Ip(pi, pj);
        }
      }
    }
  }
  return res;
}

// __________________________________________________________________
void interpolateGrid(matrix &fD, const matrix cD, const matrix st) { // NOLINT
  int t = st.shape[0]/2;
  int bound = cD.shape[0];
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
          if ((!iEven)&jEven & (bCheck(bound, pm, n))) {
            fD(i, j) += st(k-1, 0) * cD(pm, n) *(l == 1);
          }
          if (iEven&(!jEven) & (bCheck(bound, m, pn))) {
            fD(i, j) += (st(l-1, 0) * cD(m, pn)) * once;
          } else if ((!iEven)&(!jEven) & (bCheck(bound, pm, pn))) {
            fD(i, j) += st(k-1, 0) * st(l-1, 0) * cD(pm, pn);
          }
        }
        once = false;
      }
      if (bCheck(bound, (i+1)/2+t-1, (j+1)/2+t-1)) {
        fD(i, j) += (iEven&jEven)*cD((i+1)/2+t-1, (j+1)/2+t-1);
      }
    }
  }
}

// __________________________________________________________________
void secondCorrectionStep(const matrix& st,
                          double hS, const matrix& pF, matrix& cD, // NOLINT
                          const std::vector<matrix> &cCVec, int mc) {
  int shape = cD.shape[0];
  int bound = pF.shape[0];
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
          if (bCheck(bound, pi, pj)){
            if (pF(pi,pj)!=0) {
              correction += ((!iEven)&jEven) * cCVec[0](k+mc, l+mc)*pF(pi, pj);
              correction += (iEven&(!jEven)) *  cCVec[1](k+mc, l+mc)*pF(pi, pj);
              correction += ((!iEven)&(!jEven)) *  cCVec[2](k+mc, l+mc)*pF(pi, pj);
            }
          }
        }
      }
      cD(i, j) += correction;
    }
  }
}


// __________________________________________________________________
void createCorrectionArrays(std::vector<matrix> &cCVec, // NOLINT
                            const matrix &st, double hS, // NOLINT
                            double fineSize, int mc) {
  int t = st.shape[0]/2;
  int tempMc = 2*mc+2;
  matrix cC2({tempMc, tempMc});
  matrix cC3({tempMc, tempMc});
  matrix cC4({tempMc, tempMc});
  for (int i = -(mc); i <= mc; i++) {
    for (int j = -(mc); j <= mc; j++) {
      double K = calcBoussinesq(i, j, hS, hS, fineSize, fineSize);
      double res = 0, res2 = 0, res3 = 0, res4 = 0;
      for (int k = 1; k <= t*2; k++) {
        for (int l = 1; l <= t*2; l++) {
          res4 +=  st(k-1, 0)*st(l-1, 0) *
                calcBoussinesq(i+2*(k-t)-1, j+2*(l-t)-1, hS, hS,
                          fineSize, fineSize);
        }
        res2 += st(k-1, 0) * calcBoussinesq(i+2*(k-t)-1, j, hS, hS,
                                        fineSize, fineSize);
        res3 += st(k-1, 0) * calcBoussinesq(i, j+2*(k-t)-1, hS, hS,
                                        fineSize, fineSize);
      }
      cC2(i+(mc), j+(mc)) = K-res2;
      cC3(i+(mc), j+(mc)) = K-res3;
      cC4(i+(mc), j+(mc)) = K-res4;
    }
  }
  cCVec.push_back(cC2);
  cCVec.push_back(cC3);
  cCVec.push_back(cC4);
}

// __________________________________________________________________
matrix BoussinesqMlms(double size, int grid1, int t) {
  // const double size = 2;
  // const double size_p = 1;
  const double pressure{1.};

  double size_p{size/2};

  double fineSizeA{size / grid1};

  // result array
  matrix kM({grid1, grid1});
  // pressure array
  matrix Ip({grid1, grid1});

  double lower_b = (grid1)/2. - (size_p/fineSizeA)/2.;
  double upper_b = (grid1)/2. + (size_p/fineSizeA)/2.;

  initializePressureArray(Ip, lower_b, upper_b, pressure);

  // int t = 4;
  int mc = std::max(0.7*t*std::pow(grid1,1./t)-1, t*1.);

  matrix st = initializeStylusArray(t);

  // vectors of different grid sizes for pressure->pfVec / displacement -> cDvec
  std::vector<matrix> pfVec;
  std::vector<matrix> cDVec;

  // fill vectors with empty arrays for different grid sizes before calculation
  initializeStack(st, t, Ip, kM, pfVec, cDVec);
  double d{pfVec.size()-1.};

  // std::cout << "grid: " << pfVec[d-1].shape[0] << " grid2: " << pfVec[1].shape[0] << std::endl;

  double coarseSize{fineSizeA*pow(2, d)};
  std::vector<matrix> cCVec;
  calcCoarsePressure(pfVec, st);
  calc_displacement(pfVec[d], coarseSize, fineSizeA, cDVec[d]);
  for (int i = 0; i < pfVec.size()-1; i++) {
    double hS = fineSizeA*pow(2, d-i-1);
    // std::cout << hS <<std::endl;
    int temp_mc = (mc*2)+1;
    matrix cC({temp_mc, temp_mc});
    correctionSteps(cC, st, mc, t, fineSizeA, hS);
    applyCorrection(cDVec[d-i], cC, pfVec[d-i-1], t, mc);
    interpolateGrid(cDVec[d-i-1], cDVec[d-i], st);
    cCVec.reserve(3);
    createCorrectionArrays(cCVec, st, hS, fineSizeA, mc);
    secondCorrectionStep(st, hS,
                         pfVec[d-i-1], cDVec[d-i-1], cCVec, mc);
    // writeToFile(cDVec[d-i-1], "results/partial_"+std::to_string(int(d-i-1)));
  }
  return cDVec[0];
}


