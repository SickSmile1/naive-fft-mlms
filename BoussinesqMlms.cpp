/* Copyright 2023 <Ilia Fedotov @ Uni Freiburg> */
#include "BoussinesqMlms.h"
#include "Boussinesq.h"
#include <algorithm>


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
  double qLevel = std::log2(grid * grid)/2-1;
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
void old_calcCoarsePressure(std::vector<matrix>& pFVec, // NOLINT
                        const matrix& st) {
  int t = st.shape[0]/2;
  for (int level = 0; level <= pFVec.size()-2; level++) {
    matrix &pC = pFVec[level+1];
    for (int m = 0; m < pC.shape[0]; m++) {
      for (int n = 0; n < pC.shape[1]; n++) {
        // establish bigger boundaries to interpolate pressure on the
        // edges of the grid
        int i = 2*(m-t+1);
        int j = 2*(n-t+1);
        // determine loop boundaries to avoid boundary checks
        int lb_i = 1, lb_j =1;
        int ub_i = 2*t, ub_j = 2*t;
        if (i < 2*t) {
          lb_i = t*2-m;
        }
        if (i > (pFVec[level].shape[0]-2*t)) {
          ub_i = (pC.shape[0]-m);
        }
        if (j < 2*t) {
          lb_j = t*2-m;
        }
        if (j > (pFVec[level].shape[0]-2*t)) {
          ub_j = (pC.shape[0]-n);
        }
        // interpolate pressure for m, n from i, j of coarse pressure
        // array
        pC(m, n) = 0;
        double res = 0;
        if (!boundaryCheck(pFVec[level], i, j)) {
          res += 0;
        } else {
          res += pFVec[level](i, j);
        }
        if ((j > 0)&(j < pFVec[level].shape[0])) {
          for (int k = lb_i; k <= ub_i; k++) {
            int pi = i+2*(k-t)-1;
            res += st(k-1, 0)*pFVec[level](pi, j);
          }
        }
        if ((i > 0) & (i < pFVec[level].shape[1])) {
          for (int k = lb_j; k <= ub_j; k++) {
            int pj = j+2*(k-t)-1;
            res += st(k-1, 0)*pFVec[level](i, pj);
          }
        }
        for (int k = lb_i; k <= ub_i; k++) {
          for (int l = lb_j; l <= ub_j; l++) {
            int pi = i+2*(k-t)-1;
            int pj = j+2*(l-t)-1;
            res += st(k-1, 0)*st(l-1,0)*pFVec[level](pi, pj);
          }
        }
        pC(m, n) = res; // i+2*(1-t)-1;
      }
    }
  }
}

// __________________________________________________________________
void calcCoarsePressure(std::vector<matrix>& pFVec, // NOLINT
                        const matrix& st) {

  int t = st.shape[0]/2;
  for (int level = 0; level <= pFVec.size()-2; level++) {
    matrix &pC = pFVec[level+1];
    for (int m = 0; m < pC.shape[0]; m++) {
      for (int n = 0; n < pC.shape[1]; n++) {
        int i = 2*(m-t+1);
        int j = 2*(n-t+1);
        pC(m, n) = 0;
        double res = 0;
        if (!boundaryCheck(pFVec[level], i, j)) {
          res += 0;
        } else {
          res += pFVec[level](i, j);
        }
        for (int k = 1; k <= 2*t; k++) {
          int pi = i+2*(k-t)-1;
          int pj2 = j+2*(k-t)-1;
          for (int l = 1; l <= 2*t; l++) {
            int pj = j+2*(l-t)-1;
            res += boundaryCheck(pFVec[level], pi, pj) ?
                    st(k-1, 0)*st(l-1, 0)*pFVec[level](pi, pj) : 0;
          }
          if (boundaryCheck(pFVec[level], pi, j)) {
            res += st(k-1, 0)*pFVec[level](pi, j);
          }
          res += (boundaryCheck(pFVec[level], i, pj2)) ?
                  st(k-1, 0)*pFVec[level](i, pj2): 0;
        }
        pC(m, n) = res;
      }
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
void correctionSteps(matrix& cC, const matrix& st, int mc, int t, // NOLINT
    double fineSize, double halfSize) {
  int ts = t;
  /*for (int i = -mc; i <= mc; i++) {
    for (int j = -mc; j <= mc; j++) {
      cC(i+mc, j+mc) = 0;
      bool iEven = (i%2) == 0;
      bool jEven = (j%2) == 0;
      double res1 = 0, res2 = 0, res3 = 0;
      double K = calcBoussinesq(i, j, halfSize, halfSize, fineSizeA, fineSizeB);
      for (int k = 1; k <= 2*t; k++) {
        res1 += st(k-1, 0) * calcBoussinesq(i-2*(k-ts)+1, j, 
                    halfSize, halfSize, fineSizeA, fineSizeB);
        res2 += st(k-1, 0) * calcBoussinesq(i, j-2*(k-ts)+1,
                    halfSize, halfSize, fineSizeA, fineSizeB);
        for (int l = 1; l <= 2*t; l++) {
          res3 += st(k-1, 0)*st(l-1, 0) * 
            calcBoussinesq(i-2*(k-ts)+1, j-2*(l-ts)+1, halfSize,
            halfSize, fineSizeA, fineSizeB);
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
                                  fineSize, fineSize);
      double k_odd_i = calcBoussinesq(i+1, j, halfSize, halfSize,
                                  fineSize, fineSize);
      double k_odd_i_j = calcBoussinesq(i+1, j+1, halfSize, halfSize,
                                    fineSize, fineSize);
      double odd_j = 0, odd_i = 0, odd_i_j = 0;
      for (int k = 1; k <= 2*t; k++) {
        odd_j += st(k-1, 0) * calcBoussinesq(i, j+1-2*(k-ts)+1,
                                        halfSize, halfSize,
                                        fineSize, fineSize);
        odd_i += st(k-1, 0) * calcBoussinesq(i+1-2*(k-ts)+1, j,
                                        halfSize, halfSize,
                                        fineSize, fineSize);
        for (int l = 1; l <= 2*t; l++) {
          odd_i_j += st(k-1, 0) * st(l-1, 0) *
            calcBoussinesq((i+1)-2*(k-ts)+1, (j+1)-2*(l-ts)+1,
                      halfSize, halfSize, fineSize, fineSize);
        }
      }
      cC(i+mc, (j+1)+mc) += k_odd_j - odd_j;
      cC((i+1)+mc, j+mc) += k_odd_i - odd_i;
      cC((i+1)+mc, (j+1)+mc) += k_odd_i_j - odd_i_j;
    }
  }
  for (int i = -mc; i < mc; i+=2) {
    double K1 = calcBoussinesq(mc, i+1, halfSize, halfSize,
                              fineSize, fineSize);
    cC(2*mc, i+1) = 0;
    double res1 = 0;
    for (int k = 1; k <= 2*t; k++) {
      res1 += st(k-1, 0) * calcBoussinesq(mc, (i+1)-2*(k-ts)+1,
                                      halfSize, halfSize,
                                      fineSize, fineSize);
    }
    cC(2*mc, i+1) += K1 - res1;
  }
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

      double fDoddi = 0, fDoddj = 0, fDoddij = 0;
      for (int k = 1; k <= 2*t; k++) {
        int pm_1 = m_1+k-t-1;
        int pn_1 = n_1+k-t-1;

        fDoddi += st(k-1, 0) * cD(pm_1, n);
        fDoddj += (st(k-1, 0) * cD(m, pn_1));

        for (int l = 1; l <= 2*t; l++) {
          int pn_1 = n_1+l-t-1;
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
  /*
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
          if ((!iEven)&jEven & (boundaryCheck(cD, pm, n))) {
            fD(i, j) += st(k-1, 0) * cD(pm, n) *(l == 1);
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
void secondCorrectionStep(const matrix& st,
                          double hS, const matrix& pF, matrix& cD, // NOLINT
                          const std::vector<matrix> &cCVec) {
  int shape = cD.shape[0];
  bool evenGrid = false;
  if ((shape%2) == 1) {
    shape -= 1;
    evenGrid = true;
  }
  int t = st.shape[0]/2;
  int mc = 2*t;
  for (int i = 0; i < shape; i+=2) {
    for (int j = 0; j < shape; j+=2) {
      double oddi = 0, oddj = 0, oddij = 0;
      for (int k = -mc; k <= mc; k++) {
        for (int l = -mc; l <= mc; l++) {
          int pi = i-k;
          int pj = j-l;
          if (boundaryCheck(pF, pi+1, pj+1)) {
            oddij += cCVec[2](k+mc, l+mc)*pF(pi+1, pj+1);
          }
          if (boundaryCheck(pF, pi, pj+1)) {
            oddj += cCVec[1](k+mc, l+mc)*pF(pi, pj+1);
          }
          if (boundaryCheck(pF, pi+1, pj)) {
            oddi += cCVec[0](k+mc, l+mc)*pF(pi+1, pj);
          }
        }
      }
      cD(i+1, j) += oddi;
      cD(i, j+1) += oddj;
      cD(i+1, j+1) += oddij;
    }
  }
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
                            double fineSize) {
  int t = st.shape[0]/2;
  int tempMc = 2*(2*t)+1;
  matrix cC2({tempMc, tempMc});
  matrix cC3({tempMc, tempMc});
  matrix cC4({tempMc, tempMc});
  for (int i = -(2*t); i <= 2*t; i++) {
    for (int j = -(2*t); j <= 2*t; j++) {
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
      cC2(i+(2*t), j+(2*t)) = K-res2;
      cC3(i+(2*t), j+(2*t)) = K-res3;
      cC4(i+(2*t), j+(2*t)) = K-res4;
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
  const double pressure = 1.;

  double size_p = size/2;

  double fineSizeA = size / grid1;

  matrix kM({grid1, grid1});
  matrix Ip({grid1, grid1});

  double lower_b = (grid1)/2. - (size_p/fineSizeA)/2.;
  double upper_b = (grid1)/2. + (size_p/fineSizeA)/2.;

  initializePressureArray(Ip, lower_b, upper_b, pressure);

  // int t = 4;
  int mc = std::max(0.7*t*std::pow(grid1,1./t)-1,t*1.);
  
  matrix st = initializeStylusArray(t);

  std::vector<matrix> pfVec;
  std::vector<matrix> cDVec;

  initializeStack(st, t, Ip, kM, pfVec, cDVec);
  double d = pfVec.size()-1;

  double coarseSize = fineSizeA*pow(2, d);
  std::vector<matrix> cCVec;
  cCVec.reserve(3);
  createCorrectionArrays(cCVec, st, coarseSize, fineSizeA);

  calcCoarsePressure(pfVec, st);
  calc_displacement(pfVec[d], coarseSize, fineSizeA, cDVec[d]);
  for (int i = 0; i < pfVec.size()-1; i++) {
    double hS = fineSizeA*pow(2, d-i-1);
    int temp_mc = (mc*2)+1;
    matrix cC({temp_mc, temp_mc});
    correctionSteps(cC, st, mc, t, fineSizeA, hS);
    applyCorrection(cDVec[d-i], cC, pfVec[d-i-1], t);
    interpolateGrid(cDVec[d-i-1], cDVec[d-i], st);
    secondCorrectionStep(st, hS,
                         pfVec[d-i-1], cDVec[d-i-1], cCVec);
  }
  return cDVec[0];
}


