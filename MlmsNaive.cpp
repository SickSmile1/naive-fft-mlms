/* Copyright 2023 <Ilia Fedotov @ Uni Freiburg> */
#include <fstream>
#include <iostream>
#include <cmath>
#include <vector>
#include "./Mlms.h"
#include "./MlmsWriteToFile.h"

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
double calculate(double a, double b, double x, double y) {
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
double calculate(int i, int j, double dxc, double dyc,
    double dxf, double dyf) {
  return calculate(dxf/2, dyf/2, i*dxc, j*dyc);
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
      res += calculate(a, b, xj, yi) * pressure(i, j);
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
void initializeStylusArray(matrix &st, int t) { // NOLINT
  initializeMultiplicationArray(st);
  // calculate transfer stylus (8)
  for (int i = 1; i <= 2*t; i++) {
    for (int j = 1; j <= 2*t; j++) {
      if (j != i) {
        // std::cout << st(i, 0) << " : ";
        double divider = (2.0*t-2.0*j)+1;
        double divisor = (2.0*i-2.0*j);
        st(i-1, 0) *= (divider/ divisor);
        // std::cout << st(i, 0) << std::endl;
      }
    }
  }
}

// __________________________________________________________________
void calcCoarsePressure(const std::vector<int>& qs,
                        std::vector<matrix>& pFVec,
                        std::vector<matrix>& cDVec,
                        int t, const matrix& st) {
  int coarse = 0;
  for (int level = 0; level <= qs.size()-1; level++) {
    std::cout << qs[level] << "thats grid size\n";
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
    // writeToFile(pC, "./tests/coarsened_"+std::to_string(qs[level]));
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
void myBreakpoint() {}

// __________________________________________________________________
void correctionSteps(matrix& cC, const matrix& st, int mc, int t, // NOLINT
    double fineSizeA, double fineSizeB, double halfSize) {
  int ts = t;
  for (int i = -mc; i <= mc; i++) {
    for (int j = -mc; j <= mc; j++) {
      cC(i+mc, j+mc) = 0;
      bool iEven = (i%2) == 0;
      bool jEven = (j%2) == 0;
      double res1 = 0, res2 = 0, res3 = 0;
      double K = calculate(i, j, halfSize, halfSize, fineSizeA, fineSizeB);
      for (int k = 1; k <= 2*t; k++) {
        res1 += st(k-1, 0) * calculate(i-2*(k-ts)+1, j, halfSize, halfSize, fineSizeA, fineSizeB);
        res2 += st(k-1, 0) * calculate(i, j-2*(k-ts)+1, halfSize, halfSize, fineSizeA, fineSizeB);
        for (int l = 1; l <= 2*t; l++) {
          res3 += st(k-1, 0)*st(l-1, 0) * 
            calculate(i-2*(k-ts)+1, j-2*(l-ts)+1, halfSize, halfSize, fineSizeA, fineSizeB);
        }
      }
      /*for (int l = 1; l <= 2*t; l++) {
        res2 += st(l-1, 0) * calculate(i, j-2*(l-ts)+1, halfSize, halfSize, fineSizeA, fineSizeB);
      }*/
      //for (int k = 1; k <= 2*t; k++) {
      //}
      cC(i+mc, j+mc) += ((!iEven)&jEven)*(K-res1);
      cC(i+mc, j+mc) += (iEven&(!jEven))*(K-res2);
      cC(i+mc, j+mc) += ((!iEven)&(!jEven))*(K-res3);
      // std::cout << cC(i+mc, j+mc) << std:: endl;
    }
  }
  /*writeToFile(cC, "./tests/cc1");
  for (int i = -mc; i < mc; i+=2) {
    for (int j = -mc; j < mc; j+=2) {
      //std::cout << "j: " << j+1 << " i:" << i << std:: endl;
      //std::cout << "j: " << j << " i:" << i+1 << std:: endl;
      //std::cout << "j: " << j+1 << " i:" << i+1 << std:: endl;
      cC(i+mc, j+mc) = 0;
      cC(i+mc, j+1+mc) = 0;
      cC((i+1)+mc, j+mc) = 0;
      cC((i+1)+mc, (j+1)+mc) = 0;
      double k_odd_j = calculate(i, j+1, halfSize, halfSize, fineSizeA, fineSizeB);
      double k_odd_i = calculate(i+1, j, halfSize, halfSize, fineSizeA, fineSizeB);
      double k_odd_i_j = calculate(i+1, j+1, halfSize, halfSize, fineSizeA, fineSizeB);
      double odd_j = 0, odd_i = 0, odd_i_j = 0;
      
      for (int k = 1; k <= 2*t; k++) {
        odd_j += st(k-1, 0) * calculate(i-2*(k-ts)+1, j+1, halfSize, halfSize, fineSizeA, fineSizeB);
        odd_i += st(k-1, 0) * calculate(i+1, j-2*(k-ts)+1, halfSize, halfSize, fineSizeA, fineSizeB);
        for (int l = 1; l <= 2*t; l++) {
          odd_i_j += st(k-1, 0) * st(l-1, 0) *
            calculate((i+1)-2*(k-ts)+1, (j+1)-2*(l-ts)+1,
                      halfSize, halfSize, fineSizeA, fineSizeB);
        }
      }
      cC(i+mc, (j+1)+mc) += k_odd_j - odd_j;
      cC((i+1)+mc, j+mc) += k_odd_i - odd_i;
      cC((i+1)+mc, (j+1)+mc) += k_odd_i_j - odd_i_j;
      //std::cout << cC(i+1+mc, j+mc) << std:: endl;
      //std::cout << cC(i+mc, j+1+mc) << std:: endl;
      //std::cout << cC(i+1+mc, j+1+mc) << std:: endl;
    }
  }
  for (int i = -mc; i < mc; i+=2) {
    double K1 = calculate(mc, i+1, halfSize, halfSize, fineSizeA, fineSizeB);
    cC(2*mc, i+1) = 0;
    double res1 = 0;
    for (int k = 1; k <= 2*t; k++) {
      res1 += st(k-1, 0) * calculate(mc, (i+1)-2*(k-ts)+1, halfSize, halfSize, fineSizeA, fineSizeB);
    }
    cC(2*mc, i+1) += K1 - res1;
  }
  writeToFile(cC, "./tests/cc2");*/
  // cC(2*mc, 2*mc) = 0;
}

// __________________________________________________________________
void applyCorrection(matrix &cD, const matrix cC, 
                      const matrix Ip, // NOLINT
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
      // fD(i, j) += (iEven&jEven)*res;
      // fD(i, j) += (iEven&jEven)*res1;
      // fD(i, j) += (iEven&jEven)*res2;
      if (boundaryCheck(cD, (i+1)/2+t-1, (j+1)/2+t-1)) {
        fD(i, j) += (iEven&jEven)*cD((i+1)/2+t-1, (j+1)/2+t-1);
      }
    }
  }
}

// __________________________________________________________________
void secondCorrectionStep(int mc, const matrix& st, double fineSizeA,
                          double fineSizeB, double hS,
                          const matrix& pF, matrix& cD,
                          const std::vector<matrix> &cCVec) {
  int shape = cD.shape[0];
  int t = mc/2;
  std::cout << "shape cd: " << shape << " shape pf: " << pF.shape[0] << std::endl;
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
             // cC2(k+mc, l+mc)*pF(pi, pj);
            correction += (iEven&(!jEven)) *  cCVec[1](k+mc, l+mc)*pF(pi, pj);
            // cC3(k+mc, l+mc)*pF(pi, pj);
            correction += ((!iEven)&(!jEven)) *  cCVec[2](k+mc, l+mc)*pF(pi, pj);
            // cC4(k+mc, l+mc)*pF(pi, pj);
          //} else }
            //bounds = false;
          }
        }
      }
      //if(!bounds) {
       // std::cout << "i: " << i << " j: " << j << std::endl;
      //  bounds = true;
      //}
      cD(i, j) += correction;
    }
  }
}


// __________________________________________________________________
void createCorrectionArrays(std::vector<matrix> &cCVec, 
                            const matrix &st, double hS, 
                            double fineSizeA, double fineSizeB,
                            int mc) {
  int tempMc = 2*mc+1;
  int t = mc/2;
  matrix cC2({tempMc, tempMc});
  matrix cC3({tempMc, tempMc});
  matrix cC4({tempMc, tempMc});
  for (int i = -mc; i <= mc; i++) {
    for (int j = -mc; j <= mc; j++) {
      double K = calculate(i, j, hS, hS, fineSizeA, fineSizeB);
      double res = 0, res2 = 0, res3 = 0, res4 = 0;
      for (int k = 1; k <= t*2; k++) {
        for (int l = 1; l <= t*2; l++) {
          res4 +=  st(k-1, 0)*st(l-1,0) * 
                calculate(i+2*(k-t)-1, j+2*(l-t)-1, hS, hS, fineSizeA, fineSizeB);
        }
        res2 += st(k-1, 0) * calculate(i+2*(k-t)-1, j, hS, hS, fineSizeA, fineSizeB);
        res3 += st(k-1, 0) * calculate(i, j+2*(k-t)-1, hS, hS, fineSizeA, fineSizeB);
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
