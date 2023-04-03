
#include "./Mlms.h"
#include <iostream> 
#include <cmath>

// __________________________________________________________________
void initializeMultiplicationArray(matrix &Ic) {
  for (std::size_t i = 0; i < Ic.shape[0]; i++) {
    for (std::size_t j = 0; j < Ic.shape[1]; j++) {
      Ic(i, j) = 1;
    }
  }
}

// __________________________________________________________________
void calc_coarse_pressure(const matrix &fP, const matrix &st, matrix &cP,
                          std::size_t ts) {
  int t = ts;
  for (std::size_t m = 0+t; m < cP.shape[0]; m++) {
    for (std::size_t n = 0+t; n < cP.shape[1]; n++) {
      int i = m-t+1;
      int j = n-t+1;
      cP(m,n) = 0;
      double c_pressure = fP(i, j);
      for (std::size_t k = 1; k < 2 * t; k++) {
        int pi = i+2*(k-t)-1;
        if (pi > 0) {
          c_pressure += st(k, 0)*fP(pi, j);
        }
      }
      for (int l = 1; l < 2 * t; l++) {
        int pj = j+2*(l-t)-1;
        c_pressure += st(l, 0)*fP(i, pj);
      }

      /*for (std::size_t k = 1; k < 2 * t; k++) {
        for (std::size_t l = 1; l < 2 * t; l++) {
          int si = i+2*(k-t)-1;
          int sj = j+2*(l-t)-1;
          c_pressure += st(k, 0) * st(l, 0) * fP(si, sj);
        }
      }*/
      cP(m,n) += c_pressure;
    }
  }
}

int main() {

  const double PI = 3.141592653589793238463;
  
  // initial size and pressure values
  const double size = 1000;
  const std::size_t size_p = 500;
  double pressure = 1.;

  // initial grid size forinitialization
  std::size_t grid1 = 1000;
  std::size_t grid2 = 1000;

  matrix Ic({grid1, grid2});
  matrix Ip({grid1, grid2});

  double lower_b = (size)/2 - (size_p)/2;
  double upper_b = (size)/2 + (size_p)/2;

  initializePressureArray(Ip, lower_b, upper_b, pressure);

  // material moduli and v
  double v1 = 0.;
  double v2 = 0.;
  double E1 = 1.;
  double E2 = 1.;
  
  initializeDisplacementArray(Ic);
  // create K(x,y) as of (2)
  for (std::size_t i = 0; i < Ic.shape[0]; i++) {
    for (std::size_t j = 0; j < Ic.shape[1]; j++) {
        Ic(i,j) = (1-v1)/(PI*E1)+(1-v2)/(PI*E2) * 1/(sqrt(pow(i,2)+pow(j,2)));
    }
  }

  double beta = 0.84;
  double delta = 0;
  double g_old = 0;
  double min_g = std::min(grid1, grid2);
  std::size_t t = beta*log(min_g);
  double mc = 0.7* pow(min_g, 1/t)-1;
  double eps = pow(min_g, (-3/2));

  double warnings = eps + size + size_p + pressure + delta + g_old + mc;
  warnings += warnings;

  matrix st({2*t, 1});
  initializeMultiplicationArray(st);
  // calculate transfer stylus (8)
  for (std::size_t i = 1; i < 2*t; i++) {
    for (std::size_t j = 0; j < 2*t; j++) {
      if (j!=i) {
        // std::cout << st(i, 0) << " : ";
        double divider = (2*t-2*j)+1;
        double divisor = (2*i-2*j);
        st(i, 0) *= (divider/ divisor );
        // std::cout << st(i, 0) << std::endl;
      }
    }
  }
  //std::cout << t << " is t" << std::endl;
  //printarray(st);

  // condition for coarse grid and pressure calculation
  double sizeN = Ic.shape[0]*Ic.shape[1];
  double sqrt_sizeN = sqrt(sizeN);

  std::size_t newGrid1 = grid1/2 + 2*t - 1;
  std::size_t newGrid2 = grid2/2 + 2*t - 1;

  while(sizeN>sqrt_sizeN) {
    // start iterative loop
    
    matrix coarseDisplacement({newGrid1, newGrid2});
    matrix coarsePressure({newGrid1, newGrid2});

    calc_coarse_pressure(Ip, st, coarsePressure, t);
    newGrid1 = newGrid1/2 + 2*t - 1;
    newGrid2 = newGrid2/2 + 2*t - 1;
    sizeN = newGrid1 * newGrid2;
    std::cout << sizeN << std::endl;
  }

  return 0;
}