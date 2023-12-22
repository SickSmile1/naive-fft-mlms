#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Sparse>
#include <numeric>
#include "BoussinesqFft.h"
#include "Boussinesq.h"
#include "BoussinesqMlms.h" 

int main() {
  int N = 10;
  double E =.01, feastol = 1e-5, C = 0, radius = 0, delta = 1.;

  int n = N*N;
  matrix m({N,N}), r({N,N}), Q({N*N, N*N});
  m.row(0) = Eigen::VectorXd::LinSpaced(N,1, N);
  m = m.row(0).replicate(N,1); // 

  radius = delta/2;
  C = E*M_PI*radius;

  Q.setIdentity();
  Q.array() = Q.array() * C;

  auto xv = m.array();
  auto yv = m.transpose().array();

  for(int i = 0; i < n; ++i) {
    for (int j = i; j< n; ++j) {
      if (i!=j) {
        r(i,j) = std::sqrt( std::pow( (xv(j) - xv(i)) ,2) +
            std::pow( (yv(j)-yv(i)),2) ) ;
        Q(i,j) = C*std::asin( radius/r(i,j) );
        Q(j,i) = Q(i,j);
      }
    }
  }
  // std::cout << m << std::endl;
  std::cout << Q << std::endl;
  std::cout << r << std::endl;
  // std::cout <<  << std::endl;
}

// topography = sphere
// surface = sub



/*
  for(int u = 4; u < 9; u++) {
    int n = (2*u)+1;
    int oShape = (n)/2+1;
    int shape = (n/2)-1;
    std::cout << oShape << " : " << shape << ":" << n << std::endl;
    matrix s({n,n});
    s.setZero();
    for (int i = 0; i < oShape; i++) {
      for (int j = 0; j < oShape; j++) {
        s(i,j) = i+j;
      }
    }
    s.block(0,oShape,oShape,shape+1) = s.block(0,1,oShape,shape+1).rowwise().reverse();
    // std::cout << s << std::endl;
    s.block(oShape,0,shape+1,oShape+shape+1) = s.block(1,0,shape+1,oShape+shape+1).colwise().reverse();
    std::cout << s << "\n" << std::endl;
  }
  return 0;
*/