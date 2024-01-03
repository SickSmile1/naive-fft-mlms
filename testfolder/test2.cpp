#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Sparse>
#include <numeric>
#include "BoussinesqFft.h"
#include "Boussinesq.h"
#include "BoussinesqMlms.h"
#include <random>

int main() {
  int N = 10;
  double E =.01, feastol = 1e-5, C = 0, radius = 0, delta = 1.;

  int n = N*N;
  matrix m({N,N}), Q({N, N}), r({N, N}), w({N,N}),f({N,N});
  // Eigen::ArrayXd ub = Eigen::ArrayXd::Random(n);
  matrix ub = Eigen::MatrixXd::Random(N,N);
  matrix kept = Eigen::MatrixXd::Ones(N,N);
  matrix ones = Eigen::MatrixXd::Ones(N,N);

  Eigen::SparseLU<Eigen::SparseMatrix<double> > solA;

  m.row(0) = Eigen::VectorXd::LinSpaced(N,1, N);
  m = m.row(0).replicate(N,1); // 

  radius = delta/2;
  C = E*M_PI*radius;

  Q.setIdentity();
  Q.array() = Q.array() * C;

  // auto xv(Eigen::Map<Eigen::VectorXd>(m.data(), n) );
  // m.transposeInPlace();
  // auto yv = m.transpose().data();
  auto xv = m.reshaped();
  auto yv = m.transpose().reshaped();
  // std::cout << xv << std::endl;
  for(int i = 0; i < N; ++i) {
    for (int j = i; j< N; ++j) {
      if (i!=j) {
        r(i,j) = std::sqrt( std::pow( (xv[j] - xv[i]) ,2) +
            std::pow( (yv[j]-yv[i]),2) ) ;
        Q(i,j) = C*std::asin( radius/r(i,j) );
        Q(j,i) = Q(i,j);
      }
    }
  }

  Eigen::SparseMatrix<double, Eigen::ColMajor> sQ({N, N}), sub({N, N});
  // Eigen::VectorXd p = Eigen::VectorXd::Zero(n), minw;
  Eigen::VectorXd minw;
  bool keepgoing = true;
  while(keepgoing) {
    keepgoing = false;

    Q = (kept.array() > 0).select(Q,0);
    ub = (kept.array() > 0).select(ub,0);

    sQ = Q.sparseView();
    sQ.makeCompressed();
    solA.analyzePattern(sQ);
    solA.factorize(sQ);

    // auto arr_ub = ub.reshaped();
    Eigen::VectorXd p = solA.solve(ub);

    f = (p.array() < -feastol).select(ones,0);
    keepgoing = (p.array() < -feastol).any();
    kept = (f.array() > 0).select(0, ones);
    w = (Q*p.adjoint());
    w -= ub;
    minw << w.minCoeff();
    ub = sub.toDense();
  }

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