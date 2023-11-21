#include "BoussinesqFft.h"
#include "Ccg.h"
#include <eigen3/Eigen/Core>
#include "eigen3/Eigen/SparseCore"


int main() {
  matrix s({10,10});
  s.setOnes();
  matrix p({10,10});
  p.setZero();
  s = s*0.005;
  ccg(s, p, .005, .002);
  return 0;
}

int gMean(const matrix& mat, const matrix& gap) {
  int Nc = 0;
  double g_kl = 0;
  for (std::size_t k = 0; k < mat.rows(); k++) {
    for (std::size_t l = 0; l < mat.cols(); l++) {
      bool val = (mat(k,l)>0);
      double ret = val*mat(k,l);
      g_kl += ret * gap(k,l);
    }
  }
  return Nc*g_kl;
}

void ccg(matrix sub, matrix topo, double offset, double gradient) {
  int i = sub.rows();
  int j = sub.cols();
  matrix zeros = Eigen::MatrixXd::Zero(i,j);
  matrix g_ij = zeros;
  matrix t_ij = zeros;
  matrix p_ij = zeros;
  double delta = 0;
  double G_old = 1.0;
  double eps0 = 5e-5, eps = 1.;
  
  while (eps > eps0) {
    matrix nd = BoussinesqFFT(1.0, sub, p_ij);
    matrix gap = - nd - topo;
    auto g_mean = gMean(nd, gap);
    eps = 0;
  }
}
