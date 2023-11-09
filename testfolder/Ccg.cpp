#include "BoussinesqFft.h"
#include "Ccg.h"
#include "eigen3/Eigen/Core"
#include "eigen3/Eigen/Sparse"

void ccg(matrix sub, matrix topo, double offset, double gradient) {
  matrix nodDeflections = BoussinesqFFT(2.0, sub, topo);
  int i = nodDeflections.rows();
  int j = nodDeflections.cols();
  matrix g_ij({i,j}); g_ij = Eigen::MatrixXd::Zero(i,j);
  matrix t_ij({i,j}); t_ij= Eigen::MatrixXd::Zero(i,j);
  double delta = 0;
  double G_old = 1.0;
  matrix gap({i,j});
  // gap = -nodDeflections - (topo + nodDeflections);
  // double g_mean = nodDeflections.nonZeros() * gap.sum();
}

int main() {
  matrix s({33,33});
  s.array() = 1;// = Eigen::MatrixXd::Ones(32,32);
  matrix p({33,33});
  p.array() = 1;//) = Eigen::MatrixXd::Ones(32,32);
  p.leftCols(8) = Eigen::MatrixXd::Zero(33,8);
  p.rightCols(8) = Eigen::MatrixXd::Zero(33,8);
  p.bottomRows(8) = Eigen::MatrixXd::Zero(8,33);
  p.topRows(8) = Eigen::MatrixXd::Zero(8,33);
  p = p*(-0.05);
  // std::cout << -p << std::endl;
  s = s*0.005;
  // std::cout << s << std::endl;
  ccg(s, p, .005, .002);
  return 1;
}
