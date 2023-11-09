#include "BoussinesqFft.h"
#include "Ccg.h"
#include <eigen3/Eigen/Core>
#include "eigen3/Eigen/SparseCore"


int main() {
  matrix s({10,10,1.0});
  s = s*0.005;
  ccg(s, -s, .005, .002, .003 );
  return 0;
}

void ccg(matrix sub, matrix topo, double offset, double gradient, double gap) {
  matrix nd = BoussinesqFFT(1.0, sub, topo);
  Eigen::SparseMatrixBase<double> nodDeflections = nd.sparseView();
  int i = nodDeflections.rows();
  int j = nodDeflections.cols();
  /*Eigen::SparseMatrix<double> g_ij = Eigen::MatrixXd::Zero(i,j);
  Eigen::SparseMatrix<double> t_ij = Eigen::MatrixXd::Zero(i,j);
  double delta = 0;
  double G_old = 1.0;
  Eigen::SparseMatrix<double> gap = -nodDeflections - (topo + nodDeflections);
  auto g_mean = nodDeflections.nonZeros() * gap.sum();*/
}
