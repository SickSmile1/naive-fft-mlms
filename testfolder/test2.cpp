#include <eigen3/Eigen/Core>
#include <numeric>
#include "BoussinesqFft.h"
#include "Boussinesq.h"

int gMean(const matrix& mat, const matrix& gap);
matrix make_sphere(double radius, int grid, double size, double center);
void ccg(matrix& sub, matrix& topo, double offset, double gradient);

// topography = sphere
// surface = sub

int main() {
  int nx = 16;
  double size = 6., radius = 3.;
  matrix surf = make_sphere(radius,nx,size,radius);
  matrix sub({nx,nx});
  ccg(sub, surf,.1,.1e-5);
  return 0;
}

void ccg(matrix& sub, matrix& topo, double offset, double gradient) {
  int i = sub.rows();
  int j = sub.cols();
  matrix zeros = Eigen::MatrixXd::Zero(i,j);
  matrix g_ij = zeros;
  matrix t_ij = zeros;
  matrix p_ij = Eigen::MatrixXd::Zero(2*i-1,2*j-1);
  double delta = 0.;
  double G_old = 1.0;
  double eps0 = 5e-5, eps = 1.;
  // std::cout << p_ij.rows() << " : " << p_ij.cols() << std::endl;
  // std::cout << sub.rows() << " : " << sub.cols() << std::endl;
  while (eps > eps0) {
    matrix nd = BoussinesqFFT(6.0, p_ij, topo);
    std::cout << nd << std::endl;
    matrix gap = - nd - topo;
    double g_mean = gMean(nd, gap);
    std::cout << g_mean << std::endl;
    eps = 0;
  }
}

int gMean(const matrix& mat, const matrix& gap) {
  int Nc = 0;
  double g_kl = 0;
  for (std::size_t k = 0; k < mat.rows(); k++) {
    for (std::size_t l = 0; l < mat.cols(); l++) {
      bool val = (mat(k,l)>0);
      Nc += val;
      double ret = val*mat(k,l);
      g_kl += ret * gap(k,l);
    }
  }
  return Nc*g_kl;
}

matrix make_sphere(double radius, int grid, double size, double center) {
  matrix rx2({grid,grid});
  rx2.row(0) = Eigen::VectorXd::LinSpaced(grid,0, grid-1)*size/grid;
  rx2.row(0) = (rx2.row(0).array() - center).pow(2);
  rx2 = rx2.row(0).replicate(grid,1); // ;
  Eigen::VectorXd temp = rx2.row(0);
  rx2.colwise() += temp;
  rx2.array() = - rx2.array() / (2*radius);
  return rx2;
}

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