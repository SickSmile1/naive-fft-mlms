#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Sparse>
#include <numeric>
#include "BoussinesqFft.h"
#include "Boussinesq.h"
#include "BoussinesqMlms.h" 

double gMean(sMatrix& gap, std::vector<Eigen::T>& idc);
matrix make_sphere(double radius, int grid, double size, double center);
void ccg(matrix& sub, matrix& topo, double offset, double gradient);

// topography = sphere
// surface = sub

int main() {
  int nx = 16;
  double size = 6., radius = 3.;
  matrix surf = make_sphere(radius,nx,size,radius);
  matrix sub({nx,nx});
  // writeToFile(surf, "surf");
  ccg(sub, surf,.1,.1e-5);
  return 0;
}

void ccg(matrix& sub, matrix& topo, double offset, double gradient) {
  int i = sub.rows();
  int j = sub.cols();

  std::vector<Eigen::T> indices;

  matrix zeros = Eigen::MatrixXd::Zero(i,j);
  matrix g_ij = zeros, t_ij = zeros, h_ij = zeros;
  h_ij(i/2,j/2) = -5;
  h_ij(i/2+1,j/2) = -5;
  matrix p_ij = Eigen::MatrixXd::Zero(2*i-1,2*j-1);
  double delta = 0., G_old = 1.0;
  double eps0 = 5e-5, eps = 1.;
  while (eps > eps0) {
    matrix nd = BoussinesqFFT(6.0, p_ij, topo);
    // matrix nd = BoussinesqMlms(6.0,p_ij, topo, 2);
    matrix gap = nd - h_ij;
    gap = (gap.array() > 0).select(gap, 0);
    sMatrix sgap = gap.sparseView();
    double g_mean = gMean(sgap, indices);
    sgap.coeffs() -= g_mean;
    double G = sgap.coeffs().square().sum();
    std::cout << G << std::endl;
    t_ij = gap + (delta * (G/G_old) * t_ij);

    /*sMatrix temp_t_ij(i,j);
    temp_t_ij.setFromTriplets(indices.begin(), indices.end());
    // std::cout << temp_t_ij << std::endl;
    G_old = G;
    matrix r_ij = BoussinesqFFT(6., p_ij, temp_t_ij);
    // std::cout << r_ij << std::endl;
    double r_mean = gMean(nd, r_ij, indices);

    r_ij.array() = r_ij.array() - r_mean;
    double upper = temp_t_ij.cwiseProduct(gap).sum();
    double lower = temp_t_ij.cwiseProduct(r_ij).sum();
    double roh = upper/lower;
    std::cout << upper << "/"<< lower<<" roh:"<< roh << std::endl;*/
    eps = 0;
  }
}

double gMean(sMatrix& gap, std::vector<Eigen::T>& idc) {
  double g_kl = 0;
  // std::cout << gap << std::endl;
  for (std::size_t k = 0; k < gap.outerSize(); k++) {
    for (Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(gap,k);it;++it) {
      g_kl += it.value();
    }
  }
  std::cout << "NC: "<< gap.nonZeros() << " g_kl: "<<g_kl << std::endl;
  return std::pow(gap.nonZeros(),-1.)*g_kl;
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