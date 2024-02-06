#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Sparse>
#include <numeric>
#include "BoussinesqFft.h"
#include "Boussinesq.h"
#include "BoussinesqMlms.h" 
#include "Ccg.h"


int main() {
  int nx = 16;
  double size = 6., radius = 3., pixel = size/nx;
  matrix surf = make_sphere(radius,nx,size,radius);
  matrix sub({nx,nx});
  // writeToFile(surf, "surf");
  ccg(sub, surf,15.,.1e-5, pixel);
  return 0;
}

void ccg(matrix& sub, matrix& topo, double offset, double gradient, double pixel) {
  int i = sub.rows();
  int j = sub.cols();
  double ext_force = offset;
  std::cout << "Pixel: " << pixel << std::endl;
  Eigen::Index maxRow, maxCol;
  double max = topo.minCoeff(&maxRow, &maxCol);

  matrix zeros = Eigen::MatrixXd::Zero(i,j);
  matrix g_ij = zeros, t_ij = zeros, h_ij = zeros, gap = zeros;
  // p_ij(maxRow, maxCol) = -max;
  sMatrix old_p_ij, sr_ij;
  matrix p_ij = sub;
  h_ij.array() = topo.array();
  // std::cout << h_ij << std::endl;

  double P0 = 15.;// topo.sum()*pixel*pixel;
  double delta = 0., G_old = 1.0;
  double eps0 = 5e-5, eps = 3.;
  while (eps > eps0) {
    /* 21 */
    matrix nd = BoussinesqFFT(6.0, p_ij);

    // matrix nd = BoussinesqMlms(6.0,p_ij, topo, 2);
    // std::cout << nd << std::endl;
    // gap = nd - h_ij;
    gap.array() = -nd.array() - ext_force;
    std::cout << "pij: \n" << nd << std::endl;
    gap = (p_ij.array() != 0).select(gap, 0); // <- this is wrong... we need p_ij>0 or the contact area
    std::cout << "gap: \n" << gap << std::endl;
    sMatrix sgap = gap.sparseView();
    double g_mean = gMean(sgap);
    sgap.coeffs() -= g_mean;
    /* 22 */
    double G = sgap.coeffs().square().sum();
    std::cout << "G: " << G << std::endl;
    /* 23, set t_ij \in Ic by multiplication with sparse sgap */
    t_ij = sgap + (delta * (G/G_old) * t_ij);
    G_old = G;

    std::cout << "t_ij: \n" << t_ij << std::endl;
    /* 24 */
    matrix r_ij = BoussinesqFFT(6., t_ij);
    std::cout << "rij: \n" << r_ij << std::endl;

    /* 25 */
    r_ij = (r_ij.array() > 0).select(r_ij, 0);
    sr_ij = r_ij.sparseView();
    double r_mean = gMean(sr_ij);
    sr_ij.coeffs() -= r_mean;

    /* 26 */
    double upper = sgap.cwiseProduct(t_ij).sum();
    double lower = sr_ij.cwiseProduct(t_ij).sum();
    double roh = upper/lower;
    // std::cout << upper << "/"<< lower<<" roh:"<< roh << std::endl;

    /* 27 */
    old_p_ij = p_ij.sparseView();

    /* 28 */
    p_ij -= (roh * t_ij);
    p_ij = (p_ij.array() > 0).select(p_ij, 0);

    
    /* 29 Iol returns empty(true/false), triplet list holds the index/val of sgap*/
    typedef Eigen::Triplet<double> T;
    std::vector<T> tripletList;
    delta = Iol(p_ij, sgap, tripletList) ? 1 : 0;
    std::cout << "delta: " << delta << std::endl;

    /* 30 */
    for (auto it: tripletList) {
      p_ij(it.row(),it.col()) -= roh*it.value();
      // std::cout << "-roh"<< std::endl;
    }

    double P1 = p_ij.sum(); // pixel*pixel;
    if (P1 == 0) {
      p_ij.array() = p_ij.array() - (ext_force/ (i*j) );
      std::cout << "want to add some value \n" << p_ij << std::endl;
    } else {
      p_ij.array() = (P1/P0)*p_ij.array();
    }
    ext_force = p_ij.sum() / p_ij.outerSize();
    // std::cout << p_ij << std::endl;
    double eps_test = pixel*pixel *std::pow(P0,-1) * (p_ij-old_p_ij).cwiseAbs().sum();
    // std::cout << "eps: " <<  eps_test << " iter:" << eps << std::endl;
    eps -= 1;
  }
}

double gMean(sMatrix& gap) {
  double g_kl = 0;
  // std::cout << gap << std::endl;
  for (std::size_t k = 0; k < gap.outerSize(); k++) {
    for (Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(gap,k);it;++it) {
      g_kl += it.value();
    }
  }
  // std::cout << "NC: "<< gap.nonZeros() << " g_kl: "<<g_kl << std::endl;
  return std::pow(gap.nonZeros(),-1.)*g_kl;
}

bool Iol(const matrix& p_ij, const sMatrix& gap, std::vector<Eigen::T> idc) {
  bool empty = true;
  for (std::size_t k = 0; k < gap.outerSize(); k++) {
    for (Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(gap,k);it;++it) {
      if (it.value() < 0) {
        if (p_ij(it.row(),it.col())==0)
          empty = false;
          idc.push_back(Eigen::T(it.row(),it.col(),it.value()));
      }
    }
  }
  return empty;
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