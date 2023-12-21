#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Sparse>
#include <numeric>
#include "BoussinesqFft.h"
#include "Boussinesq.h"
#include "BoussinesqMlms.h" 

double gMean(sMatrix& gap);
matrix make_sphere(double radius, int grid, double size, double center);
void ccg(matrix& sub, matrix& topo, double offset, double gradient);
bool Iol(const matrix& p_ij, const sMatrix& gap, std::vector<Eigen::T> idc);

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


  matrix zeros = Eigen::MatrixXd::Zero(i,j);
  matrix g_ij = zeros, t_ij = zeros, h_ij = zeros, gap = zeros;
  sMatrix old_p_ij, sr_ij;
  h_ij(i/2,j/2) = -5;
  h_ij(i/2+1,j/2) = -5;
  matrix p_ij = Eigen::MatrixXd::Zero(2*i-1,2*j-1);
  double delta = 0., G_old = 1.0;
  double eps0 = 5e-5, eps = 5.;
  while (eps > eps0) {
    /* 21 */
    matrix nd = BoussinesqFFT(6.0, p_ij, topo);
    // matrix nd = BoussinesqMlms(6.0,p_ij, topo, 2);
    if(eps == 5) gap = nd - h_ij;
    else gap = nd - p_ij.block(0,0,topo.rows(), topo.cols());
    gap = (gap.array() > 0).select(gap, 0);
    sMatrix sgap = gap.sparseView();
    double g_mean = gMean(sgap);
    sgap.coeffs() -= g_mean;
    /* 22 */
    double G = sgap.coeffs().square().sum();

    /* 23, set t_ij \in Ic by multiplication with sparse sgap */
    t_ij = sgap + (delta * (G/G_old) * t_ij);
    G_old = G;

    /* 24 */
    matrix r_ij = BoussinesqFFT(6., p_ij, t_ij);

    /* 25 */
    r_ij = (r_ij.array() > 0).select(r_ij, 0);
    sr_ij = r_ij.sparseView();
    double r_mean = gMean(sr_ij);
    sr_ij.coeffs() -= r_mean;

    /* 26 */
    double upper = sgap.cwiseProduct(t_ij).sum();
    double lower = sr_ij.cwiseProduct(t_ij).sum();
    double roh = upper/lower;
    std::cout << upper << "/"<< lower<<" roh:"<< roh << std::endl;

    /* 27 */
    old_p_ij = p_ij.sparseView();

    /* 28 */
    p_ij.block(0,0,t_ij.rows(), t_ij.cols()) -= (roh * t_ij);
    p_ij = (p_ij.array() > 0).select(p_ij, 0);

    std::cout << p_ij.block(0,0,topo.rows(), topo.cols()) << std::endl;
    
    /* 29 Iol returns empty(true/false), triplet list holds the index/val of sgap*/
    typedef Eigen::Triplet<double> T;
    std::vector<T> tripletList;
    delta = Iol(p_ij, sgap, tripletList) ? 1 : 0;
    std::cout << "delta: " << delta << std::endl;

    /* 30 */
    for (auto it: tripletList) {
      p_ij(it.row(),it.col()) -= roh*it.value();
      std::cout << "-roh"<< std::endl;
    }
    
    eps -= 1;
    topo = p_ij.block(0,0,topo.rows(),topo.cols());
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
  std::cout << "NC: "<< gap.nonZeros() << " g_kl: "<<g_kl << std::endl;
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