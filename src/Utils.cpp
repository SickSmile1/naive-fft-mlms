#include <numeric>
#include "Boussinesq.h"

matrix make_sphere(double radius, int grid, double size, double center);

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

int main() {
  matrix ret = make_sphere(3., 33., 6.,3. );
  std::cout << ret;
  return 0;
}