#include <eigen3/Eigen/core>
#include "Boussinesq.h"
#include "eigen3/Eigen/sparseCore"

int main() {
  matrix ss({4, 4});
  matrix s({3, 3});
  s(1,1) = 1;
  s(1,2) = 2;
  auto r = s.sparseView();
  std::cout << typeid(r).name() << std::endl;
  return 0;
}
