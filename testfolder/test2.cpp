#include <eigen3/Eigen/Core>
#include "Boussinesq.h"

int main() {
  for(int n = 9; n <15; n++) {
    int oShape = (n-1)/2;
    int shape = (n/2)-1;
    std::cout << oShape << " : " << shape << std::endl;
    matrix s({shape+oShape,shape+oShape});
    s.setZero();
    for (int i = 0; i < oShape; i++) {
      for (int j = 0; j < oShape; j++) {
        s(i,j) = i+j;
      }
    }
    s.block(0,oShape,oShape,shape) = s.block(0,0,oShape,shape).rowwise().reverse();
    // std::cout << s << std::endl;
    s.block(oShape,0,shape,shape+oShape) = s.block(0,0,shape,shape+oShape).reverse();
    std::cout << s << "\n" << std::endl;
  }
  return 0;
}