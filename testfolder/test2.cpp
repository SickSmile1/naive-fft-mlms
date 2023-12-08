#include <eigen3/Eigen/Core>
#include "Boussinesq.h"

int main() {
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
}