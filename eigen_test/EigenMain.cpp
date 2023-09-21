#include <iostream>
#include <Eigen/Core>



int main() {
    Eigen::MatrixXd m({4,4}); // no operation
    std::cout << "The matrix m is of size "
              << m.rows() << "x" << m.cols() << std::endl;
}