#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Sparse>
#include <numeric>
#include "BoussinesqFft.h"
#include "Boussinesq.h"
#include "BoussinesqMlms.h" 

double gMean(sMatrix& gap);
matrix make_sphere(double radius, int grid, double size, double center);
void ccg(matrix& sub, matrix& topo, double offset, double gradient, double pixel);
bool Iol(const matrix& p_ij, const sMatrix& gap, std::vector<Eigen::T> idc);
void ccg(matrix sub,matrix topo, double offset, double gradient);
