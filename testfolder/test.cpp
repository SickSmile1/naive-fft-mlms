#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Sparse>
#include "Boussinesq.h"
#include "eigen3/Eigen/SparseCore"
#include <benchmark/benchmark.h>


int main() {
  typedef Eigen::Triplet<double> T;
  std::vector<T> tripletList;
  
  for(int k = 3; k < 6; k++) {
    tripletList.push_back(T(k, k,5));
  }
  Eigen::SparseMatrix<double,Eigen::RowMajor> t({10,10});
  t.setFromTriplets(tripletList.begin(), tripletList.end());
  // std::cout << t << std::endl;
  matrix s({10,10});
  s.setIdentity();
  matrix q = (s.array() > 0.).select(s,t.toDense());
  sMatrix sq = q.sparseView();
  std::cout << sq << std::endl;
  //matrix r = s.array().reverse();
  //Eigen::SparseMatrix<double,Eigen::RowMajor> q = r.sparseView();
  //q = q.cwiseProduct(t);
  /*for (std::size_t i = 0; i < q.rows(); i++) {
    for(std::size_t j = 0; j < q.cols();j++) {
      std::cout << q(i,j) << "\t";
    }
    std::cout << std::endl;
  }*/
  // std::cout << q << std::endl;

  return 0;
}



/*void sparseBench(size_t n) {
  matrix s = Eigen::MatrixXd::Zero(n,n);
  s(n/2,1) = 1;
  s(n-2,n-1) = 2;
  Eigen::SparseMatrix<double> sparseView = s.sparseView();
  for (int i = 0; i < sparseView.outerSize(); ++i) {
    for (Eigen::SparseMatrix<double>::InnerIterator it(sparseView, i); it; ++it) {
        int row = it.row();
        int col = it.col();
        double value = it.value();
        s(row,col) = 0;
    }
  }
  std::cout << "outerSize: " << sparseView.outerSize() << " nonzero: "<< sparseView.nonZeros() << std::endl;
}

void normalIter(size_t n) {
  matrix s = Eigen::MatrixXd::Zero(n,n);
  s(n/2,1) = 1;
  s(n-2,n-1) = 2;
  for (size_t i = 0; i < n; i++) {
    for(size_t j = 0; j < n; j++) {
      if(s(i,j)>0) s(i,j) = 0;
    }
  }
}

static void sparseBench_t(benchmark::State &state) {
  for (auto _ : state) {
    sparseBench(state.range(0));
  }
}

static void sparseBench_t2(benchmark::State &state) {
  for (auto _ : state) {
    normalIter(state.range(0));
  }
}

BENCHMARK(sparseBench_t)->RangeMultiplier(2)->Range(8,8<<10)->Unit(benchmark::kMillisecond);
BENCHMARK(sparseBench_t2)->RangeMultiplier(2)->Range(8,8<<10)->Unit(benchmark::kMillisecond);

BENCHMARK_MAIN();*/
