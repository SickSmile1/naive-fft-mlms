#ifndef MATRIX_HH
#define MATRIX_HH

#include <vector>
#include <array>
#include <cassert>
#include <complex>

template <typename T>
struct matrix_template {
  using shape_t = std::array<std::size_t, 2>;

  matrix_template(shape_t shape): shape(shape) {
    data.resize(shape[0] * shape[1]);
  }

  T& operator()(std::size_t i, std::size_t j) {
    assert(i < shape[0] and j < shape[1]);
    return data[i * shape[1] + j];
  }

  const T& operator()(std::size_t i, std::size_t j) const {
    assert(i < shape[0] and j < shape[1]);
    return data[i * shape[1] + j];
  }

  decltype(auto) begin() {
    return data.begin();
  }

  decltype(auto) end() {
    return data.end();
  }
  
  shape_t shape;
  std::vector<T> data;
};

using complex = std::complex<double>;
using matrix = matrix_template<double>;
using matrixc = matrix_template<complex>;
using vec2d = std::array<double, 2>;

template <typename T>
std::ostream& operator<<(std::ostream& stream, const matrix_template<T>& A) {
  stream << std::fixed;
  for (std::size_t i = 0; i < A.shape[0]; ++i) {
    for (std::size_t j = 0; j < A.shape[1]; ++j) {
      stream << A(i, j) << " ";
    }

    stream << '\n';
  }

  stream << std::defaultfloat;
  return stream;
}

#endif // MATRIX_HH
