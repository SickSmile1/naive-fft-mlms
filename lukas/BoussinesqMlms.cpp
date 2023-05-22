#include "BoussinesqMlms.hh"
#include "Boussinesq.hh"

#include <algorithm>
#include <numeric>
#include <cmath>
#include <iostream>

std::size_t optimal_levels(std::array<std::size_t, 2> shape) {
  double N = shape[0] * shape[1];
  return std::log2(N) / 2 - 1;
}

grid_stack::grid_stack(const matrix& fine, vec2d pixel, int t, std::size_t levels) : t(t) {
  stack.push_back(fine);
  pixel_sizes.push_back(pixel);
  stencil = make_stencil(t);

  auto one = make_stencil(2);
  auto two = make_stencil(4);
  auto three = make_stencil(6);


  for (int i = 0; i < one.size(); i++) {
      std::cout << one[i] << std::endl;
  }
  
  std::cout << "\n\n" << std::endl;
  for (int i = 0; i < two.size(); i++) {
      std::cout << two[i] << std::endl;
  }
  std::cout << "\n\n" << std::endl;
   for (int i = 0; i < three.size(); i++) {
      std::cout << three[i] << std::endl;
  }
  std::cout << "\n\n" << std::endl;
  for (std::size_t q = 1; q < levels; ++q) {
    add_level();
  }

  assert(stack.size() == levels and pixel_sizes.size() == levels);
}

void grid_stack::add_level() {
  const auto& pix = pixel_sizes.back();
  const auto& shape = stack.back().shape;

  pixel_sizes.push_back({2 * pix[0], 2 * pix[1]});
  stack.push_back(matrix({shape[0] / 2 + 2 * t - 1, shape[1] / 2 + 2 * t - 1}));
}

void grid_stack::coarsen() {
  for (std::size_t q = 1; q < stack.size(); ++q)
    interpolate(stack[q-1], stack[q], stencil);
}

std::vector<double> make_stencil(int t) {
  std::vector<double> s(2 * t);
  std::iota(s.begin(), s.end(), 1);
  std::transform(s.begin(), s.end(), s.begin(), [t](double i) {
    double res = 1;
    for (int k = 1; k <= 2 * t; ++k)
      res *= (i == k) ? 1 : (2 * (t - k) + 1) / (2 * (i - k));
    return res;
  });
  return s;
}

std::size_t correction_size(int t, std::array<std::size_t, 2> shape) {
  double Nmin = std::min(shape[0], shape[1]);
  return static_cast<std::size_t>(std::max(0.7 * t * std::pow(Nmin, 1. / t) - 1,
                                           2. * t));
}

double interpolate(const matrix& fine, std::vector<double> stencil, int m, int n) {
  int t = stencil.size() / 2;
  int i = 2 * (m - t + 1), j = 2 * (n - t + 1);
  int Nx = fine.shape[0], Ny = fine.shape[1];

  auto padded = [&fine, Nx, Ny] (int i, int j) {
    if (i >= 0 and i < Nx and j >= 0 and j < Ny)
      return fine(i, j);
    return 0.;
  };

  double res = padded(i, j);
  for (int k = 1; k <= 2 * t; ++k)
    res += stencil[k-1] * padded(i + 2*(k - t) - 1, j);

  for (int k = 1; k <= 2 * t; ++k)
    res += stencil[k-1] * padded(i, j + 2*(k - t) - 1);

  for (int k = 1; k <= 2 * t; ++k)
    for (int l = 1; l <= 2 * t; ++l)
      res+= stencil[k-1] * stencil[l-1] * padded(i + 2*(k - t) - 1, j + 2*(l - t) -1);

  return res;
}

void interpolate(const matrix& fine, matrix& coarse, std::vector<double> stencil) {
  for (std::size_t m = 0; m < coarse.shape[0]; ++m)
    for (std::size_t n = 0; n < coarse.shape[1]; ++n)
      coarse(m, n) = interpolate(fine, stencil, m, n);
}

void compute_coarse_displacements(const grid_stack& pressures, grid_stack& displacements) {
  auto&& cpressure = pressures.stack.back();
  auto&& cdisplacement = displacements.stack.back();

  naive(cpressure, cdisplacement,
        pressures.pixel_sizes.back(), pressures.pixel_sizes.front());
}

double coarse_correction(int i, int j,
                         const std::vector<double>& stencil,
                         const vec2d& metric, const vec2d& pixel) {
  if (i % 2 == 0 and j % 2 == 0)
    return 0.;

  int t = stencil.size() / 2;
  double K = boussinesq(i, j, metric, pixel);
  double K_interp = 0;

  if (j % 2 == 0)
    for (int k = 1; k <= 2 * t; ++k)
      K_interp += stencil[k-1]
        * boussinesq((i - 2 * (k - t) + 1), j, metric, pixel);
  else if (i % 2 == 0)
    for (int k = 1; k <= 2 * t; ++k)
      K_interp += stencil[k-1]
        * boussinesq(i, (j - 2 * (k - t) + 1), metric, pixel);
  else
    for (int k = 1; k <= 2 * t; ++k)
      for (int l = 1; l <= 2 * t; ++l)
        K_interp += stencil[k-1] * stencil[l-1]
          * boussinesq((i - 2 * (k - t) + 1),
                       (j - 2 * (l - t) + 1),
                       metric, pixel);
    

  return K - K_interp;
}

matrix coarse_correction(std::size_t mc, vec2d metric, vec2d pixel,
                         const std::vector<double>& stencil) {
  matrix C({2 * mc + 1, 2 * mc + 1});

  for (std::size_t a = 0; a < C.shape[0]; ++a) {
    int i = static_cast<int>(a) - static_cast<int>(mc);
    for (std::size_t b = 0; b < C.shape[1]; ++b) {
      int j = static_cast<int>(b) - static_cast<int>(mc);
      C(a, b) = coarse_correction(i, j, stencil, metric, pixel);
    }
  }

  return C;
}

std::vector<matrix> coarse_corrections(const grid_stack& pressure, std::size_t mc) {
  std::vector<matrix> corrections;
  vec2d pixel = pressure.pixel_sizes.front();

  for (auto&& metric : pressure.pixel_sizes)
    corrections.push_back(coarse_correction(mc, metric, pixel, pressure.stencil));
  return corrections;
}


double coarse_displacement_correction(int i, int j,
                                      int t,
                                      const matrix& correction,
                                      const matrix& pfine) {
  double res = 0;
  int mc = (correction.shape[0] - 1) / 2;

  auto padded = [&pfine] (int i, int j) {
    if (i >= 0 and i < static_cast<int>(pfine.shape[0])
        and j >= 0 and j < static_cast<int>(pfine.shape[1]))
      return pfine(i, j);
    return 0.;
  };

  for (int k = -mc; k <= mc; ++k)
    for (int l = -mc; l <= mc; ++l)
      res += correction(k + mc, l + mc) * padded(2 * (i - t + 1) - k, 2 * (j - t + 1) - l);
  return res;
}

void apply_coarse_correction(const matrix& correction,
                             const matrix& pfine,
                             matrix& ucoarse) {
  int t = (ucoarse.shape[0] - pfine.shape[0] / 2 + 1) / 2;
  for (std::size_t i = 0; i < ucoarse.shape[0]; ++i)
    for (std::size_t j = 0; j < ucoarse.shape[1]; ++j)
      ucoarse(i, j) += coarse_displacement_correction(i, j, t, correction, pfine);
}

double refine(int i, int j,
              const matrix& coarse, const std::vector<double>& stencil) {
  int t = stencil.size() / 2;
  int m = (i + 1) / 2 + t - 1;
  int n = (j + 1) / 2 + t - 1;

  if (i % 2 == 0 and j % 2 == 0)
    return coarse(m, n);

  double res = 0;
  if (j % 2 == 0)
    for (int k = 1; k <= 2 * t; ++k)
      res += stencil[k-1] * coarse(m + k - t - 1, n);

  else if (i % 2 == 0)
    for (int k = 1; k <= 2 * t; ++k)
      res += stencil[k-1] * coarse(m, n + k - t - 1);

  else
    for (int k = 1; k <= 2 * t; ++k)
      for (int l = 1; l <= 2 * t; ++l)
        res += stencil[k-1] * stencil[l-1]
          * coarse(m + k - t - 1, n + l - t - 1);

  return res;
}

void refine(const matrix& coarse, matrix& fine, std::vector<double> stencil) {
  for (std::size_t i = 0; i < fine.shape[0]; ++i)
    for (std::size_t j = 0; j < fine.shape[1]; ++j)
      fine(i, j) = refine(i, j, coarse, stencil);
}


double K_interp2(int i, int j, int t, const std::vector<double>& stencil,
                 const vec2d& metric, const vec2d& pixel) {
  double res = 0;
  for (int k = 1; k <= 2 * t; ++k)
    res += stencil[k-1]
      * boussinesq((i + 2 * (k - t) - 1), j, metric, pixel);
  return res;
}

double K_interp3(int i, int j, int t, const std::vector<double>& stencil,
                 const vec2d& metric, const vec2d& pixel) {
  double res = 0;
  for (int k = 1; k <= 2 * t; ++k)
    res += stencil[k-1]
      * boussinesq(i, (j + 2 * (k - t) - 1), metric, pixel);
  return res;
}

double K_interp4(int i, int j, int t, const std::vector<double>& stencil,
                 const vec2d& metric, const vec2d& pixel) {
  double res = 0;
  for (int k = 1; k <= 2 * t; ++k)
    for (int l = 1; l <= 2 * t; ++l)
      res += stencil[k-1] * stencil[l-1]
        * boussinesq((i + 2 * (k - t) - 1),
                     (j + 2 * (l - t) - 1),
                     metric, pixel);
  return res;
}

template <typename Func>
double fine_correction(int i, int j,
                       const std::vector<double>& stencil,
                       const vec2d& metric, const vec2d& pixel, Func interp) {
  double K = boussinesq(i, j, metric, pixel);
  int t = stencil.size() / 2;
  double K_interp = interp(i, j, t, stencil, metric, pixel);
  return K - K_interp;
}

template <typename Func>
matrix fine_correction(std::size_t mc, const vec2d& metric, const vec2d& pixel,
                       const std::vector<double>& stencil, Func func) {
  matrix C({2 * mc + 1, 2 * mc + 1});

  for (std::size_t a = 0; a < C.shape[0]; ++a) {
    int i = static_cast<int>(a) - static_cast<int>(mc);
    for (std::size_t b = 0; b < C.shape[1]; ++b) {
      int j = static_cast<int>(b) - static_cast<int>(mc);
      C(a, b) = fine_correction(i, j, stencil, metric, pixel, func);
    }
  }

  return C;
}

std::vector<std::vector<matrix>> fine_corrections(const grid_stack& stack, std::size_t mc) {
  std::vector<std::vector<matrix>> corrections;
  vec2d pixel = stack.pixel_sizes.front();

  for (auto&& metric : stack.pixel_sizes) {
    std::vector<matrix> level_corrections;
    level_corrections.push_back(fine_correction(mc, metric, pixel,
                                                stack.stencil, K_interp2));
    level_corrections.push_back(fine_correction(mc, metric, pixel,
                                                stack.stencil, K_interp3));
    level_corrections.push_back(fine_correction(mc, metric, pixel,
                                                stack.stencil, K_interp4));
    corrections.push_back(std::move(level_corrections));
  }

  return corrections;
}

double fine_displacement_correction(int i, int j,
                                    const matrix& correction,
                                    const matrix& pfine) {
  double res = 0;
  int mc = (correction.shape[0] - 1) / 2;

  auto padded = [&pfine] (int i, int j) {
    if (i >= 0 and i < static_cast<int>(pfine.shape[0])
        and j >= 0 and j < static_cast<int>(pfine.shape[1]))
      return pfine(i, j);
    return 0.;
  };

  for (int k = -mc; k <= mc; ++k)
    for (int l = -mc; l <= mc; ++l)
      res += correction(k + mc, l + mc) * padded(i - k, j - l);
  return res;
}

void apply_fine_correction(const std::vector<matrix>& corrections,
                           const matrix& pfine,
                           matrix& ufine) {
  for (std::size_t i = 0; i < ufine.shape[0]; ++i) {
    for (std::size_t j = 0; j < ufine.shape[1]; ++j) {
      if (i % 2 == 0 and j % 2 == 0)
        continue;

      std::size_t idx;
      if (j % 2 == 0)
        idx = 0;
      else if (i % 2 == 0)
        idx = 1;
      else
        idx = 2;

      ufine(i, j) += fine_displacement_correction(i, j, corrections[idx], pfine);
    }
  }
}

void refine(const grid_stack& pressures,
            grid_stack& displacements,
            const std::vector<matrix>& coarse_corrections,
            const std::vector<std::vector<matrix>>& fine_corrections) {
  for (int q = displacements.stack.size() - 1; q > 0; --q) {
    apply_coarse_correction(coarse_corrections[q-1],
                            pressures.stack[q-1],
                            displacements.stack[q]);
    refine(displacements.stack[q], displacements.stack[q-1],
           pressures.stencil);
    apply_fine_correction(fine_corrections[q-1],
                          pressures.stack[q-1],
                          displacements.stack[q-1]);
  }
}
