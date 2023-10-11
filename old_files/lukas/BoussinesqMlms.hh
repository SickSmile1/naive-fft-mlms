#include "BoussinesqMatrix.hh"

#include <vector>

std::size_t optimal_levels(std::array<std::size_t, 2> shape);
std::vector<double> make_stencil(int t);
std::size_t correction_size(int t, std::array<std::size_t, 2> shape);

void interpolate(const matrix& fine, matrix& coarse, std::vector<double> stencil);

struct grid_stack {
  grid_stack(const matrix& fine, vec2d pixel, int t, std::size_t levels);

  void add_level();
  void coarsen();

  int t;
  std::vector<double> stencil;
  std::vector<matrix> stack;
  std::vector<vec2d> pixel_sizes;
};

void compute_coarse_displacements(const grid_stack& pressure, grid_stack& displacement);

std::vector<matrix> coarse_corrections(const grid_stack& pressure, std::size_t mc);
std::vector<std::vector<matrix>> fine_corrections(const grid_stack& stack, std::size_t mc);

void apply_coarse_correction(const matrix& correction,
                             const matrix& pfine,
                             matrix& ucoarse);

void refine(const matrix& coarse, matrix& fine, std::vector<double> stencil);


void apply_fine_correction(const std::vector<matrix>& corrections,
                           const matrix& pfine,
                           matrix& ufine);

void refine(const grid_stack& pressures,
            grid_stack& displacements,
            const std::vector<matrix>& coarse_corrections,
            const std::vector<std::vector<matrix>>& fine_corrections);
