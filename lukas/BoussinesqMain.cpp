#include "Boussinesq.hh"
#include "BoussinesqMatrix.hh"
#include "BoussinesqMlms.hh"

#include <iostream>
#include <iomanip>

int main() {
  double Lx = 2, Ly = 2;
  std::size_t Nx = 1048, Ny = 1048;

  vec2d pixel = {Lx / Nx, Ly / Ny};

  matrix pressure({Nx, Ny}), displacement({Nx, Ny});

  // initialize pressure
  for (std::size_t i = 0; i < pressure.shape[0]; ++i) {
    for (std::size_t j = 0; j < pressure.shape[1]; ++j) {
      if (i >= Nx/4 and i < 3 * Nx / 4 and j >= Nx/4 and j < 3 * Ny / 4) {
        pressure(i, j) = 1.;
      } else {
        pressure(i, j) = 0.;
      }
    }
  }

  // compute deflection
  // naive(pressure, displacement, pixel, pixel);
  // std::cout << displacement;
  // return 0;

  int t = 4, levels = optimal_levels(pressure.shape);

  grid_stack pressure_stack(pressure, pixel, t, levels),
    displacement_stack(displacement, pixel, t, levels);
  auto mc = correction_size(t, pressure.shape);

  pressure_stack.coarsen();

  // std::cout << pressure_stack.stack.back();
  compute_coarse_displacements(pressure_stack, displacement_stack);
  // std::cout << displacement_stack.stack.back();

  auto corrections = coarse_corrections(pressure_stack, mc);
  auto fcorrections = fine_corrections(displacement_stack, mc);

  // for (auto&& C : corrections)
  //   std:: cout << C << "-----\n";


  // for (auto&& Cs : fcorrections)
  //   for (auto&& C : Cs)
  //     std:: cout << C << "-----\n";

  refine(pressure_stack, displacement_stack, corrections, fcorrections);
  // std::cout << displacement_stack.stack.front();
  return 0;
}
