#include "pybind11/pybind11.h"
#include "pybind11/eigen.h"
#include <Eigen/Core>
#include "BoussinesqFft.h"
#include "BoussinesqMlms.h"

namespace py = pybind11;

PYBIND11_MODULE(boussinesq, mod) {
  mod.doc() = "Boussinesq calculation of rough surface with FFT or MLMS";

  mod.def("FFT", py::overload_cast<double, matrix, matrix>(&BoussinesqFFT), 
      "Calculate with FFT, args: size(double), surf(2darray), topo(2darray)");
  mod.def("MLMS", py::overload_cast<double, matrix, matrix, int>(&BoussinesqMlms), 
      "Calculate with MLMS, args: size(double), surf(2darray), topo(2darray), int t");
}
