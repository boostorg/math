// Copyright 2026 Nick Thompson
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt
//
// Numerical differentiation of Python callables, and Lanczos smoothing derivatives of samples. (autodiff is
// template-only.) The C++ checks its arguments with assertions, so they are checked here first.

#include <boost/math/differentiation/finite_difference.hpp>
#include <boost/math/differentiation/lanczos_smoothing.hpp>

#include <nanobind/stl/complex.h>
#include <nanobind/stl/pair.h>

#include <cmath>
#include <complex>
#include <cstddef>
#include <string>
#include <utility>
#include <variant>
#include <vector>

#include "bindings.hpp"
#include "macros.hpp"

namespace {

namespace nb = nanobind;

namespace bd = boost::math::differentiation;

struct py_lanczos {
   double spacing;
   std::size_t n, approximation_order;
   int order;
   std::variant<bd::discrete_lanczos_derivative<double, 1>, bd::discrete_lanczos_derivative<double, 2>> d;

   static auto make(double spacing, std::size_t n, std::size_t p, int order) {
      if (!(spacing > 0) || !std::isfinite(spacing))
         throw nb::value_error("spacing must be positive and finite");
      if (order != 1 && order != 2)
         throw nb::value_error("order must be 1 or 2");
      if (n < 1 || p > 2 * n || p < std::size_t(order == 1 ? 2 : 3))
         throw nb::value_error(order == 1 ? "need n >= 1 and 2 <= approximation_order <= 2n"
                                          : "need n >= 1 and 3 <= approximation_order <= 2n for order=2");
      using V = decltype(d);
      return order == 1 ? V(std::in_place_index<0>, spacing, n, p) : V(std::in_place_index<1>, spacing, n, p);
   }
   py_lanczos(double spacing_, std::size_t n_, std::size_t p, int order_)
       : spacing(spacing_), n(n_), approximation_order(p), order(order_), d(make(spacing_, n_, p, order_)) {}

   void check_length(std::size_t size) const {
      if (size < 2 * n + 1)
         throw nb::value_error(("need at least 2n + 1 = " + std::to_string(2 * n + 1) + " samples").c_str());
   }
};

}  // namespace

void bmpy::bind_differentiation(nb::module_& m) {
   m.def("finite_difference_derivative", [](nb::callable f, double x, int order) {
      py_callable<double, double> g{f};
      double error = 0;
      auto fd = [&](auto o) { return bd::finite_difference_derivative<decltype(g), double, o()>(g, x, &error); };
      double d = 0;
      switch (order) {
      case 1: d = fd(std::integral_constant<std::size_t, 1>()); break;
      case 2: d = fd(std::integral_constant<std::size_t, 2>()); break;
      case 4: d = fd(std::integral_constant<std::size_t, 4>()); break;
      case 6: d = fd(std::integral_constant<std::size_t, 6>()); break;
      case 8: d = fd(std::integral_constant<std::size_t, 8>()); break;
      default: throw nb::value_error(("order must be 1, 2, 4, 6 or 8, got " + std::to_string(order)).c_str());
      }
      return std::make_pair(d, error);
   }, nb::arg("f"), nb::arg("x"), nb::kw_only(), nb::arg("order") = 6,
      nb::sig("def finite_difference_derivative(f, x, /, *, order=6) -> tuple[float, float]"),
      "(f'(x), error estimate) by a central finite difference of the given order (1, 2, 4, 6 or 8)." BMPY_URL("diff"));

   m.def("complex_step_derivative", [](nb::callable f, double x) {
      return bd::complex_step_derivative(py_callable<std::complex<double>, std::complex<double>>{f}, x);
   }, nb::arg("f"), nb::arg("x"), nb::sig("def complex_step_derivative(f, x, /) -> float"),
      "f'(x) = Im f(x + ih)/h, h tiny; f must be holomorphic, real on the real line, and take complex (cmath, not math)."
      BMPY_URL("diff"));

   nb::class_<py_lanczos>(m, "discrete_lanczos_derivative",
      "The Lanczos smoothing derivative (order 1 or 2) of noisy samples `spacing` apart, fitted over 2n + 1 samples."
      BMPY_URL("diff0"))
      .def(nb::init<double, std::size_t, std::size_t, int>(), nb::arg("spacing"), nb::kw_only(), nb::arg("n") = 18,
           nb::arg("approximation_order") = 3, nb::arg("order") = 1,
           nb::sig("def __init__(self, spacing, /, *, n=18, approximation_order=3, order=1) -> None"))
      .def_ro("spacing", &py_lanczos::spacing)
      .def_ro("n", &py_lanczos::n)
      .def_ro("approximation_order", &py_lanczos::approximation_order)
      .def_ro("order", &py_lanczos::order)
      .def("__call__", [](py_lanczos const& L, nb::handle v) {
         std::vector<double> x = array_to_vector(v);
         L.check_length(x.size());
         std::vector<double> w(x.size());
         std::visit([&](auto const& d) { d(x, w); }, L.d);
         return vector_to_array(std::move(w));
      }, nb::sig("def __call__(self, v, /)"), "The derivative at every sample of v, as an array.")
      .def("__call__", [](py_lanczos const& L, nb::handle v, std::size_t i) {
         std::vector<double> x = array_to_vector(v);
         L.check_length(x.size());
         if (i >= x.size())
            throw nb::index_error("sample index out of range");
         return std::visit([&](auto const& d) { return d(x, i); }, L.d);
      }, nb::sig("def __call__(self, v, i, /) -> float"), "The derivative at sample i of v.")
      .def("__repr__", [](py_lanczos const& L) {
         return "discrete_lanczos_derivative(" + py_float_repr(L.spacing) + ", n=" + std::to_string(L.n)
                + ", approximation_order=" + std::to_string(L.approximation_order) + ", order=" + std::to_string(L.order) + ")";
      });
}
