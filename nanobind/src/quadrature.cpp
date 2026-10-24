// Copyright 2026 Nick Thompson
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt
//
// The one-dimensional quadrature routines (and the unit-square cubature), in double, taking a Python
// callable. Each returns what the C++ routine reports: the integral, then its error estimate and L1 norm
// where the C++ gives them. The integrand is called with the GIL held, and an exception it raises
// propagates. naive_monte_carlo is left out: it calls the integrand from its own threads.

#include "bindings.hpp"

#include <boost/math/quadrature/exp_sinh.hpp>
#include <boost/math/quadrature/gauss.hpp>
#include <boost/math/quadrature/gauss_hermite.hpp>
#include <boost/math/quadrature/gauss_jacobi.hpp>
#include <boost/math/quadrature/gauss_kronrod.hpp>
#include <boost/math/quadrature/gauss_laguerre.hpp>
#include <boost/math/quadrature/ooura_fourier_integrals.hpp>
#include <boost/math/quadrature/sinh_sinh.hpp>
#include <boost/math/quadrature/tanh_sinh.hpp>
#include <boost/math/quadrature/trapezoidal.hpp>
#include <boost/math/quadrature/van_den_bos_unit_square.hpp>
#include <boost/math/tools/precision.hpp>

#include <nanobind/ndarray.h>
#include <nanobind/stl/pair.h>
#include <nanobind/stl/tuple.h>

#include <cstddef>
#include <limits>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace bmpy {

namespace {

namespace bq = boost::math::quadrature;

// f(x) for a Python callable f that returns something convertible to float.
struct scalar_function
{
   nb::callable f;
   double operator()(double x) const { return nb::cast<double>(f(x)); }
};

struct planar_function
{
   nb::callable f;
   double operator()(double x, double y) const { return nb::cast<double>(f(x, y)); }
};

const double root_eps = boost::math::tools::root_epsilon<double>();

template <unsigned N>
std::pair<double, double> gauss_n(scalar_function const& f, double a, double b)
{
   double L1 = 0;
   double Q = bq::gauss<double, N>::integrate(f, a, b, &L1);
   return {Q, L1};
}

template <unsigned N>
std::tuple<double, double, double> gauss_kronrod_n(scalar_function const& f, double a, double b, unsigned max_depth, double tolerance)
{
   double error = 0, L1 = 0;
   double Q = bq::gauss_kronrod<double, N>::integrate(f, a, b, max_depth, tolerance, &error, &L1);
   return {Q, error, L1};
}

template <unsigned N>
std::tuple<double, double, double> gauss_hermite_n(scalar_function const& f)
{
   double error = 0, L1 = 0;
   double Q = bq::gauss_hermite<double, N>::integrate(f, &error, &L1);
   return {Q, error, L1};
}

template <unsigned N>
std::tuple<double, double, double> gauss_laguerre_n(scalar_function const& f, double alpha)
{
   double error = 0, L1 = 0;
   double Q = bq::gauss_laguerre<double, N>(alpha).integrate(f, &error, &L1);
   return {Q, error, L1};
}

using array1d = nb::ndarray<nb::numpy, double, nb::ndim<1>>;

// A copy, so that writing to it leaves the rule alone.
array1d to_array(std::vector<double> const& v)
{
   auto* data = new std::vector<double>(v);
   nb::capsule owner(data, [](void* p) noexcept { delete static_cast<std::vector<double>*>(p); });
   return array1d(data->data(), {data->size()}, owner);
}

// gauss_jacobi keeps its parameters to itself; remember them for repr.
struct py_gauss_jacobi : bq::gauss_jacobi<double>
{
   unsigned n;
   double alpha, beta;
   py_gauss_jacobi(unsigned n_, double alpha_, double beta_) : bq::gauss_jacobi<double>(n_, alpha_, beta_), n(n_), alpha(alpha_), beta(beta_) {}
};

[[noreturn]] void bad_points(char const* name, char const* allowed, unsigned points)
{
   throw nb::value_error((std::string(name) + ": points must be one of " + allowed + ", got " + std::to_string(points)).c_str());
}

}  // namespace

void bind_quadrature(nb::module_& m)
{
   using nb::arg;

   m.def("tanh_sinh", [](nb::callable f, double a, double b, double tolerance, std::size_t max_refinements) {
      double error = 0, L1 = 0;
      bq::tanh_sinh<double> integrator(max_refinements);
      double Q = integrator.integrate(scalar_function{f}, a, b, tolerance, &error, &L1);
      return std::make_tuple(Q, error, L1);
   }, arg("f"), arg("a"), arg("b"), arg("tolerance") = root_eps, arg("max_refinements") = 15,
      "quadrature::tanh_sinh: (integral, error estimate, L1 norm) of f over [a, b], endpoint singularities allowed");

   m.def("exp_sinh", [](nb::callable f, double a, double b, double tolerance, std::size_t max_refinements) {
      double error = 0, L1 = 0;
      bq::exp_sinh<double> integrator(max_refinements);
      double Q = integrator.integrate(scalar_function{f}, a, b, tolerance, &error, &L1);
      return std::make_tuple(Q, error, L1);
   }, arg("f"), arg("a") = 0.0, arg("b") = std::numeric_limits<double>::infinity(), arg("tolerance") = root_eps,
      arg("max_refinements") = 9, "quadrature::exp_sinh: (integral, error estimate, L1 norm) over a half-infinite interval");

   m.def("sinh_sinh", [](nb::callable f, double tolerance, std::size_t max_refinements) {
      double error = 0, L1 = 0;
      bq::sinh_sinh<double> integrator(max_refinements);
      double Q = integrator.integrate(scalar_function{f}, tolerance, &error, &L1);
      return std::make_tuple(Q, error, L1);
   }, arg("f"), arg("tolerance") = root_eps, arg("max_refinements") = 9,
      "quadrature::sinh_sinh: (integral, error estimate, L1 norm) of f over the whole real line");

   m.def("trapezoidal", [](nb::callable f, double a, double b, double tolerance, std::size_t max_refinements) {
      double error = 0, L1 = 0;
      double Q = bq::trapezoidal(scalar_function{f}, a, b, tolerance, max_refinements, &error, &L1);
      return std::make_tuple(Q, error, L1);
   }, arg("f"), arg("a"), arg("b"), arg("tolerance") = root_eps, arg("max_refinements") = 12,
      "quadrature::trapezoidal: (integral, error estimate, L1 norm); for periodic integrands over a period");

   m.def("gauss", [](nb::callable f, double a, double b, unsigned points) {
      scalar_function g{f};
      switch (points)
      {
      case 7: return gauss_n<7>(g, a, b);
      case 10: return gauss_n<10>(g, a, b);
      case 15: return gauss_n<15>(g, a, b);
      case 20: return gauss_n<20>(g, a, b);
      case 25: return gauss_n<25>(g, a, b);
      case 30: return gauss_n<30>(g, a, b);
      default: bad_points("gauss", "7, 10, 15, 20, 25, 30", points);
      }
   }, arg("f"), arg("a"), arg("b"), arg("points") = 7u, "quadrature::gauss<double, points>: (integral, L1 norm)");

   m.def("gauss_kronrod", [](nb::callable f, double a, double b, unsigned points, unsigned max_depth, double tolerance) {
      scalar_function g{f};
      switch (points)
      {
      case 15: return gauss_kronrod_n<15>(g, a, b, max_depth, tolerance);
      case 21: return gauss_kronrod_n<21>(g, a, b, max_depth, tolerance);
      case 31: return gauss_kronrod_n<31>(g, a, b, max_depth, tolerance);
      case 41: return gauss_kronrod_n<41>(g, a, b, max_depth, tolerance);
      case 51: return gauss_kronrod_n<51>(g, a, b, max_depth, tolerance);
      case 61: return gauss_kronrod_n<61>(g, a, b, max_depth, tolerance);
      default: bad_points("gauss_kronrod", "15, 21, 31, 41, 51, 61", points);
      }
   }, arg("f"), arg("a"), arg("b"), arg("points") = 15u, arg("max_depth") = 15u, arg("tolerance") = root_eps,
      "quadrature::gauss_kronrod<double, points>: (integral, error estimate, L1 norm), adaptive");

   m.def("gauss_hermite", [](nb::callable f, unsigned points) {
      scalar_function g{f};
      switch (points)
      {
      case 7: return gauss_hermite_n<7>(g);
      case 10: return gauss_hermite_n<10>(g);
      case 15: return gauss_hermite_n<15>(g);
      case 20: return gauss_hermite_n<20>(g);
      case 30: return gauss_hermite_n<30>(g);
      default: bad_points("gauss_hermite", "7, 10, 15, 20, 30", points);
      }
   }, arg("f"), arg("points") = 10u, "quadrature::gauss_hermite<double, points>: (integral of exp(-x^2) f(x) over the real line, error estimate, L1 norm)");

   m.def("gauss_laguerre", [](nb::callable f, unsigned points, double alpha) {
      scalar_function g{f};
      switch (points)
      {
      case 7: return gauss_laguerre_n<7>(g, alpha);
      case 10: return gauss_laguerre_n<10>(g, alpha);
      case 15: return gauss_laguerre_n<15>(g, alpha);
      case 20: return gauss_laguerre_n<20>(g, alpha);
      case 30: return gauss_laguerre_n<30>(g, alpha);
      default: bad_points("gauss_laguerre", "7, 10, 15, 20, 30", points);
      }
   }, arg("f"), arg("points") = 10u, arg("alpha") = 0.0,
      "quadrature::gauss_laguerre<double, points>(alpha): (integral of x^alpha exp(-x) f(x) over [0, inf), error estimate, L1 norm)");

   nb::class_<py_gauss_jacobi>(m, "gauss_jacobi",
      "quadrature::gauss_jacobi<double>(n, alpha, beta): the n-point rule for (1 - x)^alpha (1 + x)^beta on [-1, 1], "
      "and Laurie's (n + 1)-point anti-Gauss rule")
      .def(nb::init<unsigned, double, double>(), arg("n"), arg("alpha"), arg("beta"))
      .def("integrate", [](py_gauss_jacobi const& rule, nb::callable f) {
         double error = 0, L1 = 0;
         double Q = rule.integrate(scalar_function{f}, &error, &L1);
         return std::make_tuple(Q, error, L1);
      }, arg("f"), "(integral, error estimate, L1 norm); the estimate is NaN when anti_gauss_is_internal is False")
      .def_ro("n", &py_gauss_jacobi::n)
      .def_ro("alpha", &py_gauss_jacobi::alpha)
      .def_ro("beta", &py_gauss_jacobi::beta)
      .def_prop_ro("abscissa", [](py_gauss_jacobi const& rule) { return to_array(rule.abscissa()); }, nb::rv_policy::move)
      .def_prop_ro("weights", [](py_gauss_jacobi const& rule) { return to_array(rule.weights()); }, nb::rv_policy::move)
      .def_prop_ro("anti_gauss_abscissa", [](py_gauss_jacobi const& rule) { return to_array(rule.anti_gauss_abscissa()); }, nb::rv_policy::move)
      .def_prop_ro("anti_gauss_weights", [](py_gauss_jacobi const& rule) { return to_array(rule.anti_gauss_weights()); }, nb::rv_policy::move)
      .def_prop_ro("anti_gauss_is_internal", [](py_gauss_jacobi const& rule) { return rule.anti_gauss_is_internal(); },
         "False when the anti-Gauss rule would have a node outside [-1, 1]; its nodes and weights are then empty")
      .def("__repr__", [](py_gauss_jacobi const& rule) {
         return "gauss_jacobi(n=" + std::to_string(rule.n) + ", alpha=" + nb::cast<std::string>(nb::repr(nb::float_(rule.alpha)))
            + ", beta=" + nb::cast<std::string>(nb::repr(nb::float_(rule.beta))) + ")";
      });

   m.def("ooura_fourier_sin", [](nb::callable f, double omega, double tolerance, std::size_t levels) {
      bq::ooura_fourier_sin<double> integrator(tolerance, levels);
      return integrator.integrate(scalar_function{f}, omega);
   }, arg("f"), arg("omega"), arg("tolerance") = root_eps, arg("levels") = 8,
      "quadrature::ooura_fourier_sin: (integral of f(t) sin(omega t) over [0, inf), relative error estimate)");

   m.def("ooura_fourier_cos", [](nb::callable f, double omega, double tolerance, std::size_t levels) {
      bq::ooura_fourier_cos<double> integrator(tolerance, levels);
      return integrator.integrate(scalar_function{f}, omega);
   }, arg("f"), arg("omega"), arg("tolerance") = root_eps, arg("levels") = 8,
      "quadrature::ooura_fourier_cos: (integral of f(t) cos(omega t) over [0, inf), relative error estimate)");

   m.def("van_den_bos_unit_square", [](nb::callable f, double tolerance) {
      double error = 0, L1 = 0;
      double Q = bq::van_den_bos_unit_square(planar_function{f}, tolerance, &error, &L1);
      return std::make_tuple(Q, error, L1);
   }, arg("f"), arg("tolerance") = root_eps,
      "quadrature::van_den_bos_unit_square: (integral of f(x, y) over [0, 1]^2, error estimate, L1 norm)");
}

}  // namespace bmpy
