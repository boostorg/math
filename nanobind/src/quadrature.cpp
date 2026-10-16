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
#include <boost/math/quadrature/gauss_kronrod.hpp>
#include <boost/math/quadrature/ooura_fourier_integrals.hpp>
#include <boost/math/quadrature/sinh_sinh.hpp>
#include <boost/math/quadrature/tanh_sinh.hpp>
#include <boost/math/quadrature/trapezoidal.hpp>
#include <boost/math/quadrature/van_den_bos_unit_square.hpp>
#include <boost/math/tools/precision.hpp>

#include <nanobind/stl/pair.h>
#include <nanobind/stl/tuple.h>

#include <cstddef>
#include <limits>
#include <string>
#include <tuple>
#include <utility>

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
std::pair<double, double> gauss_hermite_n(scalar_function const& f)
{
   double L1 = 0;
   double Q = bq::gauss_hermite<double, N>::integrate(f, &L1);
   return {Q, L1};
}

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
   }, arg("f"), arg("points") = 10u, "quadrature::gauss_hermite<double, points>: (integral of exp(-x^2) f(x) over the real line, L1 norm)");

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
