// Copyright 2026 Nick Thompson
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt
//
// The one-dimensional quadrature routines (and the unit-square cubature), in double, taking a Python
// callable. Each returns what the C++ routine reports, as a named tuple: QuadratureResult(value, error,
// l1_norm), or GaussResult(value, l1_norm) and OouraResult(value, error) where the C++ gives less. The
// integrand is called with the GIL held, and an exception it raises propagates. naive_monte_carlo is in
// naive_monte_carlo.cpp.

#include "bindings.hpp"

#include <boost/math/quadrature/exp_sinh.hpp>
#include <boost/math/quadrature/gauss.hpp>
#include <boost/math/quadrature/gauss_hermite.hpp>
#include <boost/math/quadrature/gauss_jacobi.hpp>
#include <boost/math/quadrature/gauss_kronrod.hpp>
#include <boost/math/quadrature/gauss_laguerre.hpp>
#include <boost/math/quadrature/leja.hpp>
#include <boost/math/quadrature/ooura_fourier_integrals.hpp>
#include <boost/math/quadrature/sinh_sinh.hpp>
#include <boost/math/quadrature/sparse_grid.hpp>
#include <boost/math/quadrature/tanh_sinh.hpp>
#include <boost/math/quadrature/trapezoidal.hpp>
#include <boost/math/quadrature/van_den_bos_unit_square.hpp>
#include <boost/math/tools/precision.hpp>

#include <nanobind/ndarray.h>
#include <nanobind/stl/pair.h>
#include <nanobind/stl/tuple.h>
#include <nanobind/stl/vector.h>

#include <cmath>
#include <cstddef>
#include <limits>
#include <string>
#include <tuple>
#include <type_traits>
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

// The result types, made in bind_quadrature; the module keeps them alive. The C++ types name them in the type stubs.
nb::handle quadrature_result, gauss_result, ooura_result;
bool any_object(PyObject*) { return true; }
struct QuadratureResult : nb::object { NB_OBJECT_DEFAULT(QuadratureResult, nb::object, "QuadratureResult", any_object) };
struct GaussResult : nb::object { NB_OBJECT_DEFAULT(GaussResult, nb::object, "GaussResult", any_object) };
struct OouraResult : nb::object { NB_OBJECT_DEFAULT(OouraResult, nb::object, "OouraResult", any_object) };

QuadratureResult result(double Q, double error, double L1) { return nb::borrow<QuadratureResult>(quadrature_result(Q, error, L1)); }

template <unsigned N>
GaussResult gauss_n(scalar_function const& f, double a, double b)
{
   double L1 = 0;
   double Q = bq::gauss<double, N>::integrate(f, a, b, &L1);
   return nb::borrow<GaussResult>(gauss_result(Q, L1));
}

template <unsigned N>
QuadratureResult gauss_kronrod_n(scalar_function const& f, double a, double b, unsigned max_depth, double tolerance)
{
   double error = 0, L1 = 0;
   double Q = bq::gauss_kronrod<double, N>::integrate(f, a, b, max_depth, tolerance, &error, &L1);
   return result(Q, error, L1);
}

// (integral, L1 norm) of the Gauss rule, or (integral, error estimate, L1 norm) of the averaged rule.
template <class Rule>
GaussResult gauss_rule_integrate(Rule const& rule, scalar_function const& f)
{
   double L1 = 0;
   double Q = rule.integrate(f, &L1);
   return nb::borrow<GaussResult>(gauss_result(Q, L1));
}

template <class Rule>
QuadratureResult gauss_rule_integrate_averaged(Rule const& rule, scalar_function const& f)
{
   double error = 0, L1 = 0;
   double Q = rule.integrate_averaged(f, &error, &L1);
   return result(Q, error, L1);
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

// The Leja rules with parameters keep them to themselves; remember them for repr.
struct py_leja_laguerre : bq::leja_laguerre<double>
{
   double alpha;
   explicit py_leja_laguerre(double alpha_) : bq::leja_laguerre<double>(alpha_), alpha(alpha_) {}
};

struct py_leja_jacobi : bq::leja_jacobi<double>
{
   double alpha, beta;
   py_leja_jacobi(double alpha_, double beta_) : bq::leja_jacobi<double>(alpha_, beta_), alpha(alpha_), beta(beta_) {}
};

std::string py_repr(double x) { return nb::cast<std::string>(nb::repr(nb::float_(x))); }

// The n-point Gauss rule of a family, computed afresh: (nodes ascending, weights). A symmetric family keeps only
// its nonnegative nodes, so mirror them.
template <class Family>
std::pair<array1d, array1d> gauss_rule(unsigned n, Family const& family, bool symmetric = false)
{
   BOOST_MATH_NAMESPACE::detail::orthogonal_polynomial<double, Family> rule(n, family);
   std::vector<double> x = rule.abscissa(), w = rule.weights();
   if (symmetric)
   {
      std::size_t m = x.size(), skip = n % 2;  // An odd rule has 0 once.
      std::vector<double> xs(n), ws(n);
      for (std::size_t i = 0; i < m; ++i)
      {
         xs[n - m + i] = x[i];
         ws[n - m + i] = w[i];
         if (i >= skip)
         {
            xs[m - 1 - i] = -x[i];
            ws[m - 1 - i] = w[i];
         }
      }
      x.swap(xs);
      w.swap(ws);
   }
   return {to_array(x), to_array(w)};
}

// A nested Leja rule, whose integrate adds points until the tolerance; over [a, b] for the rules on [-1, 1].
template <class Rule>
nb::class_<Rule> bind_leja(nb::module_& m, char const* name, char const* doc)
{
   nb::class_<Rule> c(m, name, doc);
   if constexpr (bq::detail::leja_bounded<Rule>::value)
      c.def("integrate", [](Rule const& rule, nb::callable f, double a, double b, double tolerance) {
         double error = 0, L1 = 0;
         double Q = rule.integrate(scalar_function{f}, a, b, tolerance, &error, &L1);
         return result(Q, error, L1);
      }, nb::arg("f"), nb::arg("a"), nb::arg("b"), nb::kw_only(), nb::arg("tolerance") = root_eps, "The QuadratureResult of f over [a, b], adding one point at a time until the error is below tolerance times the L1 norm."
      BMPY_URL("leja"));
   return c
      .def("integrate", [](Rule const& rule, nb::callable f, double tolerance) {
         double error = 0, L1 = 0;
         double Q = rule.integrate(scalar_function{f}, tolerance, &error, &L1);
         return result(Q, error, L1);
      }, nb::arg("f"), nb::kw_only(), nb::arg("tolerance") = root_eps, "The QuadratureResult of f against the weight, adding one point at a time until the error is below tolerance times the L1 norm."
      BMPY_URL("leja"))
      .def_prop_ro("max_points", [](Rule const&) { return Rule::max_points(); })
      .def_prop_ro("abscissa", [](Rule const& rule) { return to_array(rule.abscissa()); }, nb::rv_policy::move, "the points computed so far, in the order chosen")
      .def("weights", [](Rule const& rule, unsigned n) { return to_array(rule.weights(n)); }, nb::kw_only(), nb::arg("n"), nb::rv_policy::move,
         "the weights of the n-point rule, whose points are abscissa[:n]");
}

// sparse_grid(f, rule=..., ...), one overload per rule; over a box, the rule defaults to leja_legendre.
template <class Rule>
void bind_sparse_grid(nb::module_& m, nb::object default_rule = nb::object())
{
   m.def("sparse_grid", [](nb::callable f, Rule const& rule, unsigned dimension, std::size_t max_calls, double tolerance) {
      double error = 0, L1 = 0;
      double Q = bq::sparse_grid([&f](std::vector<double> const& x) { return nb::cast<double>(f(to_array(x))); },
                                 rule, dimension, max_calls, tolerance, &error, &L1);
      return result(Q, error, L1);
   }, nb::arg("f"), nb::kw_only(), nb::arg("rule"), nb::arg("dimension"), nb::arg("max_calls") = 100000, nb::arg("tolerance") = root_eps,
      "The QuadratureResult of f over R^dimension by the dimension-adaptive sparse grid of the rule; f takes a 1D array."
      BMPY_URL("sparse_grid"));
   if constexpr (bq::detail::leja_bounded<Rule>::value)
   {
      auto box = [](nb::callable f, std::vector<double> const& a, std::vector<double> const& b, Rule const& rule,
                    std::size_t max_calls, double tolerance) {
         double error = 0, L1 = 0;
         double Q = bq::sparse_grid([&f](std::vector<double> const& x) { return nb::cast<double>(f(to_array(x))); },
                                    rule, a, b, max_calls, tolerance, &error, &L1);
         return result(Q, error, L1);
      };
      char const* doc = "The QuadratureResult of f over the box [a, b] by the dimension-adaptive sparse grid; the dimension is len(a)."
      BMPY_URL("sparse_grid");
      if (default_rule)
         m.def("sparse_grid", box, nb::arg("f"), nb::arg("a"), nb::arg("b"), nb::kw_only(), nb::arg("rule") = default_rule,
               nb::arg("max_calls") = 100000, nb::arg("tolerance") = root_eps, doc);
      else
         m.def("sparse_grid", box, nb::arg("f"), nb::arg("a"), nb::arg("b"), nb::kw_only(), nb::arg("rule"),
               nb::arg("max_calls") = 100000, nb::arg("tolerance") = root_eps, doc);
   }
}

void check_exponent(char const* name, double a)
{
   if (!(a > -1) || !std::isfinite(a))
      throw nb::value_error((std::string(name) + " must be finite and > -1").c_str());
}

[[noreturn]] void bad_points(char const* name, char const* allowed, unsigned points)
{
   throw nb::value_error((std::string(name) + ": points must be one of " + allowed + ", got " + std::to_string(points)).c_str());
}

// fn(std::integral_constant<unsigned, points>()) for the point counts gauss_hermite and gauss_laguerre are bound for.
template <class Fn>
auto with_points(char const* name, unsigned points, Fn fn)
{
   switch (points)
   {
   case 7: return fn(std::integral_constant<unsigned, 7>());
   case 10: return fn(std::integral_constant<unsigned, 10>());
   case 15: return fn(std::integral_constant<unsigned, 15>());
   case 20: return fn(std::integral_constant<unsigned, 20>());
   case 30: return fn(std::integral_constant<unsigned, 30>());
   default: bad_points(name, "7, 10, 15, 20, 30", points);
   }
}

}  // namespace

void bind_quadrature(nb::module_& m)
{
   using nb::arg;

   quadrature_result = named_tuple(m, "QuadratureResult", "value error l1_norm", "An integral, its error estimate and its L1 norm.");
   gauss_result = named_tuple(m, "GaussResult", "value l1_norm", "An integral by a fixed Gauss rule, and its L1 norm.");
   ooura_result = named_tuple(m, "OouraResult", "value error", "A Fourier integral, and its relative error estimate.");

   m.def("tanh_sinh", [](nb::callable f, double a, double b, double tolerance, std::size_t max_refinements) {
      double error = 0, L1 = 0;
      bq::tanh_sinh<double> integrator(max_refinements);
      double Q = integrator.integrate(scalar_function{f}, a, b, tolerance, &error, &L1);
      return result(Q, error, L1);
   }, arg("f"), arg("a"), arg("b"), arg("tolerance") = root_eps, arg("max_refinements") = 15,
      "The QuadratureResult of f over [a, b] by tanh-sinh quadrature; endpoint singularities are allowed."
      BMPY_URL("double_exponential/de_tanh_sinh"));

   m.def("exp_sinh", [](nb::callable f, double a, double b, double tolerance, std::size_t max_refinements) {
      double error = 0, L1 = 0;
      bq::exp_sinh<double> integrator(max_refinements);
      double Q = integrator.integrate(scalar_function{f}, a, b, tolerance, &error, &L1);
      return result(Q, error, L1);
   }, arg("f"), arg("a") = 0.0, arg("b") = std::numeric_limits<double>::infinity(), arg("tolerance") = root_eps,
      arg("max_refinements") = 9, "The QuadratureResult of f over a half-infinite interval [a, inf) or (-inf, b] by exp-sinh quadrature."
      BMPY_URL("double_exponential/de_exp_sinh"));

   m.def("sinh_sinh", [](nb::callable f, double tolerance, std::size_t max_refinements) {
      double error = 0, L1 = 0;
      bq::sinh_sinh<double> integrator(max_refinements);
      double Q = integrator.integrate(scalar_function{f}, tolerance, &error, &L1);
      return result(Q, error, L1);
   }, arg("f"), arg("tolerance") = root_eps, arg("max_refinements") = 9,
      "The QuadratureResult of f over the whole real line by sinh-sinh quadrature." BMPY_URL("double_exponential/de_sinh_sinh"));

   m.def("trapezoidal", [](nb::callable f, double a, double b, double tolerance, std::size_t max_refinements) {
      double error = 0, L1 = 0;
      double Q = bq::trapezoidal(scalar_function{f}, a, b, tolerance, max_refinements, &error, &L1);
      return result(Q, error, L1);
   }, arg("f"), arg("a"), arg("b"), arg("tolerance") = root_eps, arg("max_refinements") = 12,
      "The QuadratureResult of f over [a, b] by the adaptive trapezoidal rule, exponentially convergent for periodic f over a period."
      BMPY_URL("trapezoidal"));

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
   }, arg("f"), arg("a"), arg("b"), arg("points") = 7u, "The GaussResult of f over [a, b] by the Gauss-Legendre rule of `points` points (7, 10, 15, 20, 25 or 30)."
      BMPY_URL("gauss"));

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
      "The QuadratureResult of f over [a, b] by adaptive Gauss-Kronrod quadrature, bisecting to depth max_depth."
      BMPY_URL("gauss_kronrod"));

   m.def("gauss_hermite", [](nb::callable f, unsigned points) {
      return with_points("gauss_hermite", points, [g = scalar_function{f}](auto n) {
         return gauss_rule_integrate(bq::gauss_hermite<double, decltype(n)::value>(), g);
      });
   }, arg("f"), arg("points") = 10u, "The GaussResult of the integral of exp(-x^2) f(x) over the real line by the Gauss-Hermite rule of `points` points."
      BMPY_URL("gauss_hermite"));

   m.def("gauss_hermite_averaged", [](nb::callable f, unsigned points) {
      return with_points("gauss_hermite_averaged", points, [g = scalar_function{f}](auto n) {
         return gauss_rule_integrate_averaged(bq::gauss_hermite<double, decltype(n)::value>(), g);
      });
   }, arg("f"), arg("points") = 10u,
      "The QuadratureResult of the integral of exp(-x^2) f(x) by Laurie's average of the Gauss and anti-Gauss rules."
      BMPY_URL("gauss_hermite"));

   m.def("gauss_laguerre", [](nb::callable f, unsigned points, double alpha) {
      return with_points("gauss_laguerre", points, [g = scalar_function{f}, alpha](auto n) {
         return gauss_rule_integrate(bq::gauss_laguerre<double, decltype(n)::value>(alpha), g);
      });
   }, arg("f"), arg("points") = 10u, arg("alpha") = 0.0,
      "The GaussResult of the integral of x^alpha exp(-x) f(x) over [0, inf) by the Gauss-Laguerre rule of `points` points."
      "\n\nhttps://www.boost.org/doc/libs/release/libs/math/doc/html/quadrature.html");

   m.def("gauss_laguerre_averaged", [](nb::callable f, unsigned points, double alpha) {
      return with_points("gauss_laguerre_averaged", points, [g = scalar_function{f}, alpha](auto n) {
         return gauss_rule_integrate_averaged(bq::gauss_laguerre<double, decltype(n)::value>(alpha), g);
      });
   }, arg("f"), arg("points") = 10u, arg("alpha") = 0.0,
      "The QuadratureResult of the integral of x^alpha exp(-x) f(x) by Laurie's average of the Gauss and anti-Gauss rules."
      "\n\nhttps://www.boost.org/doc/libs/release/libs/math/doc/html/quadrature.html");

   nb::class_<py_gauss_jacobi>(m, "gauss_jacobi",
      "gauss_jacobi(n, alpha, beta): the n-point Gauss rule for the weight (1 - x)^alpha (1 + x)^beta on [-1, 1], "
      "and Laurie's (n + 1)-point anti-Gauss rule; alpha, beta > -1."
      "\n\nhttps://www.boost.org/doc/libs/release/libs/math/doc/html/quadrature.html")
      .def(nb::init<unsigned, double, double>(), arg("n"), arg("alpha"), arg("beta"))
      .def("integrate", [](py_gauss_jacobi const& rule, nb::callable f) {
         return gauss_rule_integrate(rule, scalar_function{f});
      }, arg("f"), "The GaussResult of the integral of (1 - x)^alpha (1 + x)^beta f(x) over [-1, 1].")
      .def("integrate_averaged", [](py_gauss_jacobi const& rule, nb::callable f) {
         return gauss_rule_integrate_averaged(rule, scalar_function{f});
      }, arg("f"), "The QuadratureResult of Laurie's average of the Gauss and anti-Gauss rules; raises when anti_gauss_is_internal is False.")
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

   bind_leja<bq::leja_legendre<double>>(m, "leja_legendre",
      "leja_legendre(): nested Leja rules on [-1, 1], or on [a, b] by an affine map." BMPY_URL("leja"))
      .def(nb::init<>())
      .def("__repr__", [](bq::leja_legendre<double> const&) { return std::string("leja_legendre()"); });
   bind_leja<bq::leja_hermite<double>>(m, "leja_hermite",
      "leja_hermite(): nested Leja rules for the weight exp(-x^2) on the real line." BMPY_URL("leja"))
      .def(nb::init<>())
      .def("__repr__", [](bq::leja_hermite<double> const&) { return std::string("leja_hermite()"); });
   bind_leja<py_leja_laguerre>(m, "leja_laguerre",
      "leja_laguerre(alpha=0): nested Leja rules for the weight x^alpha exp(-x) on [0, inf)." BMPY_URL("leja"))
      .def(nb::init<double>(), nb::kw_only(), nb::arg("alpha") = 0.0)
      .def_ro("alpha", &py_leja_laguerre::alpha)
      .def("__repr__", [](py_leja_laguerre const& rule) { return "leja_laguerre(alpha=" + py_repr(rule.alpha) + ")"; });
   bind_leja<py_leja_jacobi>(m, "leja_jacobi",
      "leja_jacobi(alpha, beta): nested Leja rules for the weight (1 - x)^alpha (1 + x)^beta on [-1, 1]." BMPY_URL("leja"))
      .def(nb::init<double, double>(), nb::kw_only(), nb::arg("alpha"), nb::arg("beta"))
      .def_ro("alpha", &py_leja_jacobi::alpha)
      .def_ro("beta", &py_leja_jacobi::beta)
      .def("__repr__", [](py_leja_jacobi const& rule) {
         return "leja_jacobi(alpha=" + py_repr(rule.alpha) + ", beta=" + py_repr(rule.beta) + ")";
      });
   bind_sparse_grid<bq::leja_legendre<double>>(m, nb::cast(bq::leja_legendre<double>()));
   bind_sparse_grid<bq::leja_hermite<double>>(m);
   bind_sparse_grid<py_leja_laguerre>(m);
   bind_sparse_grid<py_leja_jacobi>(m);

   m.def("gauss_hermite_rule", [](unsigned n) {
      return gauss_rule(n, BOOST_MATH_NAMESPACE::detail::hermite_family<double>(n), true);
   }, arg("n"), nb::rv_policy::move,
      "(nodes, weights) of the n-point Gauss rule for exp(-x^2), computed on each call; as numpy's hermgauss."
      "\n\nhttps://www.boost.org/doc/libs/release/libs/math/doc/html/quadrature.html");

   m.def("gauss_laguerre_rule", [](unsigned n, double alpha) {
      check_exponent("alpha", alpha);
      return gauss_rule(n, BOOST_MATH_NAMESPACE::detail::laguerre_family<double>(n, alpha));
   }, arg("n"), arg("alpha") = 0.0, nb::rv_policy::move,
      "(nodes, weights) of the n-point Gauss rule for x^alpha exp(-x), computed on each call; as scipy's roots_genlaguerre."
      "\n\nhttps://www.boost.org/doc/libs/release/libs/math/doc/html/quadrature.html");

   m.def("gauss_jacobi_rule", [](unsigned n, double alpha, double beta) {
      check_exponent("alpha", alpha);
      check_exponent("beta", beta);
      return gauss_rule(n, BOOST_MATH_NAMESPACE::detail::jacobi_family<double>(n, alpha, beta));
   }, arg("n"), arg("alpha"), arg("beta"), nb::rv_policy::move,
      "(nodes, weights) of the n-point Gauss rule for (1 - x)^alpha (1 + x)^beta, computed on each call, without the "
      "anti-Gauss rule; as scipy's roots_jacobi."
      "\n\nhttps://www.boost.org/doc/libs/release/libs/math/doc/html/quadrature.html");

   m.def("ooura_fourier_sin", [](nb::callable f, double omega, double tolerance, std::size_t levels) {
      if (levels < 1 || levels > 16)  // 0 reads an empty node table; each level doubles the work.
         throw nb::value_error("ooura_fourier_sin: levels must be in [1, 16]");
      bq::ooura_fourier_sin<double> integrator(tolerance, levels);
      auto [Q, error] = integrator.integrate(scalar_function{f}, omega);
      return nb::borrow<OouraResult>(ooura_result(Q, error));
   }, arg("f"), arg("omega"), arg("tolerance") = root_eps, arg("levels") = 8,
      "The OouraResult of the integral of f(t) sin(omega t) over [0, inf), by Ooura and Mori's double-exponential rule."
      BMPY_URL("fourier_integrals"));

   m.def("ooura_fourier_cos", [](nb::callable f, double omega, double tolerance, std::size_t levels) {
      if (levels < 1 || levels > 16)  // 0 reads an empty node table; each level doubles the work.
         throw nb::value_error("ooura_fourier_cos: levels must be in [1, 16]");
      bq::ooura_fourier_cos<double> integrator(tolerance, levels);
      auto [Q, error] = integrator.integrate(scalar_function{f}, omega);
      return nb::borrow<OouraResult>(ooura_result(Q, error));
   }, arg("f"), arg("omega"), arg("tolerance") = root_eps, arg("levels") = 8,
      "The OouraResult of the integral of f(t) cos(omega t) over [0, inf), by Ooura and Mori's double-exponential rule."
      BMPY_URL("fourier_integrals"));

   m.def("van_den_bos_unit_square", [](nb::callable f, double tolerance) {
      double error = 0, L1 = 0;
      double Q = bq::van_den_bos_unit_square(planar_function{f}, tolerance, &error, &L1);
      return result(Q, error, L1);
   }, arg("f"), arg("tolerance") = root_eps,
      "The QuadratureResult of f(x, y) over the unit square [0, 1]^2 by Van den Bos's rule." BMPY_URL("van_den_bos_unit_square"));
}

}  // namespace bmpy
