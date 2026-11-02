// Copyright 2026 Nick Thompson
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt
//
// The iterative root finders of tools/roots.hpp, tools/multivariate_roots.hpp and tools/toms748_solve.hpp on a
// Python callable, and the closed-form quadratic, cubic and quartic solvers. The callable is called with the GIL
// held, and an exception it raises propagates. Boost.Math errors (no sign change, roots out of order, ...) follow
// boost.math.errstate as for a scalar. Tolerances are a number of binary digits, as in Boost.

#include "bindings.hpp"

#include <boost/math/tools/cubic_roots.hpp>
#include <boost/math/tools/multivariate_roots.hpp>
#include <boost/math/tools/quartic_roots.hpp>
#include <boost/math/tools/roots.hpp>
#include <boost/math/tools/toms748_solve.hpp>

#include <nanobind/stl/complex.h>
#include <nanobind/stl/pair.h>
#include <nanobind/stl/tuple.h>

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdint>
#include <limits>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace bmpy {

namespace {

namespace bt = boost::math::tools;
using C = std::complex<double>;

constexpr int default_digits = std::numeric_limits<double>::digits;
constexpr std::uintmax_t default_max_iterations = 1000;

// fn() with the Boost.Math errors it raises handled as for a scalar argument; see errors.hpp.
template <class Fn>
auto solve(char const* name, Fn fn)
{
   error_actions actions = current_error_actions;
   error_record errors;
   auto r = [&] {
      recording rec(errors);
      return fn();
   }();
   if (errors.any())
      report(name, errors, 1, true, false, actions);
   return r;
}

bt::eps_tolerance<double> tolerance(int digits)
{
   if (digits < 1)
      throw nb::value_error("digits must be positive");
   return bt::eps_tolerance<double>(static_cast<unsigned>(digits));
}

// f(x) -> (F, J) or (F, J, d2) for the multivariate solvers: F and d2(v) = F''(x)[v, v] of shape (n,), J of shape
// (n, n), passed on row-major.
struct py_system
{
   nb::callable f;
   std::size_t n;

   auto operator()(std::vector<double> const& x) const
   {
      recording suspended(nullptr);
      nb::object r = f(vector_to_array(x));
      std::size_t const n_ = n;
      // a as a float64 array of shape (n,), or (n, n) when square, flattened row-major.
      auto shaped = [n_](nb::handle a, bool square, char const* what) {
         nb::object arr = numpy().attr("asarray")(a, nb::arg("dtype") = numpy().attr("float64"));
         nb::object shape = arr.attr("shape");
         bool ok = nb::len(shape) == (square ? 2u : 1u);
         for (std::size_t i = 0; ok && i < nb::len(shape); ++i)
            ok = nb::cast<std::size_t>(shape[nb::int_(i)]) == n_;
         if (!ok)
         {
            std::string got = nb::str(shape).c_str();
            throw nb::value_error((std::string(what) + " has shape " + got + ", expected (" + std::to_string(n_)
                                   + (square ? ", " + std::to_string(n_) : std::string(",")) + ")").c_str());
         }
         return array_to_vector(arr);
      };
      nb::object d2 = nb::len(r) > 2 ? nb::object(r[2]) : nb::none();
      auto second = [d2, shaped](std::vector<double> const& v) {
         recording inner(nullptr);
         if (d2.is_none())
            throw nb::value_error("multivariate_halley: f(x) must return (F, J, d2)");
         return shaped(d2(vector_to_array(v)), false, "d2(v)");
      };
      return std::make_tuple(shaped(r[0], false, "F(x)"), shaped(r[1], true, "J(x)"), second);
   }
};

// The real roots, ascending; Boost pads its fixed-size result with NaN.
template <std::size_t N>
auto real_roots(std::array<double, N> const& roots)
{
   std::vector<double> v;
   std::copy_if(roots.begin(), roots.end(), std::back_inserter(v), [](double r) { return !std::isnan(r); });
   return vector_to_array(std::move(v));
}

}  // namespace

void bind_roots(nb::module_& m)
{
   using nb::arg;
   nb::object root_result = named_tuple(m, "RootResult", "root lower upper iterations",
      "The midpoint of the final bracket [lower, upper], and the iterations used (max_iterations: not converged).");
   nb::object newton_result = named_tuple(m, "NewtonResult", "root iterations",
      "The root, and the iterations used (max_iterations: not converged).");

   auto bracket = [root_result](std::pair<double, double> ab, std::uintmax_t iterations) {
      return root_result(ab.first + (ab.second - ab.first) / 2, ab.first, ab.second, iterations);
   };

   m.def("bisect", [bracket](nb::callable f, double a, double b, int digits, std::uintmax_t max_iterations) {
      auto tol = tolerance(digits);
      auto ab = solve("bisect", [&] { return bt::bisect(py_callable<double, double>{f}, a, b, tol, max_iterations); });
      return bracket(ab, max_iterations);
   }, arg(), arg(), arg(), nb::kw_only(), arg("digits") = default_digits, arg("max_iterations") = default_max_iterations,
      nb::sig("def bisect(f, a, b, /, *, digits=53, max_iterations=1000) -> RootResult"),
      "A root of f in [a, b] by bisection; f(a) and f(b) must differ in sign." BMPY_URL("roots_noderiv/bisect"));

   m.def("toms748", [bracket](nb::callable f, double a, double b, int digits, std::uintmax_t max_iterations) {
      auto tol = tolerance(digits);
      auto ab = solve("toms748", [&] { return bt::toms748_solve(py_callable<double, double>{f}, a, b, tol, max_iterations); });
      return bracket(ab, max_iterations);
   }, arg(), arg(), arg(), nb::kw_only(), arg("digits") = default_digits, arg("max_iterations") = default_max_iterations,
      nb::sig("def toms748(f, a, b, /, *, digits=53, max_iterations=1000) -> RootResult"),
      "A root of f in [a, b] by TOMS Algorithm 748; f(a) and f(b) must differ in sign." BMPY_URL("roots_noderiv/TOMS748"));

   m.def("bracket_and_solve_root", [bracket](nb::callable f, double guess, double factor, bool rising, int digits,
                                             std::uintmax_t max_iterations) {
      auto tol = tolerance(digits);
      auto ab = solve("bracket_and_solve_root", [&] {
         return bt::bracket_and_solve_root(py_callable<double, double>{f}, guess, factor, rising, tol, max_iterations);
      });
      return bracket(ab, max_iterations);
   }, arg(), arg(), arg(), arg(), nb::kw_only(), arg("digits") = default_digits, arg("max_iterations") = default_max_iterations,
      nb::sig("def bracket_and_solve_root(f, guess, factor, rising, /, *, digits=53, max_iterations=1000) -> RootResult"),
      "A root of a monotone f, increasing if rising: brackets it by scaling guess (not 0) by factor, then toms748."
      BMPY_URL("roots_noderiv/bracket_solve"));

   auto def_iterate = [&](char const* name, auto iterate, auto function, char const* sig, char const* doc) {
      m.def(name, [newton_result, name, iterate, function](nb::callable f, double guess, double min, double max, int digits,
                                                         std::uintmax_t max_iterations) {
         double root = solve(name, [&] { return iterate(decltype(function){f}, guess, min, max, digits, max_iterations); });
         return newton_result(root, max_iterations);
      }, arg(), arg(), arg(), arg(), nb::kw_only(), arg("digits") = default_digits, arg("max_iterations") = default_max_iterations,
         nb::sig(sig), doc);
   };
   using first = py_callable<std::tuple<double, double>, double>;
   using second = py_callable<std::tuple<double, double, double>, double>;
   def_iterate("newton", [](auto&&... a) { return bt::newton_raphson_iterate(a...); }, first{},
      "def newton(f, guess, min, max, /, *, digits=53, max_iterations=1000) -> NewtonResult",
      "A root of f in [min, max] by Newton-Raphson iteration from guess; f(x) returns (f(x), f'(x))." BMPY_URL("roots_deriv"));
   def_iterate("halley", [](auto&&... a) { return bt::halley_iterate(a...); }, second{},
      "def halley(f, guess, min, max, /, *, digits=53, max_iterations=1000) -> NewtonResult",
      "A root of f in [min, max] by Halley iteration from guess; f(x) returns (f(x), f'(x), f''(x))." BMPY_URL("roots_deriv"));
   def_iterate("schroder", [](auto&&... a) { return bt::schroder_iterate(a...); }, second{},
      "def schroder(f, guess, min, max, /, *, digits=53, max_iterations=1000) -> NewtonResult",
      "A root of f in [min, max] by Schroder iteration from guess; f(x) returns (f(x), f'(x), f''(x))." BMPY_URL("roots_deriv"));

   m.def("complex_newton", [](nb::callable f, C guess, int max_iterations) {
      if (max_iterations < 0)
         throw nb::value_error("complex_newton: max_iterations must be non-negative");
      return bt::complex_newton(py_callable<std::pair<C, C>, C>{f}, guess, max_iterations);
   }, arg(), arg(), nb::kw_only(), arg("max_iterations") = default_digits,
      nb::sig("def complex_newton(f, guess, /, *, max_iterations=53) -> complex"),
      "A root of the analytic f near guess by Newton's method; f(z) returns (f(z), f'(z)); nan+nanj if it fails."
      BMPY_URL("roots_deriv"));

   auto def_multivariate = [&](char const* name, auto iterate, char const* sig, char const* doc) {
      m.def(name, [newton_result, name, iterate](nb::callable f, nb::handle x0, int digits, std::uintmax_t max_iterations) {
         std::vector<double> x = array_to_vector(x0);
         if (x.empty())
            throw nb::value_error((std::string(name) + ": x0 must not be empty").c_str());
         py_system system{f, x.size()};
         auto root = solve(name, [&] { return iterate(system, x, digits, max_iterations); });
         return newton_result(vector_to_array(std::move(root)), max_iterations);
      }, arg(), arg(), nb::kw_only(), arg("digits") = default_digits, arg("max_iterations") = default_max_iterations,
         nb::sig(sig), doc);
   };
   def_multivariate("multivariate_newton", [](auto&&... a) { return bt::multivariate_newton(a...); },
      "def multivariate_newton(f, x0, /, *, digits=53, max_iterations=1000) -> NewtonResult",
      "A root of the system F(x) = 0 near x0 by Newton's method with a line search; f(x) returns (F(x), J(x)), "
      "arrays of shape (n,) and (n, n), J[i, j] = dF_i/dx_j." BMPY_URL("multivariate_roots"));
   def_multivariate("multivariate_halley", [](auto&&... a) { return bt::multivariate_halley(a...); },
      "def multivariate_halley(f, x0, /, *, digits=53, max_iterations=1000) -> NewtonResult",
      "A root of the system F(x) = 0 near x0 by Cuyt's multivariate Halley method; f(x) returns (F(x), J(x), d2), "
      "where d2(v) is F''(x)[v, v], the second derivative of F(x + t v) in t at 0." BMPY_URL("multivariate_roots"));

   m.def("quadratic_roots", [](double a, double b, double c) {
      auto [x0, x1] = bt::quadratic_roots(a, b, c);
      return real_roots(std::array<double, 2>{x0, x1});
   }, arg(), arg(), arg(), nb::sig("def quadratic_roots(a, b, c, /)"),
      "The real roots of a x^2 + b x + c, ascending, as a float64 array (repeated roots repeated)." BMPY_URL("roots_deriv"));

   m.def("cubic_roots", [](double a, double b, double c, double d) {
      return real_roots(bt::cubic_roots(a, b, c, d));
   }, arg(), arg(), arg(), arg(), nb::sig("def cubic_roots(a, b, c, d, /)"),
      "The real roots of a x^3 + b x^2 + c x + d, ascending, as a float64 array; highest degree first, unlike "
      "polynomial_roots." BMPY_URL("cubic_roots"));

   m.def("quartic_roots", [](double a, double b, double c, double d, double e) {
      return real_roots(bt::quartic_roots(a, b, c, d, e));
   }, arg(), arg(), arg(), arg(), arg(), nb::sig("def quartic_roots(a, b, c, d, e, /)"),
      "The real roots of a x^4 + b x^3 + c x^2 + d x + e, ascending, as a float64 array; highest degree first, unlike "
      "polynomial_roots." BMPY_URL("quartic_roots"));
}

}  // namespace bmpy
