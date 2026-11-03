// Copyright Nick Thompson 2026.
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0. (See accompanying file
// LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include "math_unit_test.hpp"
#include <boost/math/tools/multivariate_roots.hpp>
#include <boost/multiprecision/cpp_bin_float.hpp>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <tuple>
#include <vector>

using boost::math::tools::multivariate_halley;
using boost::math::tools::multivariate_newton;

// Each system returns (F, J, F''[v, v]); Newton ignores the third element.

// x^2 + y^2 = 4 and x = 2y, with roots +-(4, 2)/sqrt(5). J is a std::array of rows.
template <class Real>
auto circle_line(std::array<Real, 2> const& x)
{
   std::array<Real, 2> F{x[0] * x[0] + x[1] * x[1] - 4, x[0] - 2 * x[1]};
   std::array<std::array<Real, 2>, 2> J{{{2 * x[0], 2 * x[1]}, {1, -2}}};
   auto d2 = [](std::array<Real, 2> const& v) { return std::array<Real, 2>{2 * (v[0] * v[0] + v[1] * v[1]), 0}; };
   return std::make_tuple(F, J, d2);
}

// The gradient of the Rosenbrock function (1 - x)^2 + 100(y - x^2)^2, zero only at (1, 1).
template <class Real>
auto rosenbrock_gradient(std::vector<Real> const& x)
{
   Real const u = x[0];
   Real const w = x[1];
   std::vector<Real> F{-2 * (1 - u) - 400 * u * (w - u * u), 200 * (w - u * u)};
   std::vector<std::vector<Real>> J{{2 - 400 * (w - u * u) + 800 * u * u, -400 * u}, {-400 * u, Real(200)}};
   auto d2 = [u](std::vector<Real> const& v) {
      return std::vector<Real>{-800 * v[0] * v[1] + 2400 * u * v[0] * v[0], -400 * v[0] * v[0]};
   };
   return std::make_tuple(F, J, d2);
}

// Broyden's tridiagonal function (More, Garbow and Hillstrom 1981, #30); J is row-major in a std::vector.
template <class Real>
auto broyden_tridiagonal(std::vector<Real> const& x)
{
   std::size_t const n = x.size();
   std::vector<Real> F(n);
   std::vector<Real> J(n * n, Real(0));
   for (std::size_t i = 0; i < n; ++i)
   {
      Real const left = i > 0 ? x[i - 1] : Real(0);
      Real const right = i + 1 < n ? x[i + 1] : Real(0);
      F[i] = (3 - 2 * x[i]) * x[i] - left - 2 * right + 1;
      J[i * n + i] = 3 - 4 * x[i];
      if (i > 0)
         J[i * n + i - 1] = -1;
      if (i + 1 < n)
         J[i * n + i + 1] = -2;
   }
   auto d2 = [](std::vector<Real> const& v) {
      std::vector<Real> w(v.size());
      for (std::size_t i = 0; i < v.size(); ++i)
         w[i] = -4 * v[i] * v[i];
      return w;
   };
   return std::make_tuple(F, J, d2);
}

// Powell's badly scaled function (More, Garbow and Hillstrom 1981, #3), root near (1.098e-5, 9.106).
template <class Real>
auto powell_badly_scaled(std::vector<Real> const& x)
{
   using std::exp;
   Real const ex = exp(-x[0]);
   Real const ey = exp(-x[1]);
   std::vector<Real> F{10000 * x[0] * x[1] - 1, ex + ey - Real(1.0001)};
   std::vector<std::vector<Real>> J{{10000 * x[1], 10000 * x[0]}, {-ex, -ey}};
   auto d2 = [ex, ey](std::vector<Real> const& v) {
      return std::vector<Real>{20000 * v[0] * v[1], ex * v[0] * v[0] + ey * v[1] * v[1]};
   };
   return std::make_tuple(F, J, d2);
}

// x^2 = 1, y^2 = 1: J is singular wherever x = 0 or y = 0.
template <class Real>
auto squares(std::vector<Real> const& x)
{
   std::vector<Real> F{x[0] * x[0] - 1, x[1] * x[1] - 1};
   std::vector<Real> J{2 * x[0], 0, 0, 2 * x[1]};
   auto d2 = [](std::vector<Real> const& v) { return std::vector<Real>{2 * v[0] * v[0], 2 * v[1] * v[1]}; };
   return std::make_tuple(F, J, d2);
}

template <class Real>
Real residual(std::vector<Real> const& F)
{
   using std::abs;
   Real r = 0;
   for (Real f : F)
      r = (std::max)(r, Real(abs(f)));
   return r;
}

template <class Real>
void test_circle_line()
{
   int const digits = std::numeric_limits<Real>::digits;
   Real const s = 1 / std::sqrt(Real(5));
   for (bool halley : {false, true})
   {
      std::uintmax_t iters = 100;
      std::array<Real, 2> x0{1, 1};
      auto f = [](std::array<Real, 2> const& x) { return circle_line(x); };
      auto x = halley ? multivariate_halley(f, x0, digits, iters) : multivariate_newton(f, x0, digits, iters);
      CHECK_ULP_CLOSE(4 * s, x[0], 2);
      CHECK_ULP_CLOSE(2 * s, x[1], 2);
      CHECK_LE(iters, std::uintmax_t(8));
      // From the other side, the other root:
      x = halley ? multivariate_halley(f, std::array<Real, 2>{-3, 1}, digits) : multivariate_newton(f, std::array<Real, 2>{-3, 1}, digits);
      CHECK_ULP_CLOSE(-4 * s, x[0], 2);
      CHECK_ULP_CLOSE(-2 * s, x[1], 2);
   }
}

template <class Real>
void test_rosenbrock()
{
   int const digits = std::numeric_limits<Real>::digits;
   auto f = [](std::vector<Real> const& x) { return rosenbrock_gradient(x); };
   std::uintmax_t newton_iters = 1000;
   std::uintmax_t halley_iters = 1000;
   auto xn = multivariate_newton(f, std::vector<Real>{-1.2, 1}, digits, newton_iters);
   auto xh = multivariate_halley(f, std::vector<Real>{-1.2, 1}, digits, halley_iters);
   for (auto const& x : {xn, xh})
   {
      CHECK_ULP_CLOSE(Real(1), x[0], 4);
      CHECK_ULP_CLOSE(Real(1), x[1], 8);
   }
   // The line search makes Newton creep along the curved valley; Halley's steps bend with it.
   CHECK_LE(newton_iters, std::uintmax_t(300));
   CHECK_LE(halley_iters, std::uintmax_t(100));
}

template <class Real>
void test_broyden_tridiagonal()
{
   int const digits = std::numeric_limits<Real>::digits;
   auto f = [](std::vector<Real> const& x) { return broyden_tridiagonal(x); };
   std::uintmax_t newton_iters = 100;
   std::uintmax_t halley_iters = 100;
   auto xn = multivariate_newton(f, std::vector<Real>(10, -1), digits, newton_iters);
   auto xh = multivariate_halley(f, std::vector<Real>(10, -1), digits, halley_iters);
   Real const tol = 8 * std::numeric_limits<Real>::epsilon();
   CHECK_LE(residual(std::get<0>(broyden_tridiagonal(xn))), tol);
   CHECK_LE(residual(std::get<0>(broyden_tridiagonal(xh))), tol);
   for (std::size_t i = 0; i < 10; ++i)
      CHECK_ULP_CLOSE(xn[i], xh[i], 8);
   CHECK_LE(newton_iters, std::uintmax_t(8));
   // In float the last iterations are at the rounding level, and their number is noise.
   CHECK_LE(halley_iters, newton_iters + (std::is_same<Real, float>::value ? 2 : 0));
}

template <class Real>
void test_powell_badly_scaled()
{
   int const digits = std::numeric_limits<Real>::digits;
   auto f = [](std::vector<Real> const& x) { return powell_badly_scaled(x); };
   for (bool halley : {false, true})
   {
      std::uintmax_t iters = 200;
      std::vector<Real> x0{0, 1};
      auto x = halley ? multivariate_halley(f, x0, digits, iters) : multivariate_newton(f, x0, digits, iters);
      CHECK_LE(iters, std::uintmax_t(200 - 1));
      // 1.0001 is not representable, so compare with the root of the system as rounded:
      CHECK_LE(residual(std::get<0>(powell_badly_scaled(x))), 8 * std::numeric_limits<Real>::epsilon());
      if (std::is_same<Real, double>::value)
      {
         CHECK_RELATIVE_ERROR(Real(1.0981593296997149e-05L), x[0], Real(1e-12));
         CHECK_RELATIVE_ERROR(Real(9.1061467398665240L), x[1], Real(1e-12));
      }
   }
}

template <class Real>
void test_singular()
{
   int const digits = std::numeric_limits<Real>::digits;
   auto f = [](std::vector<Real> const& x) { return squares(x); };
   // J = 0 at the start:
   CHECK_THROW(multivariate_newton(f, std::vector<Real>{0, 0}, digits), boost::math::evaluation_error);
   CHECK_THROW(multivariate_halley(f, std::vector<Real>{0, 0.5}, digits), boost::math::evaluation_error);
   // Nearly singular: the Newton step is huge, and the line search pulls it back.
   for (bool halley : {false, true})
   {
      std::vector<Real> x0{Real(1) / 1024, Real(0.5)};
      auto x = halley ? multivariate_halley(f, x0, digits) : multivariate_newton(f, x0, digits);
      CHECK_ULP_CLOSE(Real(1), x[0], 1);
      CHECK_ULP_CLOSE(Real(1), x[1], 1);
   }
   // x^2 + 1 = 0 has no real root; |F| is least at 0, where J is singular.
   auto g = [](std::array<Real, 1> const& x) {
      return std::make_pair(std::array<Real, 1>{x[0] * x[0] + 1}, std::array<Real, 1>{2 * x[0]});
   };
   CHECK_THROW(multivariate_newton(g, std::array<Real, 1>{1}, digits), boost::math::evaluation_error);
}

// The max-norm errors of the iterates from (2.5, 0.5) on circle_line, the 0th to the 7th.
template <class Real>
std::vector<Real> circle_line_errors(bool halley)
{
   using std::abs;
   using std::sqrt;
   auto f = [](std::array<Real, 2> const& x) { return circle_line(x); };
   Real const s = 1 / sqrt(Real(5));
   std::vector<Real> e;
   for (std::uintmax_t k = 0; k < 8; ++k)
   {
      std::uintmax_t iters = k;
      std::array<Real, 2> x0{Real(2.5), Real(0.5)};
      auto x = halley ? multivariate_halley(f, x0, std::numeric_limits<Real>::digits, iters)
                      : multivariate_newton(f, x0, std::numeric_limits<Real>::digits, iters);
      e.push_back((std::max)(abs(x[0] - 4 * s), abs(x[1] - 2 * s)));
   }
   return e;
}

// Near the root e_{k+1} <= e_k^2 for Newton and e_k^3 for Halley (the constants are about 0.3 and 0.5), as long
// as the error is above rounding.
template <class Real>
void test_convergence_order()
{
   for (bool halley : {false, true})
   {
      auto e = circle_line_errors<Real>(halley);
      int const p = halley ? 3 : 2;
      int checked = 0;
      for (std::size_t k = 0; k + 1 < e.size(); ++k)
      {
         if (e[k] < Real(0.1) && e[k + 1] > 10 * std::numeric_limits<Real>::epsilon())
         {
            CHECK_LE(e[k + 1], Real(std::pow(e[k], p)));
            ++checked;
         }
      }
      CHECK_GE(checked, 1);
   }
}

// The order itself, log(e_{k+2}/e_{k+1}) / log(e_{k+1}/e_k), in 50 digits, from the last three errors well above
// rounding.
void test_convergence_order_multiprecision()
{
   using Real = boost::multiprecision::cpp_bin_float_50;
   for (bool halley : {false, true})
   {
      auto e = circle_line_errors<Real>(halley);
      Real const floor = 1e-40;
      std::size_t k = 0;
      while (k + 3 < e.size() && e[k + 3] > floor)
         ++k;
      Real const order = log(e[k + 2] / e[k + 1]) / log(e[k + 1] / e[k]);
      CHECK_LE(abs(order - (halley ? 3 : 2)), Real(0.1));
   }
}

// At a double or triple root Newton and Halley converge only linearly, with steps shrinking by a constant factor; they
// must not be mistaken for rounding noise. Here F = (x - 1)^m is computed exactly near the root, so the iteration
// can go all the way (with rounding errors in F, only about digits/m digits would be attainable).
template <class Real>
void test_multiple_roots()
{
   using std::abs;
   int const digits = boost::math::tools::digits<Real>();
   for (int m : {2, 3})
   {
      auto f = [m](std::vector<Real> const& x) {
         Real const d = x[0] - 1;
         Real const d2m = m == 2 ? Real(1) : d;  // d^(m-2)
         std::vector<Real> F{d * d * d2m, x[1] - 2};
         std::vector<Real> J{m * d * d2m, 0, 0, 1};
         auto d2 = [m, d2m](std::vector<Real> const& v) { return std::vector<Real>{m * (m - 1) * d2m * v[0] * v[0], 0}; };
         return std::make_tuple(F, J, d2);
      };
      for (bool halley : {false, true})
      {
         for (Real start : {Real(1.001), Real(1.5), Real(3)})
         {
            std::uintmax_t iters = 1000;
            std::vector<Real> x0{start, 0};
            auto x = halley ? multivariate_halley(f, x0, digits, iters) : multivariate_newton(f, x0, digits, iters);
            CHECK_LE(abs(x[0] - 1), 8 * boost::math::tools::epsilon<Real>());
            CHECK_ULP_CLOSE(Real(2), x[1], 1);
            CHECK_LE(iters, std::uintmax_t(1000 - 1));
         }
      }
   }
}

// cond(J) ~ 4e6, and F carries an erratic error of 1e-13 |x|, as a long computation of F might: near the root (1000,
// 1000) the Newton step is that error amplified by cond(J), about 1e-7 |x|. It is too big to be ignored and cannot
// decrease |F|, but x is as accurate as F allows, so that is convergence, not a failure.
void test_noisy_ill_conditioned()
{
   using std::abs;
   using std::sin;
   double const delta = 1e-6;
   double const c = 1000;
   auto f = [delta, c](std::vector<double> const& x) {
      double const u = x[0] - c;
      double const w = x[1] - c;
      double const p = u + w;
      std::vector<double> F{p + p * p / (10 * c), u + (1 + delta) * w};
      for (std::size_t i = 0; i < 2; ++i)
         F[i] += 1e-13 * c * sin(1e6 * (x[0] + (2 + i) * x[1]));
      std::vector<std::vector<double>> J{{1 + p / (5 * c), 1 + p / (5 * c)}, {1, 1 + delta}};
      auto d2 = [c](std::vector<double> const& v) {
         double const s = v[0] + v[1];
         return std::vector<double>{s * s / (5 * c), 0};
      };
      return std::make_tuple(F, J, d2);
   };
   for (bool halley : {false, true})
   {
      for (auto const& x0 : std::vector<std::vector<double>>{{500, 1750}, {900, 1550}, {1250, 1375}, {2000, 1000}, {1100, 1000}})
      {
         std::vector<double> x{0, 0};
         CHECK_NO_THROW(x = halley ? multivariate_halley(f, x0, 53) : multivariate_newton(f, x0, 53));
         CHECK_LE(abs(x[0] / c - 1), 1e-6);
         CHECK_LE(abs(x[1] / c - 1), 1e-6);
      }
   }
}

// |F|^2 underflows long before F does; only F = 0 exactly stops the iteration.
void test_tiny_residual()
{
   auto f = [](std::array<double, 2> const& x) {
      auto r = circle_line(x);
      auto F = std::get<0>(r);
      auto J = std::get<1>(r);
      for (std::size_t i = 0; i < 2; ++i)
      {
         F[i] *= 1e-160;
         for (std::size_t j = 0; j < 2; ++j)
            J[i][j] *= 1e-160;
      }
      return std::make_pair(F, J);
   };
   auto x = multivariate_newton(f, std::array<double, 2>{1, 1}, 53);
   CHECK_ULP_CLOSE(4 / std::sqrt(5.0), x[0], 2);
   CHECK_ULP_CLOSE(2 / std::sqrt(5.0), x[1], 2);
}

void test_bad_arguments()
{
   auto f = [](std::vector<double> const& x) {
      auto d2 = [](std::vector<double> const&) { return std::vector<double>{0.0}; };
      return std::make_tuple(std::vector<double>{x[0] * x[0] - 2, x[1] - 1}, std::vector<double>{2 * x[0], 0, 0, 1}, d2);
   };
   CHECK_THROW(multivariate_halley(f, std::vector<double>{3, 0}, 53), boost::math::evaluation_error);
   auto g = [](std::vector<double> const& x) { return std::make_pair(x, x); };
   CHECK_THROW(multivariate_newton(g, std::vector<double>{}, 53), std::domain_error);
   CHECK_THROW(multivariate_newton(g, std::vector<double>{1, 2}, 53), boost::math::evaluation_error);
}

int main()
{
   test_circle_line<float>();
   test_circle_line<double>();
   test_rosenbrock<float>();
   test_rosenbrock<double>();
   test_broyden_tridiagonal<float>();
   test_broyden_tridiagonal<double>();
   test_powell_badly_scaled<float>();
   test_powell_badly_scaled<double>();
   test_singular<float>();
   test_singular<double>();
   test_convergence_order<float>();
   test_convergence_order<double>();
   test_convergence_order_multiprecision();
   test_multiple_roots<float>();
   test_multiple_roots<double>();
   test_noisy_ill_conditioned();
   test_tiny_residual();
   test_bad_arguments();
   return boost::math::test::report_errors();
}
