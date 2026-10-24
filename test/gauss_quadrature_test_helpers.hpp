// Copyright Jacob Hass, 2026
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// Checks shared by the Gauss-Hermite, Gauss-Laguerre and Gauss-Jacobi tests.

#ifndef BOOST_MATH_TEST_GAUSS_QUADRATURE_TEST_HELPERS_HPP
#define BOOST_MATH_TEST_GAUSS_QUADRATURE_TEST_HELPERS_HPP

#include "math_unit_test.hpp"
#include <cmath>
#include <limits>
#include <valarray>

// The integrals of 1, x^2 and cos(x) against the rule's weight function, as one valarray.
template <class Rule>
void check_vector_valued(const Rule& rule, double m0, double m2, double cos_moment, unsigned ulps0, unsigned ulps2)
{
   using vector = std::valarray<double>;
   const vector zero(0.0, 3);
   auto norm = [](const vector& v) { return std::sqrt((v * v).sum()); };
   auto f = [](double x) { return vector{1.0, x * x, std::cos(x)}; };
   double L1;
   vector Q = rule.integrate(f, zero, norm, &L1);
   CHECK_ULP_CLOSE(m0, Q[0], ulps0);
   CHECK_ULP_CLOSE(m2, Q[1], ulps2);
   CHECK_ULP_CLOSE(cos_moment, Q[2], 16);
   CHECK_LE(norm(Q), L1);
   // The error is the norm of half the difference from the anti-Gauss rule, at least that of any component.
   double error, L1_averaged, cos_error;
   const vector Q_averaged = rule.integrate_averaged(f, zero, norm, &error, &L1_averaged);
   CHECK_ULP_CLOSE(m0, Q_averaged[0], ulps0);
   CHECK_ULP_CLOSE(m2, Q_averaged[1], ulps2);
   CHECK_ABSOLUTE_ERROR(cos_moment, Q_averaged[2], error);
   CHECK_LE(norm(Q_averaged), L1_averaged);
   rule.integrate_averaged([](double x) { return std::cos(x); }, &cos_error);
   CHECK_LE(cos_error, error);
   CHECK_LE(2 * std::numeric_limits<double>::epsilon() * L1_averaged, error);
}

// For smooth f, G and A bracket the integral, so the error bounds that of the average and is close to that of G
// until it reaches rounding. The error may fall short of that of the average by slack eps L1.
template <class Real, class Rule, class F>
void check_error_estimate(const Rule& rule, F f, const Real& exact, const Real& slack = 0)
{
   using std::abs;
   Real error, L1;
   const Real Q = rule.integrate_averaged(f, &error, &L1);
   const Real gauss_error = abs(rule.integrate(f) - exact);
   CHECK_LE(abs(Q), L1);
   CHECK_LE(abs(Q - exact), error + slack * std::numeric_limits<Real>::epsilon() * L1);
   CHECK_LE(abs(Q - exact), gauss_error + 2 * std::numeric_limits<Real>::epsilon() * L1);
   if (gauss_error > 1000 * std::numeric_limits<Real>::epsilon() * L1)
   {
      CHECK_LE(error, 2 * gauss_error);
      CHECK_LE(gauss_error, 2 * error);
   }
   CHECK_LE(2 * std::numeric_limits<Real>::epsilon() * L1, error);
}

// The integral of a 2x2 complex matrix polynomial, with the checks every rule shares; the callers check the moments.
// Return concrete Eigen matrices so the test does not rely on the lifetime
// of expression-template temporaries. Supply a shaped zero and a scalar norm.
template <class Matrix, class Rule>
Matrix integrate_matrix(const Rule& rule, const Matrix& zero)
{
   const auto norm = [](const Matrix& m) { return m.norm(); };
   const auto f = [&zero](double x) -> Matrix {
      Matrix result = zero;
      result(0, 0) = {1, x * x};
      result(0, 1) = {x, 2 * x};
      result(1, 0) = {-x * x, x};
      result(1, 1) = {x * x * x, -1};
      return result;
   };
   double L1;
   Matrix Q = rule.integrate(f, zero, norm, &L1);
   // The integrand is a polynomial that both rules integrate exactly, so the estimate is rounding.
   double error;
   rule.integrate_averaged(f, zero, norm, &error);
   CHECK_LE(error, 64 * std::numeric_limits<double>::epsilon() * L1);
   CHECK_EQUAL(Q.rows(), decltype(Q.rows())(2));
   CHECK_EQUAL(Q.cols(), decltype(Q.cols())(2));
   CHECK_LE(norm(Q), L1);
   return Q;
}

#endif // BOOST_MATH_TEST_GAUSS_QUADRATURE_TEST_HELPERS_HPP
