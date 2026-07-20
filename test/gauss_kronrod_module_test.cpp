// Copyright Nick Thompson, 2017
// Copyright Matt Borland, 2026
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// Module build check for boost/math/quadrature/gauss_kronrod.hpp. Building this
// against the boost.math module confirms the gauss_kronrod class template and
// its integrate() entry point are usable by a module consumer. Expected values
// are the mathematically-certain ones taken from gauss_kronrod_quadrature_test.cpp
// (a Kronrod rule integrates low-degree polynomials exactly, and 1/(1+t^2) over
// the whole real line is pi).

#ifndef BOOST_MATH_BUILD_MODULE
#include <boost/math/quadrature/gauss_kronrod.hpp>
#else
import boost.math;
#endif

#include "math_unit_test.hpp"

// The integrands below use only Real arithmetic, so no standard-library header
// is needed here: std::cout and report_errors() come from math_unit_test.hpp.

template <class Real, unsigned Points>
void test_polynomials(const char* type_name)
{
   std::cout << "Testing gauss_kronrod polynomial integration for type " << type_name << std::endl;

   using boost::math::quadrature::gauss_kronrod;

   Real error {};
   Real L1 {};

   // Linear: the integral of 5x + 7 over [0, 1] is 9.5 (exact for a Kronrod rule).
   auto f_linear = [](const Real& x) -> Real { return 5 * x + 7; };
   Real Q { gauss_kronrod<Real, Points>::integrate(f_linear, Real(0), Real(1), 0, 0, &error, &L1) };
   CHECK_ULP_CLOSE(Real(9.5), Q, 16);
   CHECK_ULP_CLOSE(Real(9.5), L1, 16);

   // Reversing the limits negates the result but leaves the L1 norm unchanged.
   Q = gauss_kronrod<Real, Points>::integrate(f_linear, Real(1), Real(0), 0, 0, &error, &L1);
   CHECK_ULP_CLOSE(Real(-9.5), Q, 16);
   CHECK_ULP_CLOSE(Real(9.5), L1, 16);

   // A degenerate interval integrates to exactly zero.
   Q = gauss_kronrod<Real, Points>::integrate(f_linear, Real(0), Real(0), 0, 0, &error, &L1);
   CHECK_EQUAL(Q, Real(0));

   // Quadratic: the integral of 5x^2 + 7x + 12 over [0, 1] is 103/6.
   auto f_quadratic = [](const Real& x) -> Real { return 5 * x * x + 7 * x + 12; };
   Q = gauss_kronrod<Real, Points>::integrate(f_quadratic, Real(0), Real(1), 0, 0, &error, &L1);
   CHECK_ULP_CLOSE(Real(103) / 6, Q, 16);
   CHECK_ULP_CLOSE(Real(103) / 6, L1, 16);

   // Exercise the default arguments of the public integrate() entry point; the
   // adaptive path also integrates this polynomial to full precision.
   Q = gauss_kronrod<Real, Points>::integrate(f_quadratic, Real(0), Real(1));
   CHECK_ULP_CLOSE(Real(103) / 6, Q, 16);
}

template <class Real, unsigned Points>
void test_infinite_limits(const char* type_name)
{
   std::cout << "Testing gauss_kronrod integration over the real line for type " << type_name << std::endl;

   using boost::math::quadrature::gauss_kronrod;

   // The integral of 1/(1+t^2) over the whole real line is pi.
   const Real pi { static_cast<Real>(3.141592653589793238462643383279502884L) };

   Real error {};
   Real L1 {};
   auto f = [](const Real& t) -> Real { return 1 / (1 + t * t); };
   Real Q { gauss_kronrod<Real, Points>::integrate(f, -boost::math::tools::max_value<Real>(),
      boost::math::tools::max_value<Real>(), 0, 0, &error, &L1) };

   // A single 15-point panel resolves this to about 6e-3 relative error.
   CHECK_ABSOLUTE_ERROR(pi, Q, Real(0.02));
   CHECK_ABSOLUTE_ERROR(pi, L1, Real(0.02));
}

int main()
{
   test_polynomials<float, 15>("float");
   test_polynomials<double, 15>("double");
   test_infinite_limits<double, 15>("double");

   return boost::math::test::report_errors();
}
