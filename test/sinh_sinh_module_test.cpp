// Copyright Nick Thompson, 2017
// Copyright Matt Borland, 2026
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// Module build check for the sinh_sinh quadrature component. It exercises the
// exported boost::math::quadrature::sinh_sinh integrator (the constructor with a
// refinement count and the integrate() entry point, including its error and L1
// outputs) from a module consumer, using known-good integrals over the whole
// real line taken from sinh_sinh_quadrature_test.cpp.

#ifndef BOOST_MATH_BUILD_MODULE
#include <boost/math/quadrature/sinh_sinh.hpp>
#include <boost/math/constants/constants.hpp>
#else
import boost.math;
#endif

#include "math_unit_test.hpp"

#ifndef BOOST_MATH_BUILD_MODULE
#include <cmath>
#endif

template <class T>
void test_spots(T, const char* type_name)
{
   std::cout << "Testing sinh_sinh quadrature for type " << type_name << std::endl;

   using boost::math::constants::pi;
   using boost::math::constants::half_pi;
   using boost::math::constants::root_pi;
   using boost::math::quadrature::sinh_sinh;

   const sinh_sinh<T> integrator{10};
   T error {};
   T L1 {};

   // The zero integrand integrates to zero over the whole real line; this also
   // exercises the error and L1 output arguments of integrate().
   auto f0 = [](T)->T { return static_cast<T>(0); };
   T Q {integrator.integrate(f0, boost::math::tools::root_epsilon<T>(), &error, &L1)};
   CHECK_ABSOLUTE_ERROR(static_cast<T>(0), Q, static_cast<T>(1e-6));
   CHECK_ABSOLUTE_ERROR(static_cast<T>(0), L1, static_cast<T>(1e-6));

   // Integral of exp(-x^2) over the real line is sqrt(pi). The integrand is
   // positive, so the L1 norm equals the integral.
   auto f2 = [](const T& x)->T { return std::exp(-x*x); };
   Q = integrator.integrate(f2, boost::math::tools::root_epsilon<T>(), &error, &L1);
   CHECK_ULP_CLOSE(root_pi<T>(), Q, 50);
   CHECK_ULP_CLOSE(root_pi<T>(), L1, 50);

   // Integral of 1/cosh(x) over the real line is pi. This uses the single
   // argument integrate() overload.
   auto f5 = [](const T& t)->T { return static_cast<T>(1)/std::cosh(t); };
   Q = integrator.integrate(f5);
   CHECK_ULP_CLOSE(pi<T>(), Q, 50);

   // Integral of cos(x)/cosh(x) over the real line is pi/cosh(pi/2).
   auto f8 = [](const T& t)->T { return std::cos(t)/std::cosh(t); };
   const T expected {pi<T>()/std::cosh(half_pi<T>())};
   Q = integrator.integrate(f8);
   CHECK_ULP_CLOSE(expected, Q, 200);
}

int main()
{
   test_spots(0.0F, "float");
   test_spots(0.0, "double");

   return boost::math::test::report_errors();
}
