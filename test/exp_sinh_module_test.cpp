// Copyright Nick Thompson, 2017
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// Module build check for exp_sinh quadrature. The integrator is an exported
// class template whose member integrate() is instantiated on a user supplied
// functor; building this against the boost.math module confirms the exported
// template, its detail impl and the constants it pulls in are all reachable
// from a module consumer. The expected values are the closed-form results of
// the smooth half-infinite integrals used in exp_sinh_quadrature_test.cpp.

#ifndef BOOST_MATH_BUILD_MODULE
#include <boost/math/quadrature/exp_sinh.hpp>
#else
import boost.math;
#endif

#ifndef BOOST_MATH_BUILD_MODULE
#include <cmath>
#include <limits>
#endif
#include "math_unit_test.hpp"

template <class Real>
void test_spots(Real, const char* type_name)
{
   std::cout << "Testing exp_sinh quadrature for type " << type_name << std::endl;

   // root_epsilon: the default convergence tolerance the integrator requests.
   const Real tol_conv {std::sqrt(std::numeric_limits<Real>::epsilon())};
   const Real inf {std::numeric_limits<Real>::infinity()};
   const int ulps {64};

   const boost::math::quadrature::exp_sinh<Real> integrator {};

   Real error {};
   Real L1 {};

   // Integral of exp(-t)/sqrt(t) over (0, inf) is sqrt(pi). Default overload.
   const auto f1 = [](Real t) -> Real { return std::exp(-t) / std::sqrt(t); };
   Real Q {integrator.integrate(f1, tol_conv, &error, &L1)};
   CHECK_ULP_CLOSE(static_cast<Real>(1.77245385090551602729816748334114518), Q, ulps);
   // Integrand is strictly positive, so the L1 norm equals the integral.
   CHECK_ULP_CLOSE(static_cast<Real>(1.77245385090551602729816748334114518), L1, ulps);

   // Integral of 1/(1+t*t) over (0, inf) is pi/2. Default overload.
   const auto f2 = [](Real t) -> Real { return static_cast<Real>(1) / (1 + t * t); };
   Q = integrator.integrate(f2, tol_conv, &error, &L1);
   CHECK_ULP_CLOSE(static_cast<Real>(1.57079632679489661923132169163975144), Q, ulps);
   CHECK_ULP_CLOSE(static_cast<Real>(1.57079632679489661923132169163975144), L1, ulps);

   // Integral of 1/(sqrt(t)*(1+t)) over (0, inf) is pi. Default overload.
   const auto f3 = [](Real t) -> Real { return static_cast<Real>(1) / (std::sqrt(t) * (1 + t)); };
   Q = integrator.integrate(f3, tol_conv, &error, &L1);
   CHECK_ULP_CLOSE(static_cast<Real>(3.14159265358979323846264338327950288), Q, ulps);
   CHECK_ULP_CLOSE(static_cast<Real>(3.14159265358979323846264338327950288), L1, ulps);

   // Integral of 1/(1+t*t) over (1, inf) is pi/4. Bounded (right-infinite) overload.
   Q = integrator.integrate(f2, static_cast<Real>(1), inf, tol_conv, &error, &L1);
   CHECK_ULP_CLOSE(static_cast<Real>(0.785398163397448309615660845819875721), Q, ulps);
   CHECK_ULP_CLOSE(static_cast<Real>(0.785398163397448309615660845819875721), L1, ulps);

   // Integral of 1/(1+t*t) over (-inf, 0) is pi/2. Bounded (left-infinite) overload.
   Q = integrator.integrate(f2, -inf, static_cast<Real>(0), tol_conv, &error, &L1);
   CHECK_ULP_CLOSE(static_cast<Real>(1.57079632679489661923132169163975144), Q, ulps);
   CHECK_ULP_CLOSE(static_cast<Real>(1.57079632679489661923132169163975144), L1, ulps);
}

int main()
{
   test_spots(0.0F, "float");
   test_spots(0.0, "double");

   return boost::math::test::report_errors();
}
