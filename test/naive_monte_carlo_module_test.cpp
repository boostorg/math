/*
 * Copyright Nick Thompson, 2018
 * Copyright Matt Borland, 2026
 * Use, modification and distribution are subject to the
 * Boost Software License, Version 1.0. (See accompanying file
 * LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
 */
//
// Module build check for boost::math::quadrature::naive_monte_carlo. When
// BOOST_MATH_BUILD_MODULE is defined the component is consumed through
// `import boost.math;` and the standard library through `import std;` (performed
// by math_unit_test.hpp). Nothing standard or Boost is included textually in
// that mode. Expected values are the exact analytic integrals of the constant
// and linear integrands taken from naive_monte_carlo_test.cpp.

#ifndef BOOST_MATH_BUILD_MODULE
#include <boost/math/quadrature/naive_monte_carlo.hpp>
#else
import boost.math;
#endif

#ifndef BOOST_MATH_BUILD_MODULE
#include <vector>
#include <utility>
#endif
#include "math_unit_test.hpp"

using boost::math::quadrature::naive_monte_carlo;

template <class Real>
void test_spots(Real, const char* type_name)
{
    std::cout << "Testing naive_monte_carlo for type " << type_name << std::endl;

    // The integral of the constant 1 over the unit square is exactly 1. A
    // constant integrand has zero variance, so the estimate and the error
    // estimate are exact regardless of how many samples are drawn.
    {
        auto g = [](std::vector<Real> const&) -> Real
        {
            return Real(1);
        };
        std::vector<std::pair<Real, Real>> bounds{{Real(0), Real(1)}, {Real(0), Real(1)}};
        naive_monte_carlo<Real, decltype(g)> mc(g, bounds, Real(0.0001), false, 1, 87);

        auto task = mc.integrate();
        Real one = task.get();
        CHECK_ABSOLUTE_ERROR(Real(1), one, Real(0.001));
        CHECK_LE(mc.current_error_estimate(), Real(1e-6));
        CHECK_GE(static_cast<double>(mc.calls()), 1000.0);
    }

    // The integral of x over [0, 1] is 1/2 and the variance of a uniform
    // variate on [0, 1] is 1/12.
    {
        auto g = [](std::vector<Real> const& x) -> Real
        {
            return x[0];
        };
        std::vector<std::pair<Real, Real>> bounds{{Real(0), Real(1)}};
        naive_monte_carlo<Real, decltype(g)> mc(g, bounds, Real(0.001), false, 1, 12341);

        auto task = mc.integrate();
        Real y = task.get();
        CHECK_ABSOLUTE_ERROR(Real(0.5), y, Real(0.01));
        CHECK_ABSOLUTE_ERROR(Real(1) / Real(12), mc.variance(), Real(0.01));
    }
}

int main()
{
    test_spots(0.0, "double");
    test_spots(0.0F, "float");

    return boost::math::test::report_errors();
}
