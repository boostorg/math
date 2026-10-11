/*
 * Copyright Nick Thompson, 2017
 * Use, modification and distribution are subject to the
 * Boost Software License, Version 1.0. (See accompanying file
 * LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
 *
 * Module build check for boost::math::quadrature::trapezoidal. It exercises the
 * exported adaptive trapezoidal entry point when the test is built against the
 * boost.math module. The double/float expected values mirror the corresponding
 * assertions in test_trapezoidal.cpp.
 */

#ifndef BOOST_MATH_BUILD_MODULE
#include <boost/math/quadrature/trapezoidal.hpp>
#include <boost/math/constants/constants.hpp>
#else
import boost.math;
#endif

#ifndef BOOST_MATH_BUILD_MODULE
#include <cmath>
#include <limits>
#endif
#include "math_unit_test.hpp"

// Constant and linear integrands: the trapezoidal rule reproduces these exactly.
template <class Real>
void test_polynomial()
{
    using boost::math::quadrature::trapezoidal;

    // Integral of the constant 1/2 over [0, 10] is 5.
    auto c = [](Real)->Real { return static_cast<Real>(0.5); };
    Real Q = trapezoidal(c, static_cast<Real>(0), static_cast<Real>(10));
    CHECK_ULP_CLOSE(static_cast<Real>(5), Q, 8);

    // Reversing the limits of integration negates the result.
    Q = trapezoidal(c, static_cast<Real>(10), static_cast<Real>(0));
    CHECK_ULP_CLOSE(static_cast<Real>(-5), Q, 8);

    // A degenerate interval integrates to exactly zero.
    Q = trapezoidal(c, static_cast<Real>(10), static_cast<Real>(10));
    CHECK_EQUAL(static_cast<Real>(0), Q);

    // Integral of x over [0, 1] is 1/2; the trapezoidal rule is exact for lines.
    auto line = [](Real x)->Real { return x; };
    Q = trapezoidal(line, static_cast<Real>(0), static_cast<Real>(1));
    CHECK_ULP_CLOSE(static_cast<Real>(0.5), Q, 8);
}

// Smooth periodic integrands: the adaptive trapezoidal rule converges rapidly.
template <class Real>
void test_periodic()
{
    using boost::math::quadrature::trapezoidal;

    const Real tol = 100 * std::numeric_limits<Real>::epsilon();

    // Integral of sin(10 x)^2 over [0, pi] is pi/2.
    auto s = [](Real x)->Real { return std::sin(10 * x) * std::sin(10 * x); };
    Real Q = trapezoidal(s, static_cast<Real>(0), boost::math::constants::pi<Real>(), tol);
    CHECK_MOLLIFIED_CLOSE(boost::math::constants::half_pi<Real>(), Q, 10 * tol);

    // Integral of 1/(5 - 4 cos x) over [0, 2 pi] is 2 pi / 3.
    auto r = [](Real x)->Real { return static_cast<Real>(1) / (static_cast<Real>(5) - static_cast<Real>(4) * std::cos(x)); };
    Q = trapezoidal(r, static_cast<Real>(0), boost::math::constants::two_pi<Real>(), tol);
    const Real expected = boost::math::constants::two_pi<Real>() / static_cast<Real>(3);
    CHECK_MOLLIFIED_CLOSE(expected, Q, 10 * tol);
}

int main()
{
    test_polynomial<double>();
    test_periodic<double>();

    test_polynomial<float>();
    test_periodic<float>();

    return boost::math::test::report_errors();
}
