/*
 * Copyright Nick Thompson, 2019
 * Copyright Matt Borland, 2026
 * Use, modification and distribution are subject to the
 * Boost Software License, Version 1.0. (See accompanying file
 * LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
 *
 * Module build check for boost::math::tools::agm. It exercises the exported
 * arithmetic-geometric mean entry point when the test is built against the
 * boost.math module. The double and float expected values mirror the
 * assertions in agm_test.cpp.
 */

#ifndef BOOST_MATH_BUILD_MODULE
#include <boost/math/tools/agm.hpp>
#else
import boost.math;
#endif

#ifndef BOOST_MATH_BUILD_MODULE
#include <cmath>
#endif
#include "math_unit_test.hpp"

using boost::math::tools::agm;
using std::sqrt;

// Gauss's constant G = 1/agm(sqrt(2), 1). http://oeis.org/A014549/constant
template <typename Real>
void test_gauss_constant()
{
    const double G_expected = 0.83462684167407318628142973279904680899399301349034700244982737010368199270952641186969116035127532412906785;
    Real G_computed = 1 / agm(Real(sqrt(Real(2))), Real(1));
    CHECK_ULP_CLOSE(G_expected, G_computed, 2);
}

// All of these are mathematically-exact properties of the agm.
template <typename Real>
void test_scaling()
{
    Real a = 2;
    Real g = 1;
    Real scale = 7;

    // Homogeneity: agm(k*x, k*y) = k*agm(x, y).
    Real expected = agm(scale * a, scale * g);
    Real computed = scale * agm(a, g);
    CHECK_ULP_CLOSE(expected, computed, 2);

    // agm(a, 0) = 0 and agm(0, 0) = 0.
    CHECK_ULP_CLOSE(Real(0), agm(a, Real(0)), 0);
    CHECK_ULP_CLOSE(Real(0), agm(Real(0), Real(0)), 0);

    // agm(x, x) = x.
    CHECK_ULP_CLOSE(Real(1), agm(Real(1), Real(1)), 0);
    CHECK_ULP_CLOSE(Real(7), agm(Real(7), Real(7)), 0);

    // Composition identity: agm(x, y) = agm((x+y)/2, sqrt(x*y)).
    // With x = 3, y = 1: (x+y)/2 = 2, sqrt(x*y) = sqrt(3).
    expected = agm(Real(3), Real(1));
    computed = agm(Real(2), Real(sqrt(Real(3))));
    CHECK_ULP_CLOSE(expected, computed, 0);
}

int main()
{
    test_gauss_constant<double>();
    test_scaling<double>();

    test_gauss_constant<float>();
    test_scaling<float>();

    return boost::math::test::report_errors();
}
