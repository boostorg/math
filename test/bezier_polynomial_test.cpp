/*
 * Copyright Nick Thompson, 2021
 * Use, modification and distribution are subject to the
 * Boost Software License, Version 1.0. (See accompanying file
 * LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
 */

#include "math_unit_test.hpp"
#include <numeric>
#include <random>
#include <array>
#include <boost/core/demangle.hpp>
#include <boost/math/interpolators/bezier_polynomial.hpp>
#ifdef BOOST_HAS_FLOAT128
#include <boost/multiprecision/float128.hpp>
using boost::multiprecision::float128;
#endif

using boost::math::interpolators::bezier_polynomial;

template<typename Real>
void test_linear()
{
    std::vector<std::array<Real, 2>> control_points(2);
    control_points[0] = {0.0, 0.0};
    control_points[1] = {1.0, 1.0};
    auto control_points_copy = control_points;
    auto bp = bezier_polynomial(std::move(control_points_copy));

    // P(0) = P_0:
    CHECK_ULP_CLOSE(control_points[0][0], bp(0)[0], 3);
    CHECK_ULP_CLOSE(control_points[0][1], bp(0)[1], 3);

    // P(1) = P_n:
    CHECK_ULP_CLOSE(control_points[1][0], bp(1)[0], 3);
    CHECK_ULP_CLOSE(control_points[1][1], bp(1)[1], 3);

    for (Real t = Real(1)/32; t < 1; t += Real(1)/32) {
        Real expected0 = (1-t)*control_points[0][0] + t*control_points[1][0];
        CHECK_ULP_CLOSE(expected0, bp(t)[0], 3);
    }
}

template<typename Real>
void test_quadratic()
{
    std::vector<std::array<Real, 2>> control_points(3);
    control_points[0] = {0.0, 0.0};
    control_points[1] = {1.0, 1.0};
    control_points[2] = {2.0, 2.0};
    auto control_points_copy = control_points;
    auto bp = bezier_polynomial(std::move(control_points_copy));

    // P(0) = P_0:
    auto computed_point = bp(0);
    CHECK_ULP_CLOSE(control_points[0][0], computed_point[0], 3);
    CHECK_ULP_CLOSE(control_points[0][1], computed_point[1], 3);

    // P(1) = P_n:
    computed_point = bp(1);
    CHECK_ULP_CLOSE(control_points[2][0], computed_point[0], 3);
    CHECK_ULP_CLOSE(control_points[2][1], computed_point[1], 3);
}

// All points on a Bezier polynomial fall into the convex hull of the control polygon.
template<typename Real>
void test_convex_hull()
{
    std::vector<std::array<Real, 2>> control_points(4);
    control_points[0] = {0.0, 0.0};
    control_points[1] = {0.0, 1.0};
    control_points[2] = {1.0, 1.0};
    control_points[3] = {1.0, 0.0};
    auto bp = bezier_polynomial(std::move(control_points));

    for (Real t = 0; t < 1; t += Real(1)/32) {
        auto p = bp(t);
        CHECK_LE(p[0], Real(1));
        CHECK_LE(Real(0), p[0]);
        CHECK_LE(p[1], Real(1));
        CHECK_LE(Real(0), p[1]);
    }
}

int main()
{
    test_linear<float>();
    test_linear<double>();
    test_quadratic<float>();
    test_quadratic<double>();
    test_convex_hull<float>();
    test_convex_hull<double>();
#ifdef BOOST_HAS_FLOAT128
    test_linear<float128>();
    test_quadratic<float128>();
    test_convex_hull<float128>();
#endif
    return boost::math::test::report_errors();
}
