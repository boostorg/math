// Copyright Nick Thompson, 2021
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// Module build check for bezier_polynomial. Consuming the interpolator through
// `import boost.math;` confirms the exported class template and its
// operator()/prime/edit_control_point members resolve from a module consumer
// (along with the std::shared_ptr / std::array members they rely on). The
// assertions exercise only mathematically-certain Bezier properties: endpoint
// interpolation, the exact straight line produced by collinear equidistant
// control points, the degree-n endpoint derivative n*(P_1 - P_0), and the
// convex-hull bound. Expected values are taken from bezier_polynomial_test.cpp.

#ifndef BOOST_MATH_BUILD_MODULE
#include <boost/math/interpolators/bezier_polynomial.hpp>
#else
import boost.math;
#endif

#include "math_unit_test.hpp"

#ifndef BOOST_MATH_BUILD_MODULE
#include <vector>
#include <array>
#include <utility>
#endif

using boost::math::interpolators::bezier_polynomial;

// A Bezier polynomial with two control points is the straight segment
// P(t) = (1-t) P_0 + t P_1, and editing a control point moves the endpoint.
template <class Real>
void test_linear(const char* type_name)
{
    std::cout << "Testing bezier_polynomial linear segment on type " << type_name << std::endl;

    std::vector<std::array<Real, 2>> control_points(2);
    control_points[0] = {Real(0), Real(0)};
    control_points[1] = {Real(1), Real(1)};
    auto control_points_copy = control_points;
    auto bp = bezier_polynomial(std::move(control_points_copy));

    // P(0) == P_0.
    CHECK_ULP_CLOSE(control_points[0][0], bp(0)[0], 3);
    CHECK_ULP_CLOSE(control_points[0][1], bp(0)[1], 3);

    // P(1) == P_1.
    CHECK_ULP_CLOSE(control_points[1][0], bp(1)[0], 3);
    CHECK_ULP_CLOSE(control_points[1][1], bp(1)[1], 3);

    // Interior points fall on the line P(t) = (t, t). 1/32 is dyadic, so the
    // sampled t values are exact.
    for (Real t = Real(1)/32; t < 1; t += Real(1)/32)
    {
        Real expected0 = (1 - t) * control_points[0][0] + t * control_points[1][0];
        CHECK_ULP_CLOSE(expected0, bp(t)[0], 3);
    }

    // Editing the endpoint control point moves P(1) to the new point.
    std::array<Real, 2> endpoint {Real(1), Real(2)};
    bp.edit_control_point(endpoint, 1);
    CHECK_ULP_CLOSE(endpoint[0], bp(1)[0], 3);
    CHECK_ULP_CLOSE(endpoint[1], bp(1)[1], 3);
}

// Three collinear, equidistant control points produce the straight segment
// from P_0 to P_2, so P(0) == P_0, P(1) == P_2, and P'(0) == 2 (P_1 - P_0).
template <class Real>
void test_quadratic(const char* type_name)
{
    std::cout << "Testing bezier_polynomial quadratic on type " << type_name << std::endl;

    std::vector<std::array<Real, 2>> control_points(3);
    control_points[0] = {Real(0), Real(0)};
    control_points[1] = {Real(1), Real(1)};
    control_points[2] = {Real(2), Real(2)};
    auto control_points_copy = control_points;
    auto bp = bezier_polynomial(std::move(control_points_copy));

    // P(0) == P_0.
    auto computed_point = bp(0);
    CHECK_ULP_CLOSE(control_points[0][0], computed_point[0], 3);
    CHECK_ULP_CLOSE(control_points[0][1], computed_point[1], 3);

    // P'(0) == n (P_1 - P_0) with degree n == 2.
    auto computed_dp = bp.prime(0);
    CHECK_ULP_CLOSE(2 * (control_points[1][0] - control_points[0][0]), computed_dp[0], 3);
    CHECK_ULP_CLOSE(2 * (control_points[1][1] - control_points[0][1]), computed_dp[1], 3);

    // P(1) == P_2.
    computed_point = bp(1);
    CHECK_ULP_CLOSE(control_points[2][0], computed_point[0], 3);
    CHECK_ULP_CLOSE(control_points[2][1], computed_point[1], 3);
}

// Every point of a Bezier polynomial lies in the convex hull of its control
// polygon; here that hull is the unit square, and the endpoints are exact.
template <class Real>
void test_convex_hull(const char* type_name)
{
    std::cout << "Testing bezier_polynomial convex hull on type " << type_name << std::endl;

    std::vector<std::array<Real, 2>> control_points(4);
    control_points[0] = {Real(0), Real(0)};
    control_points[1] = {Real(0), Real(1)};
    control_points[2] = {Real(1), Real(1)};
    control_points[3] = {Real(1), Real(0)};
    auto bp = bezier_polynomial(std::move(control_points));

    // Endpoint interpolation: P(0) == (0,0), P(1) == (1,0).
    auto p0 = bp(0);
    CHECK_ULP_CLOSE(Real(0), p0[0], 3);
    CHECK_ULP_CLOSE(Real(0), p0[1], 3);
    auto p1 = bp(1);
    CHECK_ULP_CLOSE(Real(1), p1[0], 3);
    CHECK_ULP_CLOSE(Real(0), p1[1], 3);

    for (Real t = 0; t <= 1; t += Real(1)/32)
    {
        auto p = bp(t);
        CHECK_LE(p[0], Real(1));
        CHECK_LE(Real(0), p[0]);
        CHECK_LE(p[1], Real(1));
        CHECK_LE(Real(0), p[1]);
    }
}

int main()
{
    test_linear<float>("float");
    test_linear<double>("double");

    test_quadratic<float>("float");
    test_quadratic<double>("double");

    test_convex_hull<float>("float");
    test_convex_hull<double>("double");

    return boost::math::test::report_errors();
}
