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

    // P(1) = P_n:
    std::array<Real, 2> endpoint{1,2};
    bp.edit_control_point(endpoint, 1);
    CHECK_ULP_CLOSE(endpoint[0], bp(1)[0], 3);
    CHECK_ULP_CLOSE(endpoint[1], bp(1)[1], 3);

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
    auto computed_dp = bp.prime(0);
    CHECK_ULP_CLOSE(2*(control_points[1][0] - control_points[0][0]), computed_dp[0], 3);
    CHECK_ULP_CLOSE(2*(control_points[1][1] - control_points[0][1]), computed_dp[1], 3);

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

// Reversal Symmetry: If q(t) is the Bezier polynomial which consumes the control points in reversed order from p(t),
// then p(t) = q(1-t).
template<typename Real>
void test_reversal_symmetry()
{
    std::vector<std::array<Real, 3>> control_points(2);
    std::uniform_real_distribution<Real> dis(-1,1);
    std::mt19937_64 gen;
    for (size_t i = 0; i < control_points.size(); ++i) {
        for (size_t j = 0; j < 3; ++j) {
            control_points[i][j] = dis(gen);
        }
    }

    auto control_points_copy = control_points;
    auto bp0 = bezier_polynomial(std::move(control_points_copy));

    std::reverse(control_points.begin(), control_points.end());
    auto bp1 = bezier_polynomial(std::move(control_points));

    for (Real t = 0; t <= 1; t += 1.0/16) {
        auto P0 = bp0(t);
        auto P1 = bp1(1-t);
        CHECK_ULP_CLOSE(P0[0], P1[0], 3);
    }
}

// Linear precision: If all control points lie *equidistantly* on a line, then the Bezier curve falls on a line.
// See Bezier and B-spline techniques, Section 2.8, Remark 8.
template<typename Real>
void test_linear_precision()
{
    std::vector<std::array<Real, 3>> control_points(10);
    std::array<Real, 3> P0 = {1,1,1};
    std::array<Real, 3> Pf = {2,2,2};
    control_points[0] = P0;
    control_points[9] = Pf;
    for (size_t i = 1; i < 9; ++i) {
        Real t = Real(i)/(control_points.size()-1);
        control_points[i][0] = (1-t)*P0[0] + t*Pf[0];
        control_points[i][1] = (1-t)*P0[1] + t*Pf[1];
        control_points[i][2] = (1-t)*P0[2] + t*Pf[2];
    }

    auto bp = bezier_polynomial(std::move(control_points));
    for (Real t = 0; t < 1; t += Real(1)/32) {
        std::array<Real, 3> P;
        P[0] = (1-t)*P0[0] + t*Pf[0];
        P[1] = (1-t)*P0[1] + t*Pf[1];
        P[2] = (1-t)*P0[2] + t*Pf[2];

        auto computed = bp(t);
        CHECK_ULP_CLOSE(P[0], computed[0], 3);
        CHECK_ULP_CLOSE(P[1], computed[1], 3);
        CHECK_ULP_CLOSE(P[2], computed[2], 3);

        std::array<Real, 3> dP;
        dP[0] = Pf[0] - P0[0];
        dP[1] = Pf[1] - P0[1];
        dP[2] = Pf[2] - P0[2];
        auto dpComputed = bp.prime(t);
        CHECK_ULP_CLOSE(dP[0], dpComputed[0], 5);
    }
}

template<typename Real>
void test_indefinite_integral()
{
    std::vector<std::array<Real, 3>> control_points(2);
    std::uniform_real_distribution<Real> dis(-1,1);
    std::mt19937_64 gen;
    for (size_t i = 0; i < control_points.size(); ++i) {
        for (size_t j = 0; j < 3; ++j) {
            control_points[i][j] = dis(gen);
        }
    }

    auto control_points_copy = control_points;
    auto bp = bezier_polynomial(std::move(control_points_copy));
    auto bp_int = bp.indefinite_integral();
    for (Real t = 0; t < 1; t += Real(1)/64) {
        auto expected = bp(t);
        auto computed = bp_int.prime(t);
        CHECK_ULP_CLOSE(expected[0], computed[0], 3);
        CHECK_ULP_CLOSE(expected[1], computed[1], 3);
        CHECK_ULP_CLOSE(expected[2], computed[2], 3);
    }

    auto I0 = bp_int(Real(0));
    auto I1 = bp_int(Real(1));
    std::array<Real, 3> avg;
    for (size_t j = 0; j < 3; ++j) {
        avg[j] = (control_points[0][j] + control_points[1][j])/2;
    }
    auto pnts = bp_int.control_points();
    CHECK_EQUAL(pnts.size(), control_points.size() + 1);
    /*std::cout << "I[bp](0) = {" << I0[0] << ", " << I0[1] << ", " << I0[2] << "}\n";
    std::cout << "I[bp](1) = {" << I1[0] << ", " << I1[1] << ", " << I1[2] << "}\n";
    std::cout << "avg      = {" << avg[0] << ", " << avg[1] << ", " << avg[2] << "}\n";
    std::cout << "c[2] =     {" << pnts.back()[0] << ", " << pnts.back()[1] << ", " << pnts.back()[2] << "}\n";*/
}

int main()
{
    test_linear<float>();
    test_linear<double>();
    test_quadratic<float>();
    test_quadratic<double>();
    test_convex_hull<float>();
    test_convex_hull<double>();
    test_linear_precision<float>();
    test_linear_precision<double>();
    test_reversal_symmetry<float>();
    test_reversal_symmetry<double>();
    test_indefinite_integral<float>();
    test_indefinite_integral<double>();
#ifdef BOOST_HAS_FLOAT128
    test_linear<float128>();
    test_quadratic<float128>();
    test_convex_hull<float128>();
#endif
    return boost::math::test::report_errors();
}
