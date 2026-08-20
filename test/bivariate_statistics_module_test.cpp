/*
 *  (C) Copyright Nick Thompson 2018.
 *  (C) Copyright Matt Borland 2021.
 *  Use, modification and distribution are subject to the
 *  Boost Software License, Version 1.0. (See accompanying file
 *  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
 *
 *  Module build check for boost::math::statistics::means_and_covariance,
 *  covariance and correlation_coefficient. It exercises the exported entry
 *  points when the test is built against the boost.math module. The double and
 *  float expected values mirror the assertions in bivariate_statistics_test.cpp.
 */

#ifndef BOOST_MATH_BUILD_MODULE
#include <boost/math/statistics/bivariate_statistics.hpp>
#else
import boost.math;
#endif

#ifndef BOOST_MATH_BUILD_MODULE
#include <vector>
#include <tuple>
#include <limits>
#include <cmath>
#endif
#include "math_unit_test.hpp"

// means_and_covariance returns (mean_u, mean_v, covariance); covariance() pulls
// out the third element. All expected values below are exact for these inputs.
template <typename Real>
void test_covariance()
{
    using boost::math::statistics::means_and_covariance;
    using boost::math::statistics::covariance;

    const Real tol = std::numeric_limits<Real>::epsilon();

    // Covariance of a single element is zero; the means are the elements themselves.
    std::vector<Real> u1{static_cast<Real>(8)};
    std::vector<Real> v1{static_cast<Real>(17)};
    auto temp = means_and_covariance(u1, v1);
    CHECK_ABSOLUTE_ERROR(static_cast<Real>(8), std::get<0>(temp), tol);
    CHECK_ABSOLUTE_ERROR(static_cast<Real>(17), std::get<1>(temp), tol);
    CHECK_ABSOLUTE_ERROR(static_cast<Real>(0), std::get<2>(temp), tol);

    // u = {8,4}, v = {3,7}: mean_u = 6, mean_v = 5, covariance = -4.
    std::vector<Real> u2{static_cast<Real>(8), static_cast<Real>(4)};
    std::vector<Real> v2{static_cast<Real>(3), static_cast<Real>(7)};
    temp = means_and_covariance(u2, v2);
    CHECK_ABSOLUTE_ERROR(static_cast<Real>(6), std::get<0>(temp), tol);
    CHECK_ABSOLUTE_ERROR(static_cast<Real>(5), std::get<1>(temp), tol);
    CHECK_ABSOLUTE_ERROR(static_cast<Real>(-4), std::get<2>(temp), tol);

    // v is constant, so cov(u,v) = 0 for any u; mean_u = 2, mean_v = 1.
    std::vector<Real> u3{static_cast<Real>(1), static_cast<Real>(2), static_cast<Real>(3)};
    std::vector<Real> v3{static_cast<Real>(1), static_cast<Real>(1), static_cast<Real>(1)};
    temp = means_and_covariance(u3, v3);
    CHECK_ABSOLUTE_ERROR(static_cast<Real>(2), std::get<0>(temp), tol);
    CHECK_ABSOLUTE_ERROR(static_cast<Real>(1), std::get<1>(temp), tol);
    CHECK_ABSOLUTE_ERROR(static_cast<Real>(0), std::get<2>(temp), tol);

    // covariance() must agree with the third tuple element.
    CHECK_ABSOLUTE_ERROR(static_cast<Real>(0), covariance(u3, v3), tol);
    // Covariance is symmetric: cov(v,u) = cov(u,v).
    CHECK_ABSOLUTE_ERROR(static_cast<Real>(0), covariance(v3, u3), tol);
    // cov(u,u) = variance(u) = 2/3 for {1,2,3}.
    CHECK_ABSOLUTE_ERROR(static_cast<Real>(2) / static_cast<Real>(3), covariance(u3, u3), tol);
}

// correlation_coefficient returns rho in [-1, 1], or NaN when a sample is constant.
template <typename Real>
void test_correlation_coefficient()
{
    using boost::math::statistics::correlation_coefficient;

    const Real tol = std::numeric_limits<Real>::epsilon();

    // A single sample has an undefined correlation coefficient.
    std::vector<Real> a{static_cast<Real>(1)};
    std::vector<Real> b{static_cast<Real>(1)};
    CHECK_NAN(correlation_coefficient(a, b));

    // Two constant samples are likewise undefined.
    std::vector<Real> u{static_cast<Real>(1), static_cast<Real>(1)};
    std::vector<Real> v{static_cast<Real>(1), static_cast<Real>(1)};
    CHECK_NAN(correlation_coefficient(u, v));

    // Perfectly correlated data: rho = 1.
    u = {static_cast<Real>(1), static_cast<Real>(2), static_cast<Real>(3)};
    v = {static_cast<Real>(1), static_cast<Real>(2), static_cast<Real>(3)};
    CHECK_ABSOLUTE_ERROR(static_cast<Real>(1), correlation_coefficient(u, v), tol);

    // Perfectly anti-correlated data: rho = -1, and rho is symmetric in its arguments.
    u = {static_cast<Real>(1), static_cast<Real>(2), static_cast<Real>(3)};
    v = {static_cast<Real>(-1), static_cast<Real>(-2), static_cast<Real>(-3)};
    CHECK_ABSOLUTE_ERROR(static_cast<Real>(-1), correlation_coefficient(u, v), tol);
    CHECK_ABSOLUTE_ERROR(static_cast<Real>(-1), correlation_coefficient(v, u), tol);

    // One constant dataset => correlation is undefined.
    u = {static_cast<Real>(1), static_cast<Real>(2), static_cast<Real>(3)};
    v = {static_cast<Real>(0), static_cast<Real>(0), static_cast<Real>(0)};
    CHECK_NAN(correlation_coefficient(v, u));

    // v = {0,0,3} against u = {1,2,3}: rho = sqrt(3)/2.
    u = {static_cast<Real>(1), static_cast<Real>(2), static_cast<Real>(3)};
    v = {static_cast<Real>(0), static_cast<Real>(0), static_cast<Real>(3)};
    CHECK_ABSOLUTE_ERROR(std::sqrt(static_cast<Real>(3)) / static_cast<Real>(2), correlation_coefficient(v, u), tol);
}

int main()
{
    test_covariance<double>();
    test_correlation_coefficient<double>();

    test_covariance<float>();
    test_correlation_coefficient<float>();

    return boost::math::test::report_errors();
}
