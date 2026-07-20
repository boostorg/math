/*
 * Copyright Nick Thompson, 2019
 * Copyright Matt Borland, 2026
 * Use, modification and distribution are subject to the
 * Boost Software License, Version 1.0. (See accompanying file
 * LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
 */
//
// Module build check for cardinal_quintic_b_spline. Consuming the interpolator
// through the boost.math module exercises both constructor overloads (the raw
// pointer form and the std::vector form), the estimate-derivatives path (NaN
// endpoint derivatives), and operator()/prime/double_prime/t_max from a module
// consumer. The reproduction properties used here are mathematically certain:
// a quintic B-spline reproduces polynomials up to degree five exactly. Expected
// values and tolerances are copied from cardinal_quintic_b_spline_test.cpp.

#ifndef BOOST_MATH_BUILD_MODULE
#include <boost/math/interpolators/cardinal_quintic_b_spline.hpp>
#else
import boost.math;
#endif

#ifndef BOOST_MATH_BUILD_MODULE
#include <cstddef>
#include <vector>
#include <utility>
#include <limits>
#include <cmath>
#endif
#include "math_unit_test.hpp"

using boost::math::interpolators::cardinal_quintic_b_spline;

// A quintic B-spline reproduces a straight line exactly. The raw pointer
// constructor is used with the endpoint derivatives supplied explicitly.
template <class Real>
void test_linear()
{
    using std::abs;
    Real m {Real(8.3)};
    Real b {Real(7.2)};
    Real t0 {0};
    Real h {Real(1)/Real(16)};
    std::size_t n {64};
    std::vector<Real> y(n);
    for (std::size_t i {0}; i < n; ++i)
    {
        Real t {i*h};
        y[i] = m*t + b;
    }
    std::pair<Real, Real> left_endpoint_derivatives {m, Real(0)};
    std::pair<Real, Real> right_endpoint_derivatives {m, Real(0)};
    auto qbs = cardinal_quintic_b_spline<Real>(y.data(), y.size(), t0, h, left_endpoint_derivatives, right_endpoint_derivatives);

    for (std::size_t i {0}; i < n; ++i)
    {
        Real t {t0 + i*h};
        CHECK_ULP_CLOSE(m*t + b, qbs(t), 3);
        CHECK_MOLLIFIED_CLOSE(m, qbs.prime(t), 100*abs(m*t + b)*std::numeric_limits<Real>::epsilon());
        CHECK_MOLLIFIED_CLOSE(Real(0), qbs.double_prime(t), 10000*abs(m*t + b)*std::numeric_limits<Real>::epsilon());
    }

    for (std::size_t i {0}; i < n - 1; ++i)
    {
        Real t {t0 + i*h + h/2};
        CHECK_ULP_CLOSE(m*t + b, qbs(t), 5);
        CHECK_MOLLIFIED_CLOSE(m, qbs.prime(t), 1500*std::numeric_limits<Real>::epsilon());
        t = t0 + i*h + h/4;
        CHECK_ULP_CLOSE(m*t + b, qbs(t), 5);
        CHECK_MOLLIFIED_CLOSE(m, qbs.prime(t), 3000*std::numeric_limits<Real>::epsilon());
    }
}

// A quintic B-spline reproduces a constant. The endpoint derivatives are left
// as their NaN default so the constructor estimates them by finite differences.
template <class Real>
void test_constant_estimate_derivatives()
{
    Real c {Real(7.5)};
    Real t0 {0};
    Real h {Real(1)/Real(16)};
    std::size_t n {65};
    std::vector<Real> v(n, c);
    auto qbs = cardinal_quintic_b_spline<Real>(v.data(), v.size(), t0, h);

    for (std::size_t i {0}; i < n; ++i)
    {
        Real t {t0 + i*h};
        CHECK_ULP_CLOSE(c, qbs(t), 3);
        CHECK_MOLLIFIED_CLOSE(Real(0), qbs.prime(t), 1200*std::numeric_limits<Real>::epsilon());
        CHECK_MOLLIFIED_CLOSE(Real(0), qbs.double_prime(t), 200000*std::numeric_limits<Real>::epsilon());
    }

    for (std::size_t i {0}; i < n - 1; ++i)
    {
        Real t {t0 + i*h + h/2};
        CHECK_ULP_CLOSE(c, qbs(t), 8);
        CHECK_MOLLIFIED_CLOSE(Real(0), qbs.prime(t), 1200*std::numeric_limits<Real>::epsilon());
        CHECK_MOLLIFIED_CLOSE(Real(0), qbs.double_prime(t), 80000*std::numeric_limits<Real>::epsilon());
        t = t0 + i*h + h/4;
        CHECK_ULP_CLOSE(c, qbs(t), 5);
        CHECK_MOLLIFIED_CLOSE(Real(0), qbs.prime(t), 1200*std::numeric_limits<Real>::epsilon());
        CHECK_MOLLIFIED_CLOSE(Real(0), qbs.double_prime(t), 38000*std::numeric_limits<Real>::epsilon());
    }
}

// A quintic B-spline reproduces a quadratic exactly. The std::vector
// constructor is used, and t_max is verified against the known domain endpoint.
template <class Real>
void test_quadratic()
{
    Real a {Real(1)/Real(16)};
    Real b {Real(-3.5)};
    Real c {Real(-9)};
    Real t0 {0};
    Real h {Real(1)/Real(16)};
    std::size_t n {65};
    std::vector<Real> y(n);
    for (std::size_t i {0}; i < n; ++i)
    {
        Real t {i*h};
        y[i] = a*t*t + b*t + c;
    }
    Real t_max {t0 + (n - 1)*h};
    std::pair<Real, Real> left_endpoint_derivatives {b, 2*a};
    std::pair<Real, Real> right_endpoint_derivatives {2*a*t_max + b, 2*a};
    auto qbs = cardinal_quintic_b_spline<Real>(y, t0, h, left_endpoint_derivatives, right_endpoint_derivatives);

    CHECK_ULP_CLOSE(t_max, qbs.t_max(), 2);

    for (std::size_t i {0}; i < n; ++i)
    {
        Real t {t0 + i*h};
        CHECK_ULP_CLOSE(a*t*t + b*t + c, qbs(t), 3);
    }

    for (std::size_t i {0}; i < n - 1; ++i)
    {
        Real t {t0 + i*h + h/2};
        CHECK_ULP_CLOSE(a*t*t + b*t + c, qbs(t), 5);
        t = t0 + i*h + h/4;
        CHECK_ULP_CLOSE(a*t*t + b*t + c, qbs(t), 5);
    }
}

int main()
{
    test_linear<float>();
    test_linear<double>();
    test_constant_estimate_derivatives<double>();
    test_quadratic<double>();

    return boost::math::test::report_errors();
}
