/*
 * Copyright Nick Thompson, 2020
 * Copyright Matt Borland, 2026
 * Use, modification and distribution are subject to the
 * Boost Software License, Version 1.0. (See accompanying file
 * LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
 */
//
// Module build check for the septic_hermite interpolators. Consuming the three
// exported classes through the boost.math module exercises operator(), prime(),
// double_prime() and domain(). The reproduction properties used here are
// mathematically certain: a septic Hermite interpolator reproduces polynomials
// up to degree seven exactly. Expected values and ULP tolerances are copied from
// septic_hermite_test.cpp.

#ifndef BOOST_MATH_BUILD_MODULE
#include <boost/math/interpolators/septic_hermite.hpp>
#else
import boost.math;
#endif

#ifndef BOOST_MATH_BUILD_MODULE
#include <cstddef>
#include <vector>
#include <array>
#include <utility>
#endif
#include "math_unit_test.hpp"

using boost::math::interpolators::septic_hermite;
using boost::math::interpolators::cardinal_septic_hermite;
using boost::math::interpolators::cardinal_septic_hermite_aos;

// A septic Hermite interpolator reproduces a constant with zero derivatives.
template <class Real>
void test_constant()
{
    std::vector<Real> x {0, 1, 2, 3, 9, 22, 81};
    std::vector<Real> y(x.size(), Real(7));
    std::vector<Real> dydx(x.size(), Real(0));
    std::vector<Real> d2ydx2(x.size(), Real(0));
    std::vector<Real> d3ydx3(x.size(), Real(0));
    auto sh = septic_hermite(std::move(x), std::move(y), std::move(dydx), std::move(d2ydx2), std::move(d3ydx3));

    for (Real t = 0; t <= 81; t += Real(0.25))
    {
        CHECK_ULP_CLOSE(Real(7), sh(t), 24);
        CHECK_ULP_CLOSE(Real(0), sh.prime(t), 24);
    }

    std::pair<Real, Real> sh_dom = sh.domain();
    CHECK_EQUAL(Real(0), sh_dom.first);
    CHECK_EQUAL(Real(81), sh_dom.second);

    // cardinal_septic_hermite: uniform grid, 128 nodes, x0 = 0, dx = 1.
    Real x0 {0};
    Real dx {1};
    std::vector<Real> yc(128, Real(7));
    std::vector<Real> dydxc(128, Real(0));
    std::vector<Real> d2c(128, Real(0));
    std::vector<Real> d3c(128, Real(0));
    auto csh = cardinal_septic_hermite(std::move(yc), std::move(dydxc), std::move(d2c), std::move(d3c), x0, dx);
    for (Real t = x0; t <= 127; t += Real(0.25))
    {
        CHECK_ULP_CLOSE(Real(7), csh(t), 24);
        CHECK_ULP_CLOSE(Real(0), csh.prime(t), 24);
        CHECK_ULP_CLOSE(Real(0), csh.double_prime(t), 24);
    }

    std::pair<Real, Real> csh_dom = csh.domain();
    CHECK_EQUAL(Real(0), csh_dom.first);
    CHECK_EQUAL(Real(127), csh_dom.second);

    // cardinal_septic_hermite_aos: array-of-structs storage of the same data.
    std::vector<std::array<Real, 4>> data(128);
    for (std::size_t i = 0; i < data.size(); ++i)
    {
        data[i][0] = Real(7);
        data[i][1] = Real(0);
        data[i][2] = Real(0);
        data[i][3] = Real(0);
    }
    auto csh_aos = cardinal_septic_hermite_aos(std::move(data), x0, dx);
    for (Real t = x0; t <= 127; t += Real(0.25))
    {
        CHECK_ULP_CLOSE(Real(7), csh_aos(t), 24);
        CHECK_ULP_CLOSE(Real(0), csh_aos.prime(t), 24);
        CHECK_ULP_CLOSE(Real(0), csh_aos.double_prime(t), 24);
    }
}

// A septic Hermite interpolator reproduces a straight line exactly.
template <class Real>
void test_linear()
{
    std::vector<Real> x {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
    std::vector<Real> y = x;
    std::vector<Real> dydx(x.size(), Real(1));
    std::vector<Real> d2ydx2(x.size(), Real(0));
    std::vector<Real> d3ydx3(x.size(), Real(0));
    auto sh = septic_hermite(std::move(x), std::move(y), std::move(dydx), std::move(d2ydx2), std::move(d3ydx3));

    for (Real t = 0; t <= 9; t += Real(0.25))
    {
        CHECK_ULP_CLOSE(Real(t), sh(t), 2);
        CHECK_ULP_CLOSE(Real(1), sh.prime(t), 2);
    }

    // cardinal_septic_hermite: y[i] = i on a uniform grid.
    Real x0 {0};
    Real dx {1};
    std::vector<Real> yc(10);
    std::vector<Real> dydxc(10, Real(1));
    std::vector<Real> d2c(10, Real(0));
    std::vector<Real> d3c(10, Real(0));
    for (std::size_t i = 0; i < yc.size(); ++i)
    {
        yc[i] = static_cast<Real>(i);
    }
    auto csh = cardinal_septic_hermite(std::move(yc), std::move(dydxc), std::move(d2c), std::move(d3c), x0, dx);
    for (Real t = 0; t <= 9; t += Real(0.125))
    {
        CHECK_ULP_CLOSE(t, csh(t), 15);
        CHECK_ULP_CLOSE(Real(1), csh.prime(t), 15);
        CHECK_ULP_CLOSE(Real(0), csh.double_prime(t), 15);
    }

    // cardinal_septic_hermite_aos: array-of-structs storage of the same line.
    std::vector<std::array<Real, 4>> data(10);
    for (std::size_t i = 0; i < data.size(); ++i)
    {
        data[i][0] = static_cast<Real>(i);
        data[i][1] = Real(1);
        data[i][2] = Real(0);
        data[i][3] = Real(0);
    }
    auto csh_aos = cardinal_septic_hermite_aos(std::move(data), x0, dx);
    for (Real t = 0; t <= 9; t += Real(0.125))
    {
        CHECK_ULP_CLOSE(t, csh_aos(t), 15);
        CHECK_ULP_CLOSE(Real(1), csh_aos.prime(t), 15);
        CHECK_ULP_CLOSE(Real(0), csh_aos.double_prime(t), 15);
    }
}

// A septic Hermite interpolator reproduces a quadratic exactly, exercising a
// nontrivial second derivative (double_prime == 1 for y = x*x/2).
template <class Real>
void test_quadratic()
{
    Real x0 {0};
    Real dx {1};
    std::vector<Real> y(10);
    std::vector<Real> dydx(10);
    std::vector<Real> d2ydx2(10, Real(1));
    std::vector<Real> d3ydx3(10, Real(0));
    for (std::size_t i = 0; i < y.size(); ++i)
    {
        y[i] = static_cast<Real>(i)*static_cast<Real>(i)/Real(2);
        dydx[i] = static_cast<Real>(i);
    }
    auto csh = cardinal_septic_hermite(std::move(y), std::move(dydx), std::move(d2ydx2), std::move(d3ydx3), x0, dx);
    for (Real t = x0; t <= 9; t += Real(0.125))
    {
        CHECK_ULP_CLOSE(t*t/2, csh(t), 24);
        CHECK_ULP_CLOSE(t, csh.prime(t), 24);
        CHECK_ULP_CLOSE(Real(1), csh.double_prime(t), 24);
    }

    std::vector<std::array<Real, 4>> data(10);
    for (std::size_t i = 0; i < data.size(); ++i)
    {
        data[i][0] = static_cast<Real>(i)*static_cast<Real>(i)/Real(2);
        data[i][1] = static_cast<Real>(i);
        data[i][2] = Real(1);
        data[i][3] = Real(0);
    }
    auto csh_aos = cardinal_septic_hermite_aos(std::move(data), x0, dx);
    for (Real t = x0; t <= 9; t += Real(0.125))
    {
        CHECK_ULP_CLOSE(t*t/2, csh_aos(t), 24);
        CHECK_ULP_CLOSE(t, csh_aos.prime(t), 24);
        CHECK_ULP_CLOSE(Real(1), csh_aos.double_prime(t), 24);
    }
}

int main()
{
    test_constant<float>();
    test_constant<double>();
    test_linear<float>();
    test_linear<double>();
    test_quadratic<double>();

    return boost::math::test::report_errors();
}
