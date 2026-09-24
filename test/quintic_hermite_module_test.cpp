/*
 * Copyright Nick Thompson, 2020
 * Copyright Matt Borland, 2026
 * Use, modification and distribution are subject to the
 * Boost Software License, Version 1.0. (See accompanying file
 * LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
 */
//
// Module build check for quintic_hermite. Consuming the interpolators through
// the boost.math module exercises all three public classes (quintic_hermite,
// cardinal_quintic_hermite, and the array-of-structs cardinal_quintic_hermite_aos)
// and their operator()/prime/double_prime/domain entry points from a module
// consumer. The polynomial-reproduction properties used here are mathematically
// certain: a quintic Hermite interpolator reproduces polynomials up to degree
// five exactly given matching values and derivatives at the nodes. Expected
// values and tolerances are copied from quintic_hermite_test.cpp.

#ifndef BOOST_MATH_BUILD_MODULE
#include <boost/math/interpolators/quintic_hermite.hpp>
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

using boost::math::interpolators::quintic_hermite;
using boost::math::interpolators::cardinal_quintic_hermite;
using boost::math::interpolators::cardinal_quintic_hermite_aos;

// quintic_hermite reproduces a constant exactly with zero derivatives.
template <class Real>
void test_constant()
{
    std::vector<Real> x {0, 1, 2, 3, 9, 22, 81};
    std::vector<Real> y(x.size(), Real(7));
    std::vector<Real> dydx(x.size(), Real(0));
    std::vector<Real> d2ydx2(x.size(), Real(0));

    auto qh = quintic_hermite(std::move(x), std::move(y), std::move(dydx), std::move(d2ydx2));
    for (Real t = 0; t <= 81; t += Real(1))
    {
        CHECK_ULP_CLOSE(Real(7), qh(t), 24);
        CHECK_ULP_CLOSE(Real(0), qh.prime(t), 24);
        CHECK_ULP_CLOSE(Real(0), qh.double_prime(t), 24);
    }
}

// quintic_hermite reproduces a cubic exactly: value, first and second derivative.
template <class Real>
void test_cubic()
{
    std::vector<Real> x {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
    std::vector<Real> y(x.size());
    std::vector<Real> dydx(x.size());
    std::vector<Real> d2ydx2(x.size());
    for (std::size_t i = 0; i < x.size(); ++i)
    {
        y[i] = x[i]*x[i]*x[i];
        dydx[i] = 3*x[i]*x[i];
        d2ydx2[i] = 6*x[i];
    }

    auto qh = quintic_hermite(std::move(x), std::move(y), std::move(dydx), std::move(d2ydx2));
    for (Real t = 0; t <= 9; t += Real(0.25))
    {
        CHECK_ULP_CLOSE(t*t*t, qh(t), 10);
        CHECK_ULP_CLOSE(3*t*t, qh.prime(t), 15);
        CHECK_ULP_CLOSE(6*t, qh.double_prime(t), 20);
    }
}

// cardinal_quintic_hermite reproduces a line exactly; also exercises domain().
template <class Real>
void test_cardinal_linear()
{
    std::vector<Real> y {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
    std::vector<Real> dydx(y.size(), Real(1));
    std::vector<Real> d2ydx2(y.size(), Real(0));
    Real x0 {0};
    Real dx {1};

    auto qh = cardinal_quintic_hermite(std::move(y), std::move(dydx), std::move(d2ydx2), x0, dx);
    for (Real t = 0; t <= 9; t += Real(0.25))
    {
        CHECK_ULP_CLOSE(Real(t), qh(t), 2);
        CHECK_ULP_CLOSE(Real(1), qh.prime(t), 2);
        CHECK_ULP_CLOSE(Real(0), qh.double_prime(t), 2);
    }

    // Evaluate at the reported domain endpoints.
    auto dom = qh.domain();
    Real tlo {dom.first};
    Real thi {dom.second};
    CHECK_ULP_CLOSE(tlo, qh(tlo), 2);
    CHECK_ULP_CLOSE(thi, qh(thi), 2);
    CHECK_ULP_CLOSE(Real(1), qh.prime(tlo), 2);
    CHECK_ULP_CLOSE(Real(1), qh.prime(thi), 128);
}

// cardinal_quintic_hermite_aos reproduces a line exactly from array-of-structs data.
template <class Real>
void test_cardinal_aos_linear()
{
    std::vector<std::array<Real, 3>> data(10);
    for (std::size_t i = 0; i < data.size(); ++i)
    {
        data[i][0] = Real(i);
        data[i][1] = Real(1);
        data[i][2] = Real(0);
    }
    Real x0 {0};
    Real dx {1};

    auto qh_aos = cardinal_quintic_hermite_aos(std::move(data), x0, dx);
    for (Real t = 0; t <= 9; t += Real(0.25))
    {
        CHECK_ULP_CLOSE(Real(t), qh_aos(t), 2);
        CHECK_ULP_CLOSE(Real(1), qh_aos.prime(t), 2);
        CHECK_ULP_CLOSE(Real(0), qh_aos.double_prime(t), 2);
    }
}

int main()
{
    test_constant<float>();
    test_constant<double>();
    test_cubic<float>();
    test_cubic<double>();
    test_cardinal_linear<float>();
    test_cardinal_linear<double>();
    test_cardinal_aos_linear<float>();
    test_cardinal_aos_linear<double>();

    return boost::math::test::report_errors();
}
