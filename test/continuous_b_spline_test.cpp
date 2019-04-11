/*
 * Copyright Nick Thompson, 2019
 * Use, modification and distribution are subject to the
 * Boost Software License, Version 1.0. (See accompanying file
 * LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
 */

#include "math_unit_test.hpp"
#include <numeric>
#include <utility>
#include <random>
#include <boost/core/demangle.hpp>
#include <boost/math/special_functions/b_splines.hpp>
#ifdef BOOST_HAS_FLOAT128
#include <boost/multiprecision/float128.hpp>
using boost::multiprecision::float128;
#endif

using boost::math::b_spline;

template<class Real>
void test_box()
{
    auto box = b_spline<Real>(1);


    Real t = box(1.1);
    Real expected = 0;
    CHECK_ULP_CLOSE(expected, t, 0);

    t = box(-0.1);
    expected = 0;
    CHECK_ULP_CLOSE(expected, t, 0);

    Real h = Real(1)/Real(256);
    for (Real t = h; t < 1; t += h)
    {
        expected = 1;
        CHECK_ULP_CLOSE(expected, box(t), 0);
    }

    for (Real t = -Real(1)/Real(2) + h; t < Real(1)/Real(2); t += h)
    {
        expected = 1;
        CHECK_ULP_CLOSE(expected, box.centered(t), 0);
    }

}

template<class Real>
void test_hat()
{
    auto hat = b_spline<Real>(2);


    Real t = hat(2.1);
    Real expected = 0;
    CHECK_ULP_CLOSE(expected, t, 0);

    t = hat(-0.1);
    expected = 0;
    CHECK_ULP_CLOSE(expected, t, 0);

    Real h = Real(1)/Real(256);
    for (Real t = 0; t <= 1; t += h)
    {
        expected = t;
        CHECK_ULP_CLOSE(expected, hat(t), 0);
    }

    for (Real t = 1 + h; t < 2; t += h)
    {
        expected = 2 - t;
        CHECK_ULP_CLOSE(expected, hat(t), 0);
    }

    for (Real t = -Real(1); t < 0; t += h)
    {
        expected = t + 1;
        CHECK_ULP_CLOSE(expected, hat.centered(t), 0);
    }

    for (Real t = Real(0); t < Real(1); t += h)
    {
        expected = 1 - t;
        CHECK_ULP_CLOSE(expected, hat.centered(t), 0);
    }

}

template<class Real>
void test_quadratic()
{
    auto quad = b_spline<Real>(3);


    Real t = quad(3.1);
    Real expected = 0;
    CHECK_ULP_CLOSE(expected, t, 0);

    t = quad(-0.1);
    expected = 0;
    CHECK_ULP_CLOSE(expected, t, 0);

    t = quad(Real(3)/Real(2));
    expected = Real(3)/Real(4);
    CHECK_ULP_CLOSE(expected, t, 0);

    t = quad(Real(5)/Real(2));
    expected = Real(1)/Real(8);
    CHECK_ULP_CLOSE(expected, t, 0);
}


int main()
{
    test_box<float>();
    test_box<double>();
    test_box<long double>();

    test_hat<float>();
    test_hat<double>();
    test_hat<long double>();

    test_quadratic<float>();
    test_quadratic<double>();
    test_quadratic<long double>();

#ifdef BOOST_HAS_FLOAT128
    test_box<float128>();
    test_hat<float128>();
    test_quadratic<float128>();
#endif

    return boost::math::test::report_errors();
}
