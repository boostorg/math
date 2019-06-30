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
#include <boost/math/special_functions/cardinal_b_spline.hpp>
#ifdef BOOST_HAS_FLOAT128
#include <boost/multiprecision/float128.hpp>
using boost::multiprecision::float128;
#endif

using std::abs;
using boost::math::cardinal_b_spline;
using boost::math::cardinal_b_spline_prime;
using boost::math::forward_cardinal_b_spline;

template<class Real>
void test_box()
{
    Real t = cardinal_b_spline(0,1.1);
    Real expected = 0;
    CHECK_ULP_CLOSE(expected, t, 0);

    t = cardinal_b_spline(0, -1.1);
    expected = 0;
    CHECK_ULP_CLOSE(expected, t, 0);

    Real h = Real(1)/Real(256);
    for (Real t = -Real(1)/Real(2)+h; t < Real(1)/Real(2); t += h)
    {
        expected = 1;
        CHECK_ULP_CLOSE(expected, cardinal_b_spline(0, t), 0);
    }

    for (Real t = h; t < 1; t += h)
    {
        expected = 1;
        CHECK_ULP_CLOSE(expected, forward_cardinal_b_spline(0, t), 0);
    }
}

template<class Real>
void test_hat()
{
    Real t = cardinal_b_spline(1,2.1);
    Real expected = 0;
    CHECK_ULP_CLOSE(expected, t, 0);

    t = cardinal_b_spline(1, -2.1);
    expected = 0;
    CHECK_ULP_CLOSE(expected, t, 0);

    Real h = Real(1)/Real(256);
    for (Real t = -1; t <= 1; t += h)
    {
        expected = 1-abs(t);
        if(!CHECK_ULP_CLOSE(expected, cardinal_b_spline(1, t), 0) )
        {
            std::cerr << "  Problem at t = " << t << "\n";
        }
        if (t < 0 && t > -1) {
            Real expected = 1;
            if(!CHECK_ULP_CLOSE(expected, cardinal_b_spline_prime(1, t), 0) )
            {
                std::cerr << "  Problem at t = " << t << "\n";
            }
        }
        else if (t > 0 && t < 1) {
            if(!CHECK_ULP_CLOSE(-1, cardinal_b_spline_prime(1, t), 0) )
            {
                std::cerr << "  Problem at t = " << t << "\n";
            }
        }
    }

    for (Real t = 0; t < 2; t += h)
    {
        expected = 1 - abs(t-1);
        CHECK_ULP_CLOSE(expected, forward_cardinal_b_spline(1,t), 0);
    }
}

template<class Real>
void test_quadratic()
{
    using std::abs;
    auto b2 = [](Real x) {
        Real absx = abs(x);
        if (absx >= 3/Real(2)) {
            return Real(0);
        }
        if (absx >= 1/Real(2)) {
            Real t = absx - 3/Real(2);
            return t*t/2;
        }
        Real t1 = absx - 1/Real(2);
        Real t2 = absx + 1/Real(2);
        return (2-t1*t1 -t2*t2)/2;
    };

    Real h = 1/Real(256);
    for (Real t = -5; t <= 5; t += h) {
        Real expected = b2(t);
        CHECK_ULP_CLOSE(expected, cardinal_b_spline(2,t), 0);
    }
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
