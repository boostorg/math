/*
 * Copyright Nick Thompson, 2020
 * Use, modification and distribution are subject to the
 * Boost Software License, Version 1.0. (See accompanying file
 * LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
 */

#include "math_unit_test.hpp"
#include <numeric>
#include <utility>
#include <random>
#include <boost/math/interpolators/quintic_hermite.hpp>
#include <boost/circular_buffer.hpp>
#ifdef BOOST_HAS_FLOAT128
#include <boost/multiprecision/float128.hpp>
using boost::multiprecision::float128;
#endif


using boost::math::interpolators::quintic_hermite;

template<typename Real>
void test_constant()
{

    std::vector<Real> x{0,1,2,3, 9, 22, 81};
    std::vector<Real> y(x.size());
    std::vector<Real> dydx(x.size(), 0);
    std::vector<Real> d2ydx2(x.size(), 0);
    for (auto & t : y) {
        t = 7;
    }

    auto qh = quintic_hermite(std::move(x), std::move(y), std::move(dydx), std::move(d2ydx2));

    for (Real t = 0; t <= 81; t += 0.25) {
        CHECK_ULP_CLOSE(Real(7), qh(t), 24);
        CHECK_ULP_CLOSE(Real(0), qh.prime(t), 24);
    }
}


template<typename Real>
void test_linear()
{

    std::vector<Real> x{0,1,2,3, 4,5,6,7,8,9};
    std::vector<Real> y = x;
    std::vector<Real> dydx(x.size(), 1);
    std::vector<Real> d2ydx2(x.size(), 0);

    auto qh = quintic_hermite(std::move(x), std::move(y), std::move(dydx), std::move(d2ydx2));

    for (Real t = 0; t <= 9; t += 0.25) {
        CHECK_ULP_CLOSE(Real(t), qh(t), 2);
        CHECK_ULP_CLOSE(Real(1), qh.prime(t), 2);
    }
}

template<typename Real>
void test_quadratic()
{

    std::vector<Real> x{0,1,2,3, 4,5,6,7,8,9};
    std::vector<Real> y(x.size());
    for (size_t i = 0; i < y.size(); ++i)
    {
        y[i] = x[i]*x[i]/2;
    }

    std::vector<Real> dydx(x.size());
    for (size_t i = 0; i < y.size(); ++i) {
        dydx[i] = x[i];
    }

    std::vector<Real> d2ydx2(x.size(), 1);

    auto qh = quintic_hermite(std::move(x), std::move(y), std::move(dydx), std::move(d2ydx2));

    for (Real t = 0; t <= 9; t += 0.0078125) {
        CHECK_ULP_CLOSE(Real(t*t)/2, qh(t), 2);
        CHECK_ULP_CLOSE(t, qh.prime(t), 2);
    }
}

template<typename Real>
void test_cubic()
{

    std::vector<Real> x{0,1,2,3, 4,5,6,7,8,9};
    std::vector<Real> y(x.size());
    for (size_t i = 0; i < y.size(); ++i)
    {
        y[i] = x[i]*x[i]*x[i]/6;
    }

    std::vector<Real> dydx(x.size());
    for (size_t i = 0; i < y.size(); ++i) {
        dydx[i] = x[i]*x[i]/2;
    }

    std::vector<Real> d2ydx2(x.size());
    for (size_t i = 0; i < y.size(); ++i) {
        d2ydx2[i] = x[i];
    }

    auto qh = quintic_hermite(std::move(x), std::move(y), std::move(dydx), std::move(d2ydx2));

    for (Real t = 0; t <= 9; t += 0.0078125) {
        CHECK_ULP_CLOSE(Real(t*t*t)/6, qh(t), 10);
        //CHECK_ULP_CLOSE(t, qh.prime(t), 2);
    }
}

template<typename Real>
void test_quartic()
{

    std::vector<Real> x{0,1,2,3, 4,5,6,7,8,9, 10, 11};
    std::vector<Real> y(x.size());
    for (size_t i = 0; i < y.size(); ++i)
    {
        y[i] = x[i]*x[i]*x[i]*x[i]/24;
    }

    std::vector<Real> dydx(x.size());
    for (size_t i = 0; i < y.size(); ++i) {
        dydx[i] = x[i]*x[i]*x[i]/6;
    }

    std::vector<Real> d2ydx2(x.size());
    for (size_t i = 0; i < y.size(); ++i) {
        d2ydx2[i] = x[i]*x[i]/2;
    }

    auto qh = quintic_hermite(std::move(x), std::move(y), std::move(dydx), std::move(d2ydx2));

    for (Real t = 1; t <= 11; t += 0.0078125) {
        CHECK_ULP_CLOSE(Real(t*t*t*t)/24, qh(t), 100);
    }
}


template<typename Real>
void test_interpolation_condition()
{
    for (size_t n = 4; n < 50; ++n) {
        std::vector<Real> x(n);
        std::vector<Real> y(n);
        std::vector<Real> dydx(n);
        std::vector<Real> d2ydx2(n);
        std::default_random_engine rd;
        std::uniform_real_distribution<Real> dis(0,1);
        Real x0 = dis(rd);
        x[0] = x0;
        y[0] = dis(rd);
        for (size_t i = 1; i < n; ++i) {
            x[i] = x[i-1] + dis(rd);
            y[i] = dis(rd);
            dydx[i] = dis(rd);
            d2ydx2[i] = dis(rd);
        }

        auto x_copy = x;
        auto y_copy = y;
        auto dydx_copy = dydx;
        auto d2ydx2_copy = d2ydx2;
        auto s = quintic_hermite(std::move(x_copy), std::move(y_copy), std::move(dydx_copy), std::move(d2ydx2_copy));
        //std::cout << "s = " << s << "\n";
        for (size_t i = 0; i < x.size(); ++i) {
            CHECK_ULP_CLOSE(y[i], s(x[i]), 2);
            CHECK_ULP_CLOSE(dydx[i], s.prime(x[i]), 2);
        }
    }
}


int main()
{
    test_constant<float>();
    test_linear<float>();
    test_quadratic<float>();
    test_cubic<float>();
    test_quartic<float>();
    test_interpolation_condition<float>();

    test_constant<double>();
    test_linear<double>();
    test_quadratic<double>();
    test_cubic<double>();
    test_quartic<double>();
    test_interpolation_condition<double>();

    test_constant<long double>();
    test_linear<long double>();
    test_quadratic<long double>();
    test_cubic<long double>();
    test_quartic<long double>();
    test_interpolation_condition<long double>();

#ifdef BOOST_HAS_FLOAT128
    test_constant<float128>();
    test_linear<float128>();
    test_quadratic<float128>();
    test_cubic<float128>();
    test_quartic<float128>();
#endif

    return boost::math::test::report_errors();
}
