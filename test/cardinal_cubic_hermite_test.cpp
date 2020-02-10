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
#include <vector>
#include <array>
#include <boost/math/interpolators/cardinal_cubic_hermite.hpp>
#ifdef BOOST_HAS_FLOAT128
#include <boost/multiprecision/float128.hpp>
using boost::multiprecision::float128;
#endif


using boost::math::interpolators::cardinal_cubic_hermite;
using boost::math::interpolators::cardinal_cubic_hermite_aos;

template<typename Real>
void test_constant()
{
    Real x0 = 0;
    Real dx = 2;
    std::vector<Real> y(25);
    for (auto & t : y) {
        t = 7;
    }

    std::vector<Real> dydx(y.size(), Real(0));

    auto hermite_spline = cardinal_cubic_hermite(std::move(y), std::move(dydx), x0, dx);

    for (Real t = x0; t <= x0 + 24*dx; t += 0.25) {
        CHECK_ULP_CLOSE(Real(7), hermite_spline(t), 2);
        CHECK_ULP_CLOSE(Real(0), hermite_spline.prime(t), 2);
    }

    // Array of structs:

    std::vector<std::array<Real, 2>> data(25);
    for (auto & t : data) {
        t[0] = 7;
        t[1] = 0;
    }
    auto hermite_spline_aos = cardinal_cubic_hermite_aos(std::move(data), x0, dx);

    for (Real t = x0; t <= x0 + 24*dx; t += 0.25) {
        CHECK_ULP_CLOSE(Real(7), hermite_spline_aos(t), 2);
        CHECK_ULP_CLOSE(Real(0), hermite_spline_aos.prime(t), 2);
    }
 
}


template<typename Real>
void test_linear()
{
    Real x0 = 0;
    Real dx = 1;
    std::vector<Real> y{0,1,2,3};
    std::vector<Real> dydx{1,1,1,1};
    auto y_copy = y;
    auto dydx_copy = dydx;
    auto hermite_spline = cardinal_cubic_hermite(std::move(y_copy), std::move(dydx_copy), x0, dx);

    CHECK_ULP_CLOSE(y[0], hermite_spline(0), 0);
    CHECK_ULP_CLOSE(Real(1)/Real(2), hermite_spline(Real(1)/Real(2)), 10);
    CHECK_ULP_CLOSE(y[1], hermite_spline(1), 0);
    CHECK_ULP_CLOSE(Real(3)/Real(2), hermite_spline(Real(3)/Real(2)), 10);
    CHECK_ULP_CLOSE(y[2], hermite_spline(2), 0);
    CHECK_ULP_CLOSE(Real(5)/Real(2), hermite_spline(Real(5)/Real(2)), 10);
    CHECK_ULP_CLOSE(y[3], hermite_spline(3), 0);


    y.resize(45);
    dydx.resize(45);
    for (size_t i = 0; i < y.size(); ++i) {
        y[i] = i;
        dydx[i] = 1;
    }

    hermite_spline = cardinal_cubic_hermite(std::move(y), std::move(dydx), x0, dx);
    for (Real t = 0; t < 44; t += 0.5) {
        CHECK_ULP_CLOSE(t, hermite_spline(t), 0);
        CHECK_ULP_CLOSE(Real(1), hermite_spline.prime(t), 0);
    }

    std::vector<std::array<Real, 2>> data(45);
    for (size_t i = 0; i < data.size(); ++i) {
        data[i][0] = i;
        data[i][1] = 1;
    }

    auto hermite_spline_aos = cardinal_cubic_hermite_aos(std::move(data), x0, dx);
    for (Real t = 0; t < 44; t += 0.5) {
        CHECK_ULP_CLOSE(t, hermite_spline_aos(t), 0);
        CHECK_ULP_CLOSE(Real(1), hermite_spline_aos.prime(t), 0);
    }


}


template<typename Real>
void test_quadratic()
{
    Real x0 = -1;
    Real dx = Real(1)/Real(256);

    std::vector<Real> y(50);
    std::vector<Real> dydx(y.size());
    for (size_t i = 0; i < y.size(); ++i) {
        Real x = x0 + i*dx;
        y[i] = x*x/2;
        dydx[i] = x;
    }

    auto s = cardinal_cubic_hermite(std::move(y), std::move(dydx), x0, dx);
    for (Real t = x0; t <= x0 + 49*dx; t+= 0.0125)
    {
        CHECK_ULP_CLOSE(t*t/2, s(t), 12);
        CHECK_ULP_CLOSE(t, s.prime(t), 70);
    }

    std::vector<std::array<Real, 2>> data(50);
    for (size_t i = 0; i < data.size(); ++i) {
        Real x = x0 + i*dx;
        data[i][0] = x*x/2;
        data[i][1] = x;
    }


    auto saos = cardinal_cubic_hermite_aos(std::move(data), x0, dx);
    for (Real t = x0; t <= x0 + 49*dx; t+= 0.0125)
    {
        CHECK_ULP_CLOSE(t*t/2, saos(t), 12);
        CHECK_ULP_CLOSE(t, saos.prime(t), 70);
    }

}


template<typename Real>
void test_interpolation_condition()
{
    for (size_t n = 4; n < 50; ++n) {
        std::vector<Real> y(n);
        std::vector<Real> dydx(n);
        std::default_random_engine rd;
        std::uniform_real_distribution<Real> dis(0.1,1);
        Real x0 = Real(2);
        Real dx = Real(1)/Real(128);
        for (size_t i = 0; i < n; ++i) {
            y[i] = dis(rd);
            dydx[i] = dis(rd);
        }

        auto y_copy = y;
        auto dydx_copy = dydx;
        auto s = cardinal_cubic_hermite(std::move(y_copy), std::move(dydx_copy), x0, dx);
        for (size_t i = 0; i < y.size(); ++i) {
            CHECK_ULP_CLOSE(y[i], s(x0 + i*dx), 2);
            CHECK_ULP_CLOSE(dydx[i], s.prime(x0 + i*dx), 2);
        }
    }
}

int main()
{
    test_constant<float>();
    test_constant<double>();
    test_constant<long double>();

    test_linear<float>();
    test_linear<double>();
    test_linear<long double>();

    test_quadratic<float>();
    test_quadratic<double>();
    test_quadratic<long double>();

    test_interpolation_condition<float>();
    test_interpolation_condition<double>();
    test_interpolation_condition<long double>();
#ifdef BOOST_HAS_FLOAT128
    test_constant<float128>();
    test_linear<float128>();
    test_quadratic<float128>();
#endif

    return boost::math::test::report_errors();
}
