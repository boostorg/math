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
#include <boost/math/interpolators/cubic_hermite.hpp>
#include <boost/circular_buffer.hpp>
#ifdef BOOST_HAS_FLOAT128
#include <boost/multiprecision/float128.hpp>
using boost::multiprecision::float128;
#endif


using boost::math::interpolators::cubic_hermite;

template<typename Real>
void test_constant()
{

    std::vector<Real> x{0,1,2,3, 9, 22, 81};
    std::vector<Real> y(x.size());
    for (auto & t : y) {
        t = 7;
    }

    std::vector<Real> dydx(x.size(), Real(0));
    auto x_copy = x;
    auto y_copy = y;
    auto dydx_copy = dydx;
    auto hermite_spline = cubic_hermite(std::move(x_copy), std::move(y_copy), std::move(dydx_copy));

    for (Real t = x[0]; t <= x.back(); t += 0.25) {
        CHECK_ULP_CLOSE(Real(7), hermite_spline(t), 2);
        CHECK_ULP_CLOSE(Real(0), hermite_spline.prime(t), 2);
    }

    boost::circular_buffer<Real> x_buf(x.size());
    for (auto & t : x) {
        x_buf.push_back(t);
    }

    boost::circular_buffer<Real> y_buf(x.size());
    for (auto & t : y) {
        y_buf.push_back(t);
    }

    boost::circular_buffer<Real> dydx_buf(x.size());
    for (auto & t : dydx) {
        dydx_buf.push_back(t);
    }

    auto circular_hermite_spline = cubic_hermite(std::move(x_buf), std::move(y_buf), std::move(dydx_buf));

    for (Real t = x[0]; t <= x.back(); t += 0.25) {
        CHECK_ULP_CLOSE(Real(7), circular_hermite_spline(t), 2);
        CHECK_ULP_CLOSE(Real(0), circular_hermite_spline.prime(t), 2);
    }

    circular_hermite_spline.push_back(x.back() + 1, 7, 0);
    CHECK_ULP_CLOSE(Real(0), circular_hermite_spline.prime(x.back()+1), 2);

}

template<typename Real>
void test_linear()
{
    std::vector<Real> x{0,1,2,3};
    std::vector<Real> y{0,1,2,3};
    std::vector<Real> dydx{1,1,1,1};

    auto x_copy = x;
    auto y_copy = y;
    auto dydx_copy = dydx;
    auto hermite_spline = cubic_hermite(std::move(x_copy), std::move(y_copy), std::move(dydx_copy));

    CHECK_ULP_CLOSE(y[0], hermite_spline(x[0]), 0);
    CHECK_ULP_CLOSE(Real(1)/Real(2), hermite_spline(Real(1)/Real(2)), 10);
    CHECK_ULP_CLOSE(y[1], hermite_spline(x[1]), 0);
    CHECK_ULP_CLOSE(Real(3)/Real(2), hermite_spline(Real(3)/Real(2)), 10);
    CHECK_ULP_CLOSE(y[2], hermite_spline(x[2]), 0);
    CHECK_ULP_CLOSE(Real(5)/Real(2), hermite_spline(Real(5)/Real(2)), 10);
    CHECK_ULP_CLOSE(y[3], hermite_spline(x[3]), 0);

    x.resize(45);
    y.resize(45);
    dydx.resize(45);
    for (size_t i = 0; i < x.size(); ++i) {
        x[i] = i;
        y[i] = i;
        dydx[i] = 1;
    }

    x_copy = x;
    y_copy = y;
    dydx_copy = dydx;
    hermite_spline = cubic_hermite(std::move(x_copy), std::move(y_copy), std::move(dydx_copy));
    for (Real t = 0; t < x.back(); t += 0.5) {
        CHECK_ULP_CLOSE(t, hermite_spline(t), 0);
        CHECK_ULP_CLOSE(Real(1), hermite_spline.prime(t), 0);
    }

    boost::circular_buffer<Real> x_buf(x.size());
    for (auto & t : x) {
        x_buf.push_back(t);
    }

    boost::circular_buffer<Real> y_buf(x.size());
    for (auto & t : y) {
        y_buf.push_back(t);
    }

    boost::circular_buffer<Real> dydx_buf(x.size());
    for (auto & t : dydx) {
        dydx_buf.push_back(t);
    }

    auto circular_hermite_spline = cubic_hermite(std::move(x_buf), std::move(y_buf), std::move(dydx_buf));

    for (Real t = x[0]; t <= x.back(); t += 0.25) {
        CHECK_ULP_CLOSE(t, circular_hermite_spline(t), 2);
        CHECK_ULP_CLOSE(Real(1), circular_hermite_spline.prime(t), 2);
    }

    circular_hermite_spline.push_back(x.back() + 1, y.back()+1, 1);

    CHECK_ULP_CLOSE(Real(y.back() + 1), circular_hermite_spline(Real(x.back()+1)), 2);
    CHECK_ULP_CLOSE(Real(1), circular_hermite_spline.prime(Real(x.back()+1)), 2);

}

template<typename Real>
void test_interpolation_condition()
{
    for (size_t n = 4; n < 50; ++n) {
        std::vector<Real> x(n);
        std::vector<Real> y(n);
        std::vector<Real> dydx(n);
        std::default_random_engine rd;
        std::uniform_real_distribution<Real> dis(0,1);
        Real x0 = dis(rd);
        x[0] = x0;
        y[0] = dis(rd);
        for (size_t i = 1; i < n; ++i) {
            x[i] = x[i-1] + dis(rd);
            y[i] = dis(rd);
            dydx[i] = dis(rd);
        }

        auto x_copy = x;
        auto y_copy = y;
        auto dydx_copy = dydx;
        auto s = cubic_hermite(std::move(x_copy), std::move(y_copy), std::move(dydx_copy));
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
    test_interpolation_condition<float>();


    test_constant<double>();
    test_linear<double>();
    test_interpolation_condition<double>();

    test_constant<long double>();
    test_linear<long double>();
    test_interpolation_condition<long double>();

#ifdef BOOST_HAS_FLOAT128
    test_constant<float128>();
    test_linear<float128>();
#endif

    return boost::math::test::report_errors();
}
