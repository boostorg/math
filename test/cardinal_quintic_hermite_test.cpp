/*
 * Copyright Nick Thompson, 2020
 * Use, modification and distribution are subject to the
 * Boost Software License, Version 1.0. (See accompanying file
 * LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
 */

#include "math_unit_test.hpp"
#include <numeric>
#include <utility>
#include <boost/random/uniform_real.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/math/interpolators/cardinal_quintic_hermite.hpp>
#ifdef BOOST_HAS_FLOAT128
#include <boost/multiprecision/float128.hpp>
using boost::multiprecision::float128;
#endif


using boost::math::interpolators::cardinal_quintic_hermite;
using boost::math::interpolators::cardinal_quintic_hermite_aos;

template<typename Real>
void test_constant()
{

    std::vector<Real> y(25);
    std::vector<Real> dydx(y.size(), 0);
    std::vector<Real> d2ydx2(y.size(), 0);
    for (auto & t : y) {
        t = 7;
    }
    Real x0 = 4;
    Real dx = Real(1)/Real(8);

    auto qh = cardinal_quintic_hermite(std::move(y), std::move(dydx), std::move(d2ydx2), x0, dx);

    for (Real t = x0; t <= x0 + 24*dx; t += 0.25) {
        CHECK_ULP_CLOSE(Real(7), qh(t), 24);
        //CHECK_ULP_CLOSE(Real(0), qh.prime(t), 24);
    }


}


template<typename Real>
void test_linear()
{
    std::vector<Real> y{0,1,2,3, 4,5,6,7,8,9};
    Real x0 = 0;
    Real dx = 1;
    std::vector<Real> dydx(y.size(), 1);
    std::vector<Real> d2ydx2(y.size(), 0);

    auto qh = cardinal_quintic_hermite(std::move(y), std::move(dydx), std::move(d2ydx2), x0, dx);

    for (Real t = 0; t <= 9; t += 0.25) {
        CHECK_ULP_CLOSE(Real(t), qh(t), 2);
        //CHECK_ULP_CLOSE(Real(1), qh.prime(t), 2);
    }
}

template<typename Real>
void test_quadratic()
{
    Real x0 = 0;
    Real dx = 1;
    std::vector<Real> y(10);
    for (size_t i = 0; i < y.size(); ++i)
    {
        y[i] = i*i/2;
    }

    std::vector<Real> dydx(y.size());
    for (size_t i = 0; i < y.size(); ++i) {
        dydx[i] = i;
    }

    std::vector<Real> d2ydx2(y.size(), 1);

    auto qh = cardinal_quintic_hermite(std::move(y), std::move(dydx), std::move(d2ydx2), x0, dx);

    for (Real t = 0; t <= 9; t += 0.0078125) {
        CHECK_ULP_CLOSE(Real(t*t)/2, qh(t), 2);
        //CHECK_ULP_CLOSE(t, qh.prime(t), 2);
    }

}

template<typename Real>
void test_cubic()
{
    Real x0 = 0;
    Real dx = 1;
    std::vector<Real> y(10);
    for (size_t i = 0; i < y.size(); ++i)
    {
        y[i] = i*i*i/6;
    }

    std::vector<Real> dydx(y.size());
    for (size_t i = 0; i < y.size(); ++i) {
        dydx[i] = i*i/2;
    }

    std::vector<Real> d2ydx2(y.size());
    for (size_t i = 0; i < y.size(); ++i) {
        d2ydx2[i] = i;
    }

    auto qh = cardinal_quintic_hermite(std::move(y), std::move(dydx), std::move(d2ydx2), x0, dx);

    for (Real t = 0; t <= 9; t += 0.0078125) {
        CHECK_ULP_CLOSE(Real(t*t*t)/6, qh(t), 10);
    }
}

template<typename Real>
void test_quartic()
{

    Real x0 = 0;
    Real dx = 1;
    std::vector<Real> y(11);
    for (size_t i = 0; i < y.size(); ++i)
    {
        y[i] = i*i*i*i/24;
    }

    std::vector<Real> dydx(y.size());
    for (size_t i = 0; i < y.size(); ++i) {
        dydx[i] = i*i*i/6;
    }

    std::vector<Real> d2ydx2(y.size());
    for (size_t i = 0; i < y.size(); ++i) {
        d2ydx2[i] = i*i/2;
    }

    auto qh = cardinal_quintic_hermite(std::move(y), std::move(dydx), std::move(d2ydx2), x0, dx);

    for (Real t = 1; t <= 11; t += 0.0078125) {
        CHECK_ULP_CLOSE(Real(t*t*t*t)/24, qh(t), 100);
    }
}


int main()
{
    /*test_constant<float>();
    test_linear<float>();
    test_quadratic<float>();
    test_cubic<float>();
    test_quartic<float>();*/

    test_constant<double>();
    test_linear<double>();
    test_quadratic<double>();
    /*test_cubic<double>();
    test_quartic<double>();

    test_constant<long double>();
    test_linear<long double>();
    test_quadratic<long double>();
    test_cubic<long double>();
    test_quartic<long double>();

#ifdef BOOST_HAS_FLOAT128
    test_constant<float128>();
    test_linear<float128>();
    test_quadratic<float128>();
    test_cubic<float128>();
    test_quartic<float128>();
#endif*/

    return boost::math::test::report_errors();
}
