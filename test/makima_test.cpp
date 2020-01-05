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
#include <boost/math/interpolators/makima.hpp>
#ifdef BOOST_HAS_FLOAT128
#include <boost/multiprecision/float128.hpp>
using boost::multiprecision::float128;
#endif


using boost::math::interpolators::makima;

template<typename Real>
void test_constant()
{

    std::vector<Real> x{0,1,2,3, 9, 22, 81};
    std::vector<Real> y(x.size());
    for (auto & t : y) {
        t = 7;
    }

    auto x_copy = x;
    auto y_copy = y;
    auto akima = makima(std::move(x_copy), std::move(y_copy));

    std::cout << "Constant Makima: " << akima << "\n";
    for (Real t = x[0]; t <= x[x.size()-1]; t += 0.25) {
        //CHECK_ULP_CLOSE(7, akima(t), 2);
    }
}

template<typename Real>
void test_linear()
{
    std::vector<Real> x{0,1,2};
    std::vector<Real> y{0,1,2};

    auto x_copy = x;
    auto y_copy = y;
    auto akima = makima(std::move(x_copy), std::move(y_copy));

    CHECK_ULP_CLOSE(y[0], akima(x[0]), 0);
    CHECK_ULP_CLOSE(Real(1)/Real(2), akima(Real(1)/Real(2)), 10);
    CHECK_ULP_CLOSE(y[1], akima(x[1]), 0);
    CHECK_ULP_CLOSE(Real(3)/Real(2), akima(Real(3)/Real(2)), 10);
    CHECK_ULP_CLOSE(y[2], akima(x[2]), 0);

    x.resize(45);
    y.resize(45);
    for (size_t i = 0; i < x.size(); ++i) {
        x[i] = i;
        y[i] = i;
    }

    x_copy = x;
    y_copy = y;
    akima = makima(std::move(x_copy), std::move(y_copy));
    std::cout << "Akima = " << akima << "\n";
    for (Real t = 0; t < x.size() - 1; t += 0.5) {
        std::cout << "t = " << t << ", akima(t) = " << akima(t) << "\n";
       // CHECK_ULP_CLOSE(t, akima(t), 0);
    }
}

template<typename Real>
void test_quadratic()
{
    std::vector<Real> x(25);
    std::vector<Real> y(25);
    for(size_t i = 0; i < x.size(); ++i) {
        x[i] = i;
        y[i] = Real(i*i)/Real(2);
    }

    auto x_copy = x;
    auto y_copy = y;
    auto akima = makima(std::move(x_copy), std::move(y_copy));

    std::cout << "Quadratic Akima: " << akima << "\n";
}

int main()
{
    //test_constant<double>();
    test_linear<double>();
    test_quadratic<double>();
    return boost::math::test::report_errors();
}
