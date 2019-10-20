/*
 * Copyright Nick Thompson, 2019
 * Use, modification and distribution are subject to the
 * Boost Software License, Version 1.0. (See accompanying file
 * LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
 */

#include "math_unit_test.hpp"
#include <vector>
#include <random>
#include <boost/math/statistics/linear_regression.hpp>

using boost::math::statistics::simple_ordinary_least_squares;
using boost::math::statistics::simple_ordinary_least_squares_with_R_squared;

template<typename Real>
void test_line()
{
    std::vector<Real> x(128);
    std::vector<Real> y(128);
    Real expected_c0 = 7;
    Real expected_c1 = 12;
    for (size_t i = 0; i < x.size(); ++i) {
        x[i] = i;
        y[i] = expected_c0 + expected_c1*x[i];
    }

    auto [computed_c0, computed_c1] = simple_ordinary_least_squares(x, y);

    CHECK_ULP_CLOSE(expected_c0, computed_c0, 0);
    CHECK_ULP_CLOSE(expected_c1, computed_c1, 0);

    auto [computed_c0_R, computed_c1_R, Rsquared] = simple_ordinary_least_squares_with_R_squared(x, y);

    Real expected_Rsquared = 1;
    CHECK_ULP_CLOSE(expected_c0, computed_c0, 0);
    CHECK_ULP_CLOSE(expected_c1, computed_c1, 0);
    CHECK_ULP_CLOSE(expected_Rsquared, Rsquared, 0);

}

template<typename Real>
void test_constant()
{
    std::vector<Real> x(128);
    std::vector<Real> y(128);
    Real expected_c0 = 7;
    Real expected_c1 = 0;
    for (size_t i = 0; i < x.size(); ++i) {
        x[i] = i;
        y[i] = expected_c0 + expected_c1*x[i];
    }

    auto [computed_c0, computed_c1] = simple_ordinary_least_squares(x, y);

    CHECK_ULP_CLOSE(expected_c0, computed_c0, 0);
    CHECK_ULP_CLOSE(expected_c1, computed_c1, 0);

    auto [computed_c0_R, computed_c1_R, Rsquared] = simple_ordinary_least_squares_with_R_squared(x, y);

    Real expected_Rsquared = 1;
    CHECK_ULP_CLOSE(expected_c0, computed_c0, 0);
    CHECK_ULP_CLOSE(expected_c1, computed_c1, 0);
    CHECK_ULP_CLOSE(expected_Rsquared, Rsquared, 0);

}


int main()
{
    test_line<float>();
    test_line<double>();
    test_line<long double>();

    test_constant<float>();
    test_constant<double>();
    test_constant<long double>();
    return boost::math::test::report_errors();
}
