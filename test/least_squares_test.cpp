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

using boost::math::statistics::ordinary_least_squares;

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

    auto [computed_c0, computed_c1] = ordinary_least_squares(x, y);

    CHECK_ULP_CLOSE(expected_c0, computed_c0, 10);
    CHECK_ULP_CLOSE(expected_c1, computed_c1, 10);
}

int main()
{
    test_line<double>();
    return boost::math::test::report_errors();
}
