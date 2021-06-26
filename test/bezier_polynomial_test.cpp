/*
 * Copyright Nick Thompson, 2021
 * Use, modification and distribution are subject to the
 * Boost Software License, Version 1.0. (See accompanying file
 * LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
 */

#include "math_unit_test.hpp"
#include <numeric>
#include <random>
#include <array>
#include <boost/core/demangle.hpp>
#include <boost/math/interpolators/bezier_polynomial.hpp>
#ifdef BOOST_HAS_FLOAT128
#include <boost/multiprecision/float128.hpp>
using boost::multiprecision::float128;
#endif

using boost::math::interpolators::bezier_polynomial;

template<typename Real>
void test_linear()
{
    std::vector<std::array<Real, 2>> control_points(2);
    control_points[0] = {0.0, 0.0};
    control_points[1] = {1.0, 1.0};
    auto control_points_copy = control_points;
    auto bp = bezier_polynomial(std::move(control_points_copy));

    CHECK_ULP_CLOSE(bp(0)[0], control_points[0][0], 3);
    CHECK_ULP_CLOSE(bp(0)[0], control_points[0][1], 3);
}

int main()
{
    test_linear<double>();
    return boost::math::test::report_errors();
}
