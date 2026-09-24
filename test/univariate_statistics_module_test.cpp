/*
 *  (C) Copyright Nick Thompson 2018.
 *  (C) Copyright Matt Borland 2020.
 *  Use, modification and distribution are subject to the
 *  Boost Software License, Version 1.0. (See accompanying file
 *  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
 *
 *  Module build check for boost::math::statistics univariate statistics. The
 *  component is consumed through `import boost.math;` and every expected value
 *  below is copied from the non-module univariate_statistics_test.cpp. Only the
 *  container entry points are exercised (double and float).
 */

#ifndef BOOST_MATH_BUILD_MODULE
#include <boost/math/statistics/univariate_statistics.hpp>
#else
import boost.math;
#endif

#include "math_unit_test.hpp"

#ifndef BOOST_MATH_BUILD_MODULE
#include <vector>
#include <limits>
#include <tuple>
#endif

// Exercises the main public entry points on small, hand-checkable samples.
template <class Real>
void test_univariate_statistics()
{
    using namespace boost::math::statistics;

    const Real tol {std::numeric_limits<Real>::epsilon()};

    // mean of {1,2,3,4,5} is 3; mean of {1,2,3} is 2.
    std::vector<Real> v {1, 2, 3, 4, 5};
    CHECK_ULP_CLOSE(3, mean(v), 4);
    std::vector<Real> v3 {1, 2, 3};
    CHECK_ULP_CLOSE(2, mean(v3), 4);

    // variance of a constant sample is zero.
    std::vector<Real> c {1, 1, 1, 1, 1, 1};
    CHECK_ABSOLUTE_ERROR(0, variance(c), 20 * tol);
    CHECK_ABSOLUTE_ERROR(0, sample_variance(c), 20 * tol);

    // variance of the alternating sample {0,1,0,1,...} is 1/4.
    std::vector<Real> w {0, 1, 0, 1, 0, 1, 0, 1};
    CHECK_ABSOLUTE_ERROR(0.25, variance(w), 20 * tol);

    // {0,0,0,0,5}: mu = 1, sigma^2 = 4, skewness = 3/2, kurtosis = 13/4.
    std::vector<Real> s {0, 0, 0, 0, 5};
    CHECK_ABSOLUTE_ERROR(1.5, skewness(s), 20 * tol);
    CHECK_ABSOLUTE_ERROR(3.25, kurtosis(s), 20 * tol);

    // {1,2,3,4,5}: kurtosis = 17/10 and the first four moments are (3,2,0,34/5).
    CHECK_ABSOLUTE_ERROR(1.7, kurtosis(v), 20 * tol);
    auto [M1, M2, M3, M4] = first_four_moments(v);
    CHECK_ABSOLUTE_ERROR(3, M1, 20 * tol);
    CHECK_ABSOLUTE_ERROR(2, M2, 20 * tol);
    CHECK_ABSOLUTE_ERROR(0, M3, 20 * tol);
    CHECK_ABSOLUTE_ERROR(6.8, M4, 20 * tol);

    // median extracts (or averages) exact sample values.
    std::vector<Real> m_odd {1, 2, 3, 4, 5, 6, 7};
    CHECK_EQUAL(static_cast<Real>(4), median(m_odd));
    std::vector<Real> m_even {2, 4};
    CHECK_EQUAL(static_cast<Real>(3), median(m_even));

    // Gini coefficient of {1,0,0} is 2/3.
    std::vector<Real> g {1, 0, 0};
    CHECK_ABSOLUTE_ERROR(2.0 / 3.0, gini_coefficient(g), 20 * tol);

    // Interquartile range: Wikipedia's worked example gives 88; {1,..,5} gives 3.
    std::vector<Real> iqr_data {7, 7, 31, 31, 47, 75, 87, 115, 116, 119, 119, 155, 177};
    CHECK_EQUAL(static_cast<Real>(88), interquartile_range(iqr_data));
    std::vector<Real> iqr_small {1, 2, 3, 4, 5};
    CHECK_EQUAL(static_cast<Real>(3), interquartile_range(iqr_small));
}

int main()
{
    test_univariate_statistics<double>();
    test_univariate_statistics<float>();

    return boost::math::test::report_errors();
}
