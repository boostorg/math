//  (C) Copyright Matt Borland 2021.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// Module build check for boost::math::statistics::one_sample_z_test and
// two_sample_z_test. It exercises the exported z-test entry points when the
// test is built against the boost.math module. The double/float expected
// values mirror the corresponding assertions in test_z_test.cpp.

#ifndef BOOST_MATH_BUILD_MODULE
#include <boost/math/statistics/z_test.hpp>
#else
import boost.math;
#endif

#ifndef BOOST_MATH_BUILD_MODULE
#include <limits>
#include <utility>
#include <vector>
#include <cmath>
#endif
#include "math_unit_test.hpp"

// Summary-statistic overload: statistic = (sample_mean - assumed_mean) / (sample_variance / sqrt(n)).
template <typename Real>
void test_one_sample_z()
{
    // Equal means give a zero statistic; the resulting p-value underflows to zero.
    std::pair<Real, Real> temp {boost::math::statistics::one_sample_z_test(static_cast<Real>(10), static_cast<Real>(2), static_cast<Real>(100), static_cast<Real>(10))};
    Real computed_statistic {std::get<0>(temp)};
    Real computed_pvalue {std::get<1>(temp)};
    CHECK_ULP_CLOSE(static_cast<Real>(0), computed_statistic, 5);
    CHECK_MOLLIFIED_CLOSE(static_cast<Real>(0), computed_pvalue, 5 * std::numeric_limits<Real>::epsilon());

    // (10 - 5) / (2 / sqrt(100)) = 25.
    temp = boost::math::statistics::one_sample_z_test(static_cast<Real>(10), static_cast<Real>(2), static_cast<Real>(100), static_cast<Real>(5));
    Real computed_statistic_2 {std::get<0>(temp)};
    CHECK_ULP_CLOSE(static_cast<Real>(25), computed_statistic_2, 5);

    // (1/2 - 1/3) / (10 / sqrt(100)) = 1/6.
    temp = boost::math::statistics::one_sample_z_test(static_cast<Real>(1) / 2, static_cast<Real>(10), static_cast<Real>(100), static_cast<Real>(1) / 3);
    Real computed_statistic_3 {std::get<0>(temp)};
    CHECK_ULP_CLOSE(static_cast<Real>(1) / 6, computed_statistic_3, 5);
}

// Container overload: two unit-shifted samples of equal variance give a statistic of exactly 1.
template <typename Real>
void test_two_sample_z()
{
    std::vector<Real> set_1 {1, 2, 3, 4, 5};
    std::vector<Real> set_2 {2, 3, 4, 5, 6};

    std::pair<Real, Real> temp {boost::math::statistics::two_sample_z_test(set_2, set_1)};
    Real computed_statistic {std::get<0>(temp)};
    Real computed_pvalue {std::get<1>(temp)};
    CHECK_ULP_CLOSE(static_cast<Real>(1), computed_statistic, 5);
    CHECK_MOLLIFIED_CLOSE(static_cast<Real>(0), computed_pvalue, std::sqrt(std::numeric_limits<Real>::epsilon()));
}

int main()
{
    test_one_sample_z<double>();
    test_two_sample_z<double>();

    test_one_sample_z<float>();
    test_two_sample_z<float>();

    return boost::math::test::report_errors();
}
