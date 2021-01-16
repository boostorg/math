//  (C) Copyright Matt Borland 2021.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include "math_unit_test.hpp"
#include <boost/math/statistics/z_test.hpp>
#include <boost/math/statistics/univariate_statistics.hpp>
#include <boost/multiprecision/cpp_bin_float.hpp>
#include <limits>
#include <vector>
#include <random>

using quad = boost::multiprecision::cpp_bin_float_quad;

template<typename Real>
void test_one_sample_z()
{
    auto [computed_statistic, computed_pvalue] = boost::math::statistics::one_sample_z_test(Real(10), Real(2), Real(100), Real(10));
    CHECK_ULP_CLOSE(Real(0), computed_statistic, 5);
    CHECK_MOLLIFIED_CLOSE(Real(0), computed_pvalue, 5*std::numeric_limits<Real>::epsilon());

    auto [computed_statistic_2, computed_pvalue_2] = boost::math::statistics::one_sample_z_test(Real(10), Real(2), Real(100), Real(5));
    CHECK_ULP_CLOSE(Real(25), computed_statistic_2, 5);

    auto [computed_statistic_3, computed_pvalue_3] = boost::math::statistics::one_sample_z_test(Real(1.0/2.0), Real(10), Real(100), Real(1.0/3.0));
    CHECK_ULP_CLOSE(Real(1)/6, computed_statistic_3, 5);
}

template<typename Z>
void test_integer_one_sample_z()
{
    auto [computed_statistic, computed_pvalue] = boost::math::statistics::one_sample_z_test(Z(10), Z(2), Z(100), Z(10));
    CHECK_ULP_CLOSE(0, computed_statistic, 5);
    CHECK_MOLLIFIED_CLOSE(0, computed_pvalue, 5*std::numeric_limits<double>::epsilon());

    auto [computed_statistic_2, computed_pvalue_2] = boost::math::statistics::one_sample_z_test(Z(10), Z(2), Z(100), Z(5));
    CHECK_ULP_CLOSE(25, computed_statistic_2, 5);
}

template<typename Real>
void test_two_sample_z()
{
    std::vector<Real> set_1 {1,2,3,4,5};
    std::vector<Real> set_2 {2,3,4,5,6};

    auto [computed_statistic, computed_pvalue] = boost::math::statistics::two_sample_z_test(set_2, set_1);
    CHECK_ULP_CLOSE(Real(1), computed_statistic, 5);
    CHECK_MOLLIFIED_CLOSE(Real(0), computed_pvalue, 5*std::numeric_limits<Real>::epsilon());
}

template<typename Z>
void test_integer_two_sample_z()
{
    std::vector<Z> set_1 {1,2,3,4,5};
    std::vector<Z> set_2 {2,3,4,5,6};

    auto [computed_statistic, computed_pvalue] = boost::math::statistics::two_sample_z_test(set_2, set_1);
    CHECK_ULP_CLOSE(1, computed_statistic, 5);
    CHECK_MOLLIFIED_CLOSE(0, computed_pvalue, 5*std::numeric_limits<double>::epsilon());
}

int main()
{
    test_one_sample_z<float>();
    test_one_sample_z<double>();
    test_one_sample_z<quad>();

    test_integer_one_sample_z<int>();
    test_integer_one_sample_z<int32_t>();
    test_integer_one_sample_z<int64_t>();
    test_integer_one_sample_z<uint32_t>();

    test_two_sample_z<float>();
    test_two_sample_z<double>();
    test_two_sample_z<quad>();

    test_integer_two_sample_z<int>();
    test_integer_two_sample_z<int32_t>();
    test_integer_two_sample_z<int64_t>();
    test_integer_two_sample_z<uint32_t>();
}
