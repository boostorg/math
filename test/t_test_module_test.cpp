//  (C) Copyright Nick Thompson 2019.
//  (C) Copyright Matt Borland 2021.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// Module build check for boost::math::statistics t-test entry points. The exported
// one_sample_t_test / two_sample_t_test / paired_samples_t_test forward to the
// detail impls, which pull mean_and_sample_variance from univariate_statistics and
// the cdf of students_t_distribution. Exercising them here forces those exported
// pieces to resolve from a module consumer. Expected values are copied from
// test_t_test.cpp; the remaining checks use mathematically-certain results.

#ifndef BOOST_MATH_BUILD_MODULE
#include <boost/math/statistics/t_test.hpp>
#else
import boost.math;
#endif

#ifndef BOOST_MATH_BUILD_MODULE
#include <vector>
#include <utility>
#endif
#include "math_unit_test.hpp"

// One-sample test, double precision. Statistic and p-value are the known-good
// values reproduced with Mathematica in test_agreement_with_mathematica.
void test_one_sample()
{
    std::vector<double> v {0.7304375676969546, 3.227250635039257, 1.01821954205186};

    const double expected_statistic {2.103013485037935};
    const double expected_pvalue {0.1701790440880712};

    const auto temp {boost::math::statistics::one_sample_t_test(v, 0.0)};
    CHECK_ULP_CLOSE(expected_statistic, std::get<0>(temp), 2);
    CHECK_ULP_CLOSE(expected_pvalue, std::get<1>(temp), 7);
}

// Two-sample test, double precision. A set compared against itself has zero
// mean difference, so the statistic is exactly 0 and the p-value exactly 1.
void test_two_sample()
{
    std::vector<double> set_1 {301, 298, 295, 297, 304, 305, 309, 298, 291, 299, 293, 304};

    const auto temp {boost::math::statistics::two_sample_t_test(set_1, set_1)};
    CHECK_ULP_CLOSE(0.0, std::get<0>(temp), 5);
    CHECK_ULP_CLOSE(1.0, std::get<1>(temp), 5);
}

// Paired-samples test, double precision. For {2,4} vs {1,2} the paired deltas
// are {1,2}, giving mean 1.5 and sample std dev sqrt(0.5); the statistic is
// exactly 1.5 / (sqrt(0.5)/sqrt(2)) = 3.
void test_paired_samples()
{
    std::vector<double> set_1 {2, 4};
    std::vector<double> set_2 {1, 2};

    const auto temp {boost::math::statistics::paired_samples_t_test(set_1, set_2)};
    CHECK_ULP_CLOSE(3.0, std::get<0>(temp), 5);
}

// Exercise the float instantiation of the exported paired-samples template. Same
// exact result (3) as the double case above.
void test_paired_samples_float()
{
    std::vector<float> set_1 {2, 4};
    std::vector<float> set_2 {1, 2};

    const auto temp {boost::math::statistics::paired_samples_t_test(set_1, set_2)};
    CHECK_ULP_CLOSE(3.0f, std::get<0>(temp), 5);
}

int main()
{
    test_one_sample();
    test_two_sample();
    test_paired_samples();
    test_paired_samples_float();

    return boost::math::test::report_errors();
}
