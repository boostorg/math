// Copyright Matt Borland, 2023
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// 
// See: https://github.com/boostorg/math/issues/935

#include "math_unit_test.hpp"
#include <boost/math/distributions/negative_binomial.hpp>
#include <cmath>

int main (void)
{
    auto dist = boost::math::negative_binomial(5.0, 0.5);
    constexpr double supp = 6;

    double cdf_supp = boost::math::cdf(dist, supp);
    CHECK_ULP_CLOSE(0.7255859375, cdf_supp, 1);

    double cdf_supp0 = cdf_supp - 10.0 * (std::nextafter(cdf_supp, 1.0) - cdf_supp);
    CHECK_ULP_CLOSE(0.7255859374999989, cdf_supp0, 1);

    double quantile_supp = boost::math::quantile(dist, cdf_supp0);
    CHECK_ULP_CLOSE(quantile_supp, supp, 1); // Both should be 6

    return boost::math::test::report_errors();
}
