//  (C) Copyright Matt Borland 2024.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
//  See: https://github.com/boostorg/math/issues/1203

#include "math_unit_test.hpp"
#include <boost/math/distributions/poisson.hpp>

int main()
{
    constexpr auto mean = 2719.3;

    const auto dist = boost::math::poisson_distribution<>(mean);

    // Generated from: Quantile[PoissonDistribution[2719.13], p]

    // Passing before issue (p >= .50):
    CHECK_ULP_CLOSE(boost::math::quantile(dist, 0.54), 2724.0, 10);
    CHECK_ULP_CLOSE(boost::math::quantile(dist, 0.53), 2723.0, 10);
    CHECK_ULP_CLOSE(boost::math::quantile(dist, 0.52), 2722.0, 10);
    CHECK_ULP_CLOSE(boost::math::quantile(dist, 0.51), 2720.0, 10);
    CHECK_ULP_CLOSE(boost::math::quantile(dist, 0.50), 2719.0, 10);

    // Failing before issue (p < .50):
    CHECK_ULP_CLOSE(boost::math::quantile(dist, 0.49), 2718.0, 10);
    CHECK_ULP_CLOSE(boost::math::quantile(dist, 0.48), 2716.0, 10);
    CHECK_ULP_CLOSE(boost::math::quantile(dist, 0.47), 2715.0, 10);

    CHECK_ULP_CLOSE(boost::math::quantile(dist, 0.04), 2628.0, 10);
    CHECK_ULP_CLOSE(boost::math::quantile(dist, 0.03), 2621.0, 10);
    CHECK_ULP_CLOSE(boost::math::quantile(dist, 0.02), 2613.0, 10);
    CHECK_ULP_CLOSE(boost::math::quantile(dist, 0.01), 2599.0, 10);

    // Check with smaller mean
    constexpr auto mean_2 = 12.0;
    const auto dist_2 = boost::math::poisson_distribution<>(mean_2);

    CHECK_ULP_CLOSE(boost::math::quantile(dist_2, 0.54), 12.0, 10);
    CHECK_ULP_CLOSE(boost::math::quantile(dist_2, 0.53), 12.0, 10);
    CHECK_ULP_CLOSE(boost::math::quantile(dist_2, 0.52), 12.0, 10);
    CHECK_ULP_CLOSE(boost::math::quantile(dist_2, 0.40), 11.0, 10);
    CHECK_ULP_CLOSE(boost::math::quantile(dist_2, 0.30), 10.0, 10);
    CHECK_ULP_CLOSE(boost::math::quantile(dist_2, 0.20), 09.0, 10);

    return boost::math::test::report_errors();
}
