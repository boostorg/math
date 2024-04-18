//  Copyright Matt Borland 2024.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
//  See: https://github.com/scipy/scipy/issues/20124
//  See: https://github.com/boostorg/math/issues/1120

#include <boost/math/distributions/skew_normal.hpp>
#include "math_unit_test.hpp"
#include <limits>

int main()
{
    using scipy_policy = boost::math::policies::policy<boost::math::policies::promote_double<false>>;
    constexpr auto q = 0.012533469508013;

    boost::math::skew_normal_distribution<double, scipy_policy> dist(0, 1, 500);
    CHECK_ULP_CLOSE(0.01, boost::math::cdf(dist, q), 100);
    CHECK_ULP_CLOSE(q, boost::math::quantile(dist, 0.01), 100);

    return boost::math::test::report_errors();
}
