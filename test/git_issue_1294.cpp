// Copyright 2025 Matt Borland
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// See: https://github.com/boostorg/math/issues/1294

#include <boost/math/distributions/logistic.hpp>
#include "math_unit_test.hpp"

int main()
{
    using namespace boost::math::policies;
    using boost::math::logistic_distribution;

    typedef policy<
        promote_float<true>,
        promote_double<true>
    > with_promotion;

    constexpr double p = 2049.0/4096;
    constexpr double ref = 9.76562577610225755e-04;

    logistic_distribution<double, with_promotion> dist_promote;
    const double x = quantile(dist_promote, p);

    // Previously we had: 9.76562577610170027e-04
    // Which is an ULP distance of 256
    CHECK_ULP_CLOSE(x, ref, 1);

    return boost::math::test::report_errors();
}
