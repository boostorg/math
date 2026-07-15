// Copyright 2026 Matt Borland
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt
//
// Smoke test for the boost.math module: exercises at least one entity from
// every component exported by module/math.cxx.

#include <boost/core/lightweight_test.hpp>

// Import the STL if we can
#if defined(__cpp_lib_modules) && __cpp_lib_modules >= 202207L
import std;
#else
#include <iostream>
#include <cmath>
#endif

import boost.math;

int main()
{
    // constants
    BOOST_TEST_GT(boost::math::constants::pi<double>(), 3.14);
    BOOST_TEST_LT(boost::math::constants::pi<double>(), 3.15);
    BOOST_TEST_EQ(boost::math::double_constants::half, 0.5);
    static_assert(boost::math::float_constants::half == 0.5F, "float_constants");

    // policies
    const auto pol {boost::math::policies::make_policy(
        boost::math::policies::domain_error<boost::math::policies::errno_on_error>())};
    static_cast<void>(pol);
    static_assert(boost::math::policies::is_policy<decltype(pol)>::value, "is_policy");
    BOOST_TEST_GT((boost::math::policies::digits<double, boost::math::policies::policy<>>()), 50);

    // error handling types are visible through the import
    static_assert(boost::math::policies::is_noexcept_error_policy<
        boost::math::policies::policy<>>::value == false, "default policy throws");

    // special functions
    BOOST_TEST_EQ(boost::math::tgamma(4.0), 6.0);
    BOOST_TEST_EQ(boost::math::sign(-3.5), -1);
    BOOST_TEST(boost::math::isnan(std::sqrt(-1.0)));
    BOOST_TEST_EQ(boost::math::round(2.5), 3.0);
    BOOST_TEST_EQ(boost::math::factorial<double>(5), 120.0);
    BOOST_TEST_EQ(boost::math::fpclassify(1.0), FP_NORMAL);
    {
        const auto erf_value {boost::math::erf(1.0)};
        BOOST_TEST_GT(erf_value, 0.84);
        BOOST_TEST_LT(erf_value, 0.85);
        const auto bessel_value {boost::math::cyl_bessel_j(0, 1.0)};
        BOOST_TEST_GT(bessel_value, 0.76);
        BOOST_TEST_LT(bessel_value, 0.77);
    }

    // distributions
    {
        const boost::math::normal_distribution<> dist {};
        BOOST_TEST_EQ(boost::math::cdf(dist, 0.0), 0.5);
        BOOST_TEST_EQ(boost::math::median(dist), 0.0);
        BOOST_TEST_EQ(boost::math::cdf(boost::math::complement(dist, 0.0)), 0.5);
        const boost::math::students_t st {5.0};
        BOOST_TEST_GT(boost::math::pdf(st, 0.0), 0.37);
        const auto q {boost::math::quantile(boost::math::binomial(50, 0.5), 0.5)};
        BOOST_TEST_GT(q, 24.0);
        BOOST_TEST_LT(q, 26.0);
    }

    return boost::report_errors();
}
