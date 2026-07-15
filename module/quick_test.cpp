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

    return boost::report_errors();
}
