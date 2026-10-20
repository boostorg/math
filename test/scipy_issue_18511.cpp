// Copyright Matt Borland, 2023
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0. (See accompanying file
// LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// https://github.com/scipy/scipy/issues/18511

#ifndef BOOST_MATH_BUILD_MODULE
#include <boost/math/distributions/hypergeometric.hpp>
#include <boost/math/policies/policy.hpp>
#include <limits>
#include <cstdint>
#else
import boost.math;
// UINT64_C is a preprocessor macro and is not provided by import std, so pull
// it from the pure-C header (not the <cstdint> wrapper) in module mode.
#include <stdint.h>
#endif

#include "math_unit_test.hpp"

template <typename T>
void test()
{
    //https://www.wolframalpha.com/input?i=%2812000+*+370000%29+%2F+390000
    auto dist = boost::math::hypergeometric_distribution<T>(UINT64_C(12'000), UINT64_C(370'000), UINT64_C(390'000));
    auto hm = boost::math::mean(dist);
    CHECK_ULP_CLOSE(hm, static_cast<T>(11384.615384615384615384615384615384615384615384615384615384615384L), 1);

    // Same result as above because both the numerator and denom have been multiplied by 100 (Exceeds UINT32_MAX)
    auto dist2 = boost::math::hypergeometric_distribution<T>(UINT64_C(12'000), UINT64_C(37'000'000), UINT64_C(39'000'000)); 
    auto hm2 = boost::math::mean(dist2);
    CHECK_ULP_CLOSE(hm2, static_cast<T>(11384.615384615384615384615384615384615384615384615384615384615384L), 1);
}

int main()
{
    test<float>();
    test<double>();

    return boost::math::test::report_errors();
}
