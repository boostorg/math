// Copyright Matt Borland, 2022
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0. (See accompanying file
// LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#define BOOST_MATH_OVERFLOW_ERROR_POLICY ignore_error

#ifndef BOOST_MATH_BUILD_MODULE
#include <boost/math/special_functions/powm1.hpp>
#else
import boost.math;
#endif

#ifndef BOOST_MATH_BUILD_MODULE
#include <cmath>
#endif
#include "math_unit_test.hpp"

template <typename T>
void test()
{
    CHECK_EQUAL(std::pow(T(0), T(2)) - 1, boost::math::powm1(T(0), T(2)));
    CHECK_EQUAL(std::pow(T(0), T(0.1)) - 1, boost::math::powm1(T(0), T(0.1)));

    // These return infinity only under the ignore_error overflow policy; the macro
    // cannot reach the precompiled module (which throws), so keep them textual.
    #ifndef BOOST_MATH_BUILD_MODULE
    CHECK_EQUAL(std::pow(T(0), T(-2)) - 1, boost::math::powm1(T(0), T(-2)));
    CHECK_EQUAL(std::pow(T(0), T(-0.1)) - 1, boost::math::powm1(T(0), T(-0.1)));
    #endif
}

int main()
{
    test<float>();
    test<double>();

    // long double is excluded from the module build (double and float only).
    #if !defined(BOOST_MATH_BUILD_MODULE) && !defined(BOOST_MATH_NO_LONG_DOUBLE_MATH_FUNCTIONS)
    test<long double>();
    #endif

    return boost::math::test::report_errors();
}
