// Copyright Matt Borland, 2023
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0. (See accompanying file
// LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// See: https://github.com/scipy/scipy/issues/14901

#ifndef BOOST_MATH_BUILD_MODULE
#include <boost/math/distributions/non_central_f.hpp>
#else
import boost.math;
#endif

// FE_* are macros, so they must come from a textual C header even in module
// mode (import std cannot export macros); <fenv.h> is safe to leave unguarded.
#include <fenv.h>

// std::fetestexcept comes from import std in module mode; the non-module build
// gets it from the <cfenv> wrapper.
#ifndef BOOST_MATH_BUILD_MODULE
#include <cfenv>
#endif

#include "math_unit_test.hpp"

#pragma STDC FENV_ACCESS ON

int main()
{
    auto ncf1 = boost::math::non_central_f(1, 1, 1);
    const auto cdf1 = boost::math::cdf(ncf1, 2);

    CHECK_ULP_CLOSE(cdf1, 1 - 0.5230393170884924, 20);

    if (std::fetestexcept(FE_INVALID) || std::fetestexcept(FE_DIVBYZERO))
    {
        return 1;
    }

    return boost::math::test::report_errors();
}
