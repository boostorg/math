// Copyright Matt Borland, 2023
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0. (See accompanying file
// LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// See: https://godbolt.org/z/Ev4ManrsW

#include <boost/math/special_functions/round.hpp>
#include <cstdint>
#include "math_unit_test.hpp"

double x = 9223372036854775807.0;  // can't be represented as double, will have a different value at runtime. 
int main()
{
    int64_t result = boost::math::llround(x);
    CHECK_EQUAL(result, INT64_C(9223372036854775807));

    return boost::math::test::report_errors();
}
