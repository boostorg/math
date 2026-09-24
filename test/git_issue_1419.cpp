//  Copyright Matt Borland 2026
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
//  See: https://github.com/boostorg/math/issues/1419

#include <boost/math/statistics/univariate_statistics.hpp>
#include <algorithm>
#include <vector>
#include "math_unit_test.hpp"

int main()
{
    std::vector<int> v = {3, 2, 1, 2};
    const auto modes_unsorted = boost::math::statistics::mode(v);

    std::sort(v.begin(), v.end());
    const auto modes_sorted = boost::math::statistics::mode(v);

    CHECK_EQUAL(modes_unsorted.size(), modes_sorted.size());
    CHECK_EQUAL(modes_unsorted.front(), modes_sorted.front());
    CHECK_EQUAL(modes_unsorted.front(), 2);

    return boost::math::test::report_errors();
}
