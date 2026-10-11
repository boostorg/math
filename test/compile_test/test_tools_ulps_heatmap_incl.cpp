//  Copyright Nick Thompson 2026.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// Basic sanity check that header
// #includes all the files that it needs to.
//
#include <boost/math/tools/ulps_heatmap.hpp>
//
// Note this header includes no other headers, this is
// important if this test is to be meaningful:
//
#include "test_compile_result.hpp"

void compile_and_link_test()
{
    auto f = [](double x, double y) { return x*y; };
    auto g = [](float x, float y) { return x*y; };
    boost::math::tools::ulps_heatmap<decltype(f), double, float> test(f, 1.0f, 2.0f, 1.0f, 2.0f, 4, 4);
    auto dummy = test.add_fn(g, "g").summary();
}
