//  (C) Copyright Matt Borland 2022.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <boost/math/tools/color_maps.hpp>
#include "test_compile_result.hpp"

void compile_and_link_test()
{
    boost::math::tools::smooth_cool_warm_color_map cw;
    check_result<double>(cw(0.5));

    boost::math::tools::plasma_color_map<float> plasma;
    check_result<float>(plasma(0.5f));

    boost::math::tools::viridis_color_map<double> viridis;
    check_result<double>(viridis(0.5));

#ifndef BOOST_MATH_NO_LONG_DOUBLE_MATH_FUNCTIONS
    boost::math::tools::inferno_color_map<long double> inferno;
    check_result<long double>(inferno(0.5l));
#endif
}
