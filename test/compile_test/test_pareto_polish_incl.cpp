//  Copyright Nick Thompson 2026.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// Basic sanity check that header
// #includes all the files that it needs to.
//
#include <boost/math/optimization/pareto_polish.hpp>
//
// Note this header includes no other headers, this is
// important if this test is to be meaningful:
//
#include "test_compile_result.hpp"

void compile_and_link_test()
{
    auto func = [](auto const & v) { return std::vector<typename std::decay<decltype(v[0])>::type>{v[0]*v[0], (v[0] - 1)*(v[0] - 1)}; };
    boost::math::optimization::pareto_polish_parameters<std::vector<double>> params;
    params.lower_bounds = {0.0};
    params.upper_bounds = {1.0};
    std::vector<std::pair<std::vector<double>, std::vector<double>>> front{{{0.5}, {0.25, 0.25}}};
    auto v = boost::math::optimization::pareto_polish(func, front, params);
}
