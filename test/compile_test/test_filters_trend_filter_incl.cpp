//  Copyright Nick Thompson 2026.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// Basic sanity check that header <boost/math/filters/trend_filter.hpp>
// #includes all the files that it needs to.
//
#include <boost/math/filters/trend_filter.hpp>
#include <vector>

void compile_and_link_test()
{
   using boost::math::filters::trend_filter;
   const std::vector<double> y{1, 2, 4, 3, 5, 6};
   std::vector<double> x = trend_filter(y, 1, 0.5);
   x = trend_filter(y);
   const std::vector<float> yf{1, 2, 4, 3, 5, 6};
   std::vector<float> xf = trend_filter(yf, 0, 0.5f);
   xf = trend_filter(yf);
}
