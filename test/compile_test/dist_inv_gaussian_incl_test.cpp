//  Copyright John Maddock 2021.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// Basic sanity check that header <boost/math/distributions/inverse_gaussian.hpp>
// #includes all the files that it needs to.
//
#include <boost/math/distributions/inverse_gaussian.hpp>
//
// Note this header includes no other headers, this is
// important if this test is to be meaningful:
//
#include "test_compile_result.hpp"

void compile_and_link_test()
{
   TEST_DIST_FUNC(inverse_gaussian)
}

template class boost::math::inverse_gaussian_distribution<float, boost::math::policies::policy<> >;
template class boost::math::inverse_gaussian_distribution<double, boost::math::policies::policy<> >;
#ifndef BOOST_MATH_NO_LONG_DOUBLE_MATH_FUNCTIONS
template class boost::math::inverse_gaussian_distribution<long double, boost::math::policies::policy<> >;
#endif
