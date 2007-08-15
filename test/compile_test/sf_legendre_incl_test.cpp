//  Copyright John Maddock 2006.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// Basic sanity check that header <boost/math/special_functions/legendre.hpp>
// #includes all the files that it needs to.
//
#include <boost/math/special_functions/legendre.hpp>

template float boost::math::legendre_p<float>(int, float);
template double boost::math::legendre_p<double>(int, double);
template long double boost::math::legendre_p<long double>(int, long double);

template float boost::math::legendre_p<float>(int, int, float);
template double boost::math::legendre_p<double>(int, int, double);
template long double boost::math::legendre_p<long double>(int, int, long double);

template float boost::math::legendre_q<float>(unsigned, float);
template double boost::math::legendre_q<double>(unsigned, double);
template long double boost::math::legendre_q<long double>(unsigned, long double);

