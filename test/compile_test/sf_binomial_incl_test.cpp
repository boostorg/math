//  Copyright John Maddock 2006.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// Basic sanity check that header <boost/math/special_functions/binomial.hpp>
// #includes all the files that it needs to.
//
#include <boost/math/special_functions/binomial.hpp>

template float boost::math::binomial_coefficient<float>(unsigned, unsigned);
template double boost::math::binomial_coefficient<double>(unsigned, unsigned);
template long double boost::math::binomial_coefficient<long double>(unsigned, unsigned);
