//  Copyright John Maddock 2006.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// Basic sanity check that header <boost/math/special_functions/factorials.hpp>
// #includes all the files that it needs to.
//
#include <boost/math/special_functions/factorials.hpp>

template float boost::math::factorial<float>(unsigned);
template double boost::math::factorial<double>(unsigned);
template long double boost::math::factorial<long double>(unsigned);

template float boost::math::double_factorial<float>(unsigned);
template double boost::math::double_factorial<double>(unsigned);
template long double boost::math::double_factorial<long double>(unsigned);

template float boost::math::rising_factorial<float>(float, int);
template double boost::math::rising_factorial<double>(double, int);
template long double boost::math::rising_factorial<long double>(long double, int);

template float boost::math::falling_factorial<float>(float, unsigned int);
template double boost::math::falling_factorial<double>(double, unsigned int);
template long double boost::math::falling_factorial<long double>(long double, unsigned int);


