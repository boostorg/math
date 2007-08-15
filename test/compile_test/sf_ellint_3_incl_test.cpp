//  Copyright John Maddock 2006.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// Basic sanity check that header <boost/math/special_functions/ellint_3.hpp>
// #includes all the files that it needs to.
//
#include <boost/math/special_functions/ellint_3.hpp>

template float boost::math::ellint_3<float>(float, float);
template double boost::math::ellint_3<double>(double, double);
template long double boost::math::ellint_3<long double>(long double, long double);

template float boost::math::ellint_3<float>(float, float, float);
template double boost::math::ellint_3<double>(double, double, double);
template long double boost::math::ellint_3<long double>(long double, long double, long double);
