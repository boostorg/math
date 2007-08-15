//  Copyright John Maddock 2006.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// Basic sanity check that header <boost/math/special_functions/sign.hpp>
// #includes all the files that it needs to.
//
#include <boost/math/special_functions/sign.hpp>

template int boost::math::sign<float>(const float&);
template int boost::math::sign<double>(const double&);
template int boost::math::sign<long double>(const long double&);

template int boost::math::signbit<float>(const float&);
template int boost::math::signbit<double>(const double&);
template int boost::math::signbit<long double>(const long double&);

template float boost::math::copysign<float>(const float&, const float&);
template double boost::math::copysign<double>(const double&, const double&);
template long double boost::math::copysign<long double>(const long double&, const long double&);

