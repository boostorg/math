//  (C) Copyright John Maddock 2006.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_TOOLS_SIGN_HPP
#define BOOST_MATH_TOOLS_SIGN_HPP

#include <cmath>
#include <cstdlib>

namespace boost{ namespace math{ 

template <class T>
inline int sign(const T& z)
{
   return (z == 0) ? 0 : (z < 0) ? -1 : 1;
}

template <class T>
inline int signbit(const T& z)
{
   return (z < 0) ? 1 : 0;
}

template <class T>
inline T copysign(const T& x, const T& y)
{
   return fabs(x) * boost::math::sign(y);
}

} // namespace math
} // namespace boost


#endif // BOOST_MATH_TOOLS_SIGN_HPP

