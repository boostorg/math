//  (C) Copyright John Maddock 2006.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_SQRT1PM1
#define BOOST_MATH_SQRT1PM1

#include <boost/math/special_functions/log1p.hpp>
#include <boost/math/special_functions/expm1.hpp>

//
// This algorithm computes sqrt(1+x)-1 for small x:
//

namespace boost{ namespace math{

template <class T>
inline T sqrt1pm1(const T& val)
{
   using namespace std;

   if(fabs(val) > 0.75)
      return sqrt(1 + val) - 1;
   return boost::math::expm1(boost::math::log1p(val) / 2);
}

} // namespace math
} // namespace boost

#endif // BOOST_MATH_SQRT1PM1



