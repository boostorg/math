//  (C) Copyright John Maddock 2005.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_HYPOT_INCLUDED
#define BOOST_MATH_HYPOT_INCLUDED

#include <cmath>
#include <boost/limits.hpp>
#include <algorithm> // swap

#ifndef BOOST_NO_LIMITS_COMPILE_TIME_CONSTANTS
#  include <boost/static_assert.hpp>
#else
#  include <boost/assert.hpp>
#endif

#if BOOST_WORKAROUND(__BORLANDC__, BOOST_TESTED_AT(0x580))
#include <boost/type_traits/remove_const.hpp>
#define BOOST_RCT(T) typename remove_const<T>::type
#else
#define BOOST_RCT(T) T
#endif

#ifdef BOOST_NO_STDC_NAMESPACE
namespace std{ using ::sqrt; using ::fabs; }
#endif


namespace boost{ namespace math{

template <class T>
T hypot(T x, T y)
{
#ifndef BOOST_NO_LIMITS_COMPILE_TIME_CONSTANTS
   BOOST_STATIC_ASSERT(::std::numeric_limits<BOOST_RCT(T)>::is_specialized);
#else
   BOOST_ASSERT(std::numeric_limits<BOOST_RCT(T)>::is_specialized);
#endif

   //
   // normalize x and y, so that both are positive and x >= y:
   //
   x = (std::fabs)(x);
   y = (std::fabs)(y);

   // special case, see C99 Annex F:
   if(std::numeric_limits<BOOST_RCT(T)>::has_infinity
      && ((x == std::numeric_limits<BOOST_RCT(T)>::infinity())
      || (y == std::numeric_limits<BOOST_RCT(T)>::infinity())))
      return std::numeric_limits<BOOST_RCT(T)>::infinity();

   if(y > x) 
      (std::swap)(x, y);
   //
   // figure out overflow and underflow limits:
   //
   T safe_upper = (std::sqrt)((std::numeric_limits<BOOST_RCT(T)>::max)()) / 2;
   T safe_lower = (std::sqrt)((std::numeric_limits<BOOST_RCT(T)>::min)());
   static const T one = 1;
   //
   // Now handle special cases:
   //
   if(x >= safe_upper)
   {
      if(y <= one)
      {
         // y is neligible:
         return x;
      }
      return (std::sqrt)(x) * (std::sqrt)(y) * (std::sqrt)(x/y + y/x);
   }
   else if(y <= safe_lower)
   {
      if((x >= one) || (y == 0))
      {
         // y is negligible:
         return x;
      }
      return (std::sqrt)(x) * (std::sqrt)(y) * (std::sqrt)(x/y + y/x);
   }
   //
   // If we get here then x^2+y^2 will not overflow or underflow:
   //
   return (std::sqrt)(x*x + y*y);
}

} } // namespaces

#undef BOOST_RCT

#endif // BOOST_MATH_HYPOT_INCLUDED
