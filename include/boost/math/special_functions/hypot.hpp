//  (C) Copyright John Maddock 2005-2006.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_HYPOT_INCLUDED
#define BOOST_MATH_HYPOT_INCLUDED

#include <boost/math/tools/config.hpp>
#include <boost/math/tools/precision.hpp>
#include <cmath>
#include <algorithm> // for swap

#ifdef BOOST_NO_STDC_NAMESPACE
namespace std{ using ::sqrt; using ::fabs; }
#endif

namespace boost{ namespace math{

template <class T>
T hypot(T x, T y)
{
   //
   // Normalize x and y, so that both are positive and x >= y:
   //
   using namespace std; // ADL of std names

   x = fabs(x);
   y = fabs(y);

#ifdef BOOST_MSVC
#pragma warning(push) 
#pragma warning(disable: 4127)
#endif
   // special case, see C99 Annex F:
   if(std::numeric_limits<T>::has_infinity
      && ((x == std::numeric_limits<T>::infinity())
      || (y == std::numeric_limits<T>::infinity())))
      return std::numeric_limits<T>::infinity();
#ifdef BOOST_MSVC
#pragma warning(pop)
#endif

   if(y > x)
      (std::swap)(x, y);
   //
   // Figure out overflow and underflow limits,
   // we could make these constants static to save
   // a few cycles, but the code would then not be
   // thread safe :-(
   //
   T safe_upper = sqrt(tools::max_value<T>()) / 2;
   T safe_lower = sqrt(tools::min_value<T>());
   static const T one = 1;
   //
   // Now handle special cases:
   //
   if(x >= safe_upper)
   {
      if(y <= one)
      {
         // y is negligible:
         return x;
      }
      T a = sqrt(x) * sqrt(y);
      T b = sqrt(x/y + y/x);
      // We may have overflow:
      if(tools::max_value<T>() /a < b)
         return tools::overflow_error<T>(BOOST_CURRENT_FUNCTION, 0);
      return a * b;
   }
   else if(y <= safe_lower)
   {
      if((x >= one) || (y == 0))
      {
         // y is negligible:
         return x;
      }
      return sqrt(x) * sqrt(y) * sqrt(x/y + y/x);
   }
   //
   // If we get here then x^2+y^2 will not overflow or underflow:
   //
   return sqrt(x*x + y*y);
} // template <class T> T hypot(T x, T y)


} // namespace math
} // namespace boost

#endif // BOOST_MATH_HYPOT_INCLUDED

