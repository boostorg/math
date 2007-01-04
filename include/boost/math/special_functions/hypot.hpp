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

   if(x * tools::epsilon<T>() >= y)
      return x;

   T rat = y / x;
   return x * sqrt(1 + rat*rat);
} // template <class T> T hypot(T x, T y)


} // namespace math
} // namespace boost

#endif // BOOST_MATH_HYPOT_INCLUDED

