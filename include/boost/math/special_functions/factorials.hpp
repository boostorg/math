//  Copyright John Maddock 2006.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_SP_FACTORIALS_HPP
#define BOOST_MATH_SP_FACTORIALS_HPP

#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/special_functions/detail/unchecked_factorial.hpp>
#include <boost/array.hpp>
#ifdef BOOST_MSVC
#pragma warning(push) // Temporary until lexical cast fixed.
#pragma warning(disable: 4127 4701)
#endif
#include <boost/lexical_cast.hpp>
#ifdef BOOST_MSVC
#pragma warning(pop)
#endif
#include <cmath>

namespace boost { namespace math
{

template <class T>
T factorial(unsigned i)
{
   using namespace std; // Aid ADL for floor.

   if(i <= max_factorial<T>::value)
      return unchecked_factorial<T>(i);
   T result = boost::math::tgamma(static_cast<T>(i+1));
   if(result > tools::max_value<T>())
      return result; // Overflowed value! (But tgamma will have signalled the error already).
   return floor(result + 0.5);
}

template<>
inline float factorial<float>(unsigned i)
{
   if(i <= max_factorial<float>::value)
      return unchecked_factorial<float>(i);
   return tools::overflow_error<float>(BOOST_CURRENT_FUNCTION, 0);
}

template<>
inline double factorial<double>(unsigned i)
{
   if(i <= max_factorial<double>::value)
      return unchecked_factorial<double>(i);
   return tools::overflow_error<double>(BOOST_CURRENT_FUNCTION, 0);
}

template <class T>
T double_factorial(unsigned i)
{
   using namespace std;  // ADL lookup of std names
   if(i & 1)
   {
      // odd i:
      if(i < max_factorial<T>::value)
      {
         unsigned n = (i - 1) / 2;
         return ceil(unchecked_factorial<T>(i) / (ldexp(T(1), (int)n) * unchecked_factorial<T>(n)) - 0.5f);
      }
      //
      // Fallthrough: i is too large to use table lookup, try the 
      // gamma function instead.
      //
      T result = boost::math::tgamma(static_cast<T>(i) / 2 + 1) / sqrt(constants::pi<T>());
      if(ldexp(tools::max_value<T>(), -static_cast<int>(i+1) / 2) > result)
         return ceil(result * ldexp(T(1), (i+1) / 2) - 0.5f);
   }
   else
   {
      // even i:
      unsigned n = i / 2;
      T result = factorial<T>(n);
      if(ldexp(tools::max_value<T>(), -(int)n) > result)
         return result * ldexp(T(1), (int)n);
   }
   //
   // If we fall through to here then the result is infinite:
   //
   return tools::overflow_error<T>(BOOST_CURRENT_FUNCTION, 0);
}

template <class T>
T falling_factorial(T x, unsigned n);

template <class T>
T rising_factorial(T x, int n)
{
   if(x < 0)
   {
      //
      // For x less than zero, we really have a falling
      // factorial, modulo a possible change of sign.
      //
      // Note that the falling factorial isn't defined
      // for negative n, so we'll get rid of that case
      // first:
      //
      bool inv = false;
      if(n < 0)
      {
         x += n;
         n = -n;
         inv = true;
      }
      T result = ((n&1) ? -1 : 1) * falling_factorial(-x, n);
      if(inv)
         result = 1 / result;
      return result;
   }
   if(n == 0)
      return 1;
   //
   // We don't optimise this for small n, because
   // tgamma_delta_ratio is alreay optimised for that
   // use case:
   //
   return 1 / tgamma_delta_ratio(x, static_cast<T>(n));
}

template <class T>
inline T falling_factorial(T x, unsigned n)
{
   if(x == 0)
      return 0;
   if(x < 0)
   {
      //
      // For x < 0 we really have a rising factorial
      // modulo a possible change of sign:
      //
      return (n&1 ? -1 : 1) * rising_factorial(-x, n);
   }
   if(n == 0)
      return 1;
   if(x < n-1)
   {
      //
      // x+1-n will be negative and tgamma_delta_ratio won't
      // handle it, split the product up into three parts:
      //
      T xp1 = x + 1;
      unsigned n2 = tools::real_cast<unsigned>(floor(xp1));
      if(n2 == xp1)
         return 0;
      T result = tgamma_delta_ratio(xp1, -static_cast<T>(n2));
      x -= n2;
      result *= x;
      ++n2;
      if(n2 < n)
         result *= falling_factorial(x - 1, n - n2);
      return result;
   }
   //
   // Simple case: just the ratio of two
   // (positive argument) gamma functions.
   // Note that we don't optimise this for small n, 
   // because tgamma_delta_ratio is alreay optimised
   // for that use case:
   //
   return tgamma_delta_ratio(x + 1, -static_cast<T>(n));
}

} // namespace math
} // namespace boost

#endif // BOOST_MATH_SP_FACTORIALS_HPP
