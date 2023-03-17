//  Copyright John Maddock 2007.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_ROUND_HPP
#define BOOST_MATH_ROUND_HPP

#ifdef _MSC_VER
#pragma once
#endif

#include <boost/math/tools/config.hpp>
#include <boost/math/policies/error_handling.hpp>
#include <boost/math/special_functions/math_fwd.hpp>
#include <boost/math/special_functions/fpclassify.hpp>
#include <type_traits>
#include <limits>
#include <cmath>

namespace boost{ namespace math{

namespace detail{

// https://stackoverflow.com/questions/8905246/how-to-check-if-float-can-be-exactly-represented-as-an-integer/17822304#17822304
template <typename T, typename ResultType>
inline ResultType float_to_int(T x)
{
   BOOST_MATH_STD_USING

   constexpr int sign = std::is_signed<ResultType>::value ? -1 : 0;

   T y = floor(x);
   if (y < static_cast<T>(0.0L))
   {
      return 0;
   }
   
   if (y >= ldexp(static_cast<T>(1.0L), (sizeof(ResultType) * CHAR_BIT) + sign) + sign)
   {
      return (std::numeric_limits<ResultType>::max)();
   }

   return static_cast<ResultType>(y);
}

template <typename T, typename TargetType>
inline bool is_representable(T x)
{
   BOOST_MATH_STD_USING
   return (floor(x) == x && x >= static_cast<T>(0.0L) && x < ldexp(static_cast<T>(1.0L), sizeof(TargetType) * CHAR_BIT));
}

template <class T, class Policy>
inline tools::promote_args_t<T> round(const T& v, const Policy& pol, const std::false_type&)
{
   BOOST_MATH_STD_USING
   using result_type = tools::promote_args_t<T>;

   if(!(boost::math::isfinite)(v))
   {
      return policies::raise_rounding_error("boost::math::round<%1%>(%1%)", nullptr, static_cast<result_type>(v), static_cast<result_type>(v), pol);
   }
   //
   // The logic here is rather convoluted, but avoids a number of traps,
   // see discussion here https://github.com/boostorg/math/pull/8
   //
   if (-0.5 < v && v < 0.5)
   {
      // special case to avoid rounding error on the direct
      // predecessor of +0.5 resp. the direct successor of -0.5 in
      // IEEE floating point types
      return static_cast<result_type>(0);
   }
   else if (v > 0)
   {
      // subtract v from ceil(v) first in order to avoid rounding
      // errors on largest representable integer numbers
      result_type c(ceil(v));
      return 0.5 < c - v ? c - 1 : c;
   }
   else
   {
      // see former branch
      result_type f(floor(v));
      return 0.5 < v - f ? f + 1 : f;
   }
}
template <class T, class Policy>
inline tools::promote_args_t<T> round(const T& v, const Policy&, const std::true_type&)
{
   return v;
}

} // namespace detail

template <class T, class Policy>
inline tools::promote_args_t<T> round(const T& v, const Policy& pol)
{
   return detail::round(v, pol, std::integral_constant<bool, detail::is_integer_for_rounding<T>::value>());
}
template <class T>
inline tools::promote_args_t<T> round(const T& v)
{
   return round(v, policies::policy<>());
}
//
// The following functions will not compile unless T has an
// implicit conversion to the integer types.  For user-defined
// number types this will likely not be the case.  In that case
// these functions should either be specialized for the UDT in
// question, or else overloads should be placed in the same
// namespace as the UDT: these will then be found via argument
// dependent lookup.  See our concept archetypes for examples.
//
// Non-standard numeric limits syntax "(std::numeric_limits<int>::max)()"
// is to avoid macro substiution from MSVC
// https://stackoverflow.com/questions/27442885/syntax-error-with-stdnumeric-limitsmax
//
template <class T, class Policy>
inline int iround(const T& v, const Policy& pol)
{
   BOOST_MATH_STD_USING
   using result_type = tools::promote_args_t<T>;

   T r = boost::math::round(v, pol);
   if(r > static_cast<result_type>((std::numeric_limits<int>::max)()) || r < static_cast<result_type>((std::numeric_limits<int>::min)()))
   {
      return static_cast<int>(policies::raise_rounding_error("boost::math::iround<%1%>(%1%)", nullptr, v, 0, pol));
   }
   return static_cast<int>(r);
}
template <class T>
inline int iround(const T& v)
{
   return iround(v, policies::policy<>());
}

template <class T, class Policy>
inline long lround(const T& v, const Policy& pol)
{
   BOOST_MATH_STD_USING
   using result_type = tools::promote_args_t<T>;
   T r = boost::math::round(v, pol);
   if(r > static_cast<result_type>((std::numeric_limits<long>::max)()) || r < static_cast<result_type>((std::numeric_limits<long>::min)()))
   {
      return static_cast<long int>(policies::raise_rounding_error("boost::math::lround<%1%>(%1%)", nullptr, v, 0L, pol));
   }
   return static_cast<long int>(r);
}
template <class T>
inline long lround(const T& v)
{
   return lround(v, policies::policy<>());
}

template <class T, class Policy>
inline long long llround(const T& v, const Policy& pol)
{
   BOOST_MATH_STD_USING
   using result_type = tools::promote_args_t<T>;

   T r = boost::math::round(v, pol);
   long long return_val = boost::math::detail::float_to_int<result_type, long long>(r);
   bool representable = boost::math::detail::is_representable<result_type, long long>(r);

   if ((return_val == (std::numeric_limits<long long>::max)() && !representable) ||
       r < static_cast<result_type>((std::numeric_limits<long long>::min)()) ||
       r > static_cast<result_type>((std::numeric_limits<long long>::max)()))
   {
      return static_cast<long long>(policies::raise_rounding_error("boost::math::llround<%1%>(%1%)", nullptr, v, static_cast<long long>(0), pol));
   }

   if (r < 0)
   {
      return_val = static_cast<long long>(r);
   }

   return return_val;
}
template <class T>
inline long long llround(const T& v)
{
   return llround(v, policies::policy<>());
}

}} // namespaces

#endif // BOOST_MATH_ROUND_HPP
