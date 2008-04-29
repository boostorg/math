//  (C) Copyright John Maddock 2008.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_SPECIAL_NEXT_HPP
#define BOOST_MATH_SPECIAL_NEXT_HPP

#ifdef _MSC_VER
#pragma once
#endif

#include <boost/math/policies/error_handling.hpp>
#include <boost/math/special_functions/fpclassify.hpp>
#include <boost/math/special_functions/sign.hpp>

#ifdef BOOST_MSVC
#include <float.h>
#endif

namespace boost{ namespace math{

namespace detail{

template <class T>
inline T get_smallest_value(mpl::true_ const&)
{
   return std::numeric_limits<T>::denorm_min();
}

template <class T>
inline T get_smallest_value(mpl::false_ const&)
{
   return tools::min_value<T>();
}

template <class T>
inline T get_smallest_value()
{
   return get_smallest_value<T>(mpl::bool_<std::numeric_limits<T>::is_specialized && std::numeric_limits<T>::has_denorm>());
}

}

template <class T, class Policy>
T next_greater(const T& val, const Policy& pol)
{
   BOOST_MATH_STD_USING
   int expon;
   static const char* function = "next_greater<%1%>(%1%)";

   if(!(boost::math::isfinite)(val))
      return policies::raise_domain_error<T>(
         function,
         "Argument must be finite, but got %1%", val, pol);

   if(val >= tools::max_value<T>())
      return policies::raise_overflow_error<T>(function, 0, pol);

   if(val == 0)
      return detail::get_smallest_value<T>();

   if(-0.5f == frexp(val, &expon))
      --expon; // reduce exponent when val is a power of two, and negative.
   T diff = ldexp(T(1), expon - tools::digits<T>());
   if(diff == 0)
      diff = detail::get_smallest_value<T>();
   return val + diff;
}

#ifdef BOOST_MSVC
template <class Policy>
inline double next_greater(const double& val, const Policy& pol)
{
   static const char* function = "next_greater<%1%>(%1%)";

   if(!(boost::math::isfinite)(val))
      return policies::raise_domain_error<double>(
         function,
         "Argument must be finite, but got %1%", val, pol);

   if(val >= tools::max_value<double>())
      return policies::raise_overflow_error<double>(function, 0, pol);

   return ::_nextafter(val, tools::max_value<double>());
}
#endif

template <class T>
inline T next_greater(const T& val)
{
   return next_greater(val, policies::policy<>());
}

template <class T, class Policy>
T next_less(const T& val, const Policy& pol)
{
   BOOST_MATH_STD_USING
   int expon;
   static const char* function = "next_less<%1%>(%1%)";

   if(!(boost::math::isfinite)(val))
      return policies::raise_domain_error<T>(
         function,
         "Argument must be finite, but got %1%", val, pol);

   if(val <= -tools::max_value<T>())
      return -policies::raise_overflow_error<T>(function, 0, pol);

   if(val == 0)
      return -detail::get_smallest_value<T>();

   T remain = frexp(val, &expon);
   if(remain == 0.5)
      --expon; // when val is a power of two we must reduce the exponent
   T diff = ldexp(T(1), expon - tools::digits<T>());
   if(diff == 0)
      diff = detail::get_smallest_value<T>();
   return val - diff;
}

#ifdef BOOST_MSVC
template <class Policy>
inline double next_less(const double& val, const Policy& pol)
{
   static const char* function = "next_less<%1%>(%1%)";

   if(!(boost::math::isfinite)(val))
      return policies::raise_domain_error<double>(
         function,
         "Argument must be finite, but got %1%", val, pol);

   if(val <= -tools::max_value<double>())
      return -policies::raise_overflow_error<double>(function, 0, pol);

   return ::_nextafter(val, -tools::max_value<double>());
}
#endif

template <class T>
inline T next_less(const T& val)
{
   return next_less(val, policies::policy<>());
}

template <class T, class Policy>
inline T nextafter(const T& val, const T& direction, const Policy& pol)
{
   return val < direction ? boost::math::next_greater(val, pol) : val == direction ? val : boost::math::next_less(val, pol);
}

template <class T>
inline T nextafter(const T& val, const T& direction)
{
   return nextafter(val, direction, policies::policy<>());
}

template <class T, class Policy>
T edit_distance(const T& a, const T& b, const Policy& pol)
{
   BOOST_MATH_STD_USING
   //
   // Error handling:
   //
   static const char* function = "edit_distance<%1%>(%1%, %1%)";
   if(!(boost::math::isfinite)(a))
      return policies::raise_domain_error<T>(
         function,
         "Argument a must be finite, but got %1%", a, pol);
   if(!(boost::math::isfinite)(b))
      return policies::raise_domain_error<T>(
         function,
         "Argument b must be finite, but got %1%", b, pol);
   //
   // Special cases:
   //
   if(a == b)
      return 0;
   if(a == 0)
      return 1 + edit_distance(boost::math::sign(b) * detail::get_smallest_value<T>(), b, pol);
   if(b == 0)
      return 1 + edit_distance(boost::math::sign(a) * detail::get_smallest_value<T>(), a, pol);
   if(boost::math::sign(a) != boost::math::sign(b))
      return 2 + edit_distance(boost::math::sign(b) * detail::get_smallest_value<T>(), b, pol)
         + edit_distance(boost::math::sign(a) * detail::get_smallest_value<T>(), a, pol);

   if((std::min)(fabs(a), fabs(b)) / (std::max)(fabs(a), fabs(b)) < 2 * tools::epsilon<T>())
   {
      bool biga = fabs(a) > fabs(b);
      T split = ldexp(biga ? b : a, tools::digits<T>() - 2);
      return edit_distance(a, split, pol) + edit_distance(split, b, pol);
   }

   BOOST_MATH_STD_USING
   int expon;
   //
   // We're going to left shift the result by the exponent of the 
   // smaller of the two values (irrespective of sign):
   //
   T mv = (std::min)(fabs(a), fabs(b));
   //
   // Note that if mv is a denorm then the usual formula fails
   // because we actually have fewer than tools::digits<T>()
   // significant bits in the representation:
   //
   frexp((boost::math::fpclassify(mv) == FP_SUBNORMAL) ? tools::min_value<T>() : mv, &expon);
   expon = tools::digits<T>() - expon;
   //
   // Use compensated double-double addition to avoid rounding 
   // errors in the subtraction, note this will still fail if
   // the two values differ by many orders of magnitute:
   //
   T mb = -b;
   T x = a + mb;
   T z = x - a;
   T y = (a - (x - z)) + (mb - z);
   if(x < 0)
   {
      x = -x;
      y = -y;
   }
   
   T result = ldexp(x, expon) + ldexp(y, expon);
   //
   // Result must be an integer:
   //
   BOOST_ASSERT(result == floor(result));
   return result;
}

template <class T>
T edit_distance(const T& a, const T& b)
{
   return boost::math::edit_distance(a, b, policies::policy<>());
}

}} // namespaces

#endif // BOOST_MATH_SPECIAL_NEXT_HPP

