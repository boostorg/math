// Copyright John Maddock 2008.

// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_SPECIAL_ROUND_FWD_HPP
#define BOOST_MATH_SPECIAL_ROUND_FWD_HPP

#include <boost/math/tools/config.hpp>
#include <boost/math/tools/promotion.hpp>

#ifdef _MSC_VER
#pragma once
#endif

namespace boost
{
   namespace math
   { 

   template <class T, class Policy>
   constexpr typename tools::promote_args<T>::type trunc(const T& v, const Policy& pol);
   template <class T>
   constexpr typename tools::promote_args<T>::type trunc(const T& v);
   template <class T, class Policy>
   constexpr int itrunc(const T& v, const Policy& pol);
   template <class T>
   constexpr int itrunc(const T& v);
   template <class T, class Policy>
   constexpr long ltrunc(const T& v, const Policy& pol);
   template <class T>
   constexpr long ltrunc(const T& v);
   template <class T, class Policy>
   constexpr long long lltrunc(const T& v, const Policy& pol);
   template <class T>
   constexpr long long lltrunc(const T& v);
   template <class T, class Policy>
   constexpr typename tools::promote_args<T>::type round(const T& v, const Policy& pol);
   template <class T>
   constexpr typename tools::promote_args<T>::type round(const T& v);
   template <class T, class Policy>
   constexpr int iround(const T& v, const Policy& pol);
   template <class T>
   constexpr int iround(const T& v);
   template <class T, class Policy>
   constexpr long lround(const T& v, const Policy& pol);
   template <class T>
   constexpr long lround(const T& v);
   template <class T, class Policy>
   constexpr long long llround(const T& v, const Policy& pol);
   template <class T>
   constexpr long long llround(const T& v);
   template <class T, class Policy>
   T modf(const T& v, T* ipart, const Policy& pol);
   template <class T>
   T modf(const T& v, T* ipart);
   template <class T, class Policy>
   T modf(const T& v, int* ipart, const Policy& pol);
   template <class T>
   T modf(const T& v, int* ipart);
   template <class T, class Policy>
   T modf(const T& v, long* ipart, const Policy& pol);
   template <class T>
   T modf(const T& v, long* ipart);
   template <class T, class Policy>
   T modf(const T& v, long long* ipart, const Policy& pol);
   template <class T>
   T modf(const T& v, long long* ipart);
   }
}

#undef BOOST_MATH_STD_USING
#define BOOST_MATH_STD_USING BOOST_MATH_STD_USING_CORE\
   using boost::math::round;\
   using boost::math::iround;\
   using boost::math::lround;\
   using boost::math::trunc;\
   using boost::math::itrunc;\
   using boost::math::ltrunc;\
   using boost::math::modf;


#endif // BOOST_MATH_SPECIAL_ROUND_FWD_HPP

