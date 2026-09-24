// Copyright John Maddock 2008.
// Copyright Matt Borland 2024

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

   BOOST_MATH_EXPORT template <class T, class Policy>
   BOOST_MATH_GPU_ENABLED typename tools::promote_args<T>::type trunc(const T& v, const Policy& pol);
   BOOST_MATH_EXPORT template <class T>
   BOOST_MATH_GPU_ENABLED typename tools::promote_args<T>::type trunc(const T& v);
   BOOST_MATH_EXPORT template <class T, class Policy>
   BOOST_MATH_GPU_ENABLED int itrunc(const T& v, const Policy& pol);
   BOOST_MATH_EXPORT template <class T>
   BOOST_MATH_GPU_ENABLED int itrunc(const T& v);
   BOOST_MATH_EXPORT template <class T, class Policy>
   BOOST_MATH_GPU_ENABLED long ltrunc(const T& v, const Policy& pol);
   BOOST_MATH_EXPORT template <class T>
   BOOST_MATH_GPU_ENABLED long ltrunc(const T& v);
   BOOST_MATH_EXPORT template <class T, class Policy>
   BOOST_MATH_GPU_ENABLED long long lltrunc(const T& v, const Policy& pol);
   BOOST_MATH_EXPORT template <class T>
   BOOST_MATH_GPU_ENABLED long long lltrunc(const T& v);
   BOOST_MATH_EXPORT template <class T, class Policy>
   BOOST_MATH_GPU_ENABLED typename tools::promote_args<T>::type round(const T& v, const Policy& pol);
   BOOST_MATH_EXPORT template <class T>
   BOOST_MATH_GPU_ENABLED typename tools::promote_args<T>::type round(const T& v);
   BOOST_MATH_EXPORT template <class T, class Policy>
   BOOST_MATH_GPU_ENABLED int iround(const T& v, const Policy& pol);
   BOOST_MATH_EXPORT template <class T>
   BOOST_MATH_GPU_ENABLED int iround(const T& v);
   BOOST_MATH_EXPORT template <class T, class Policy>
   BOOST_MATH_GPU_ENABLED long lround(const T& v, const Policy& pol);
   BOOST_MATH_EXPORT template <class T>
   BOOST_MATH_GPU_ENABLED long lround(const T& v);
   BOOST_MATH_EXPORT template <class T, class Policy>
   BOOST_MATH_GPU_ENABLED long long llround(const T& v, const Policy& pol);
   BOOST_MATH_EXPORT template <class T>
   BOOST_MATH_GPU_ENABLED long long llround(const T& v);
   BOOST_MATH_EXPORT template <class T, class Policy>
   BOOST_MATH_GPU_ENABLED T modf(const T& v, T* ipart, const Policy& pol);
   BOOST_MATH_EXPORT template <class T>
   BOOST_MATH_GPU_ENABLED T modf(const T& v, T* ipart);
   BOOST_MATH_EXPORT template <class T, class Policy>
   BOOST_MATH_GPU_ENABLED T modf(const T& v, int* ipart, const Policy& pol);
   BOOST_MATH_EXPORT template <class T>
   BOOST_MATH_GPU_ENABLED T modf(const T& v, int* ipart);
   BOOST_MATH_EXPORT template <class T, class Policy>
   BOOST_MATH_GPU_ENABLED T modf(const T& v, long* ipart, const Policy& pol);
   BOOST_MATH_EXPORT template <class T>
   BOOST_MATH_GPU_ENABLED T modf(const T& v, long* ipart);
   BOOST_MATH_EXPORT template <class T, class Policy>
   BOOST_MATH_GPU_ENABLED T modf(const T& v, long long* ipart, const Policy& pol);
   BOOST_MATH_EXPORT template <class T>
   BOOST_MATH_GPU_ENABLED T modf(const T& v, long long* ipart);
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

