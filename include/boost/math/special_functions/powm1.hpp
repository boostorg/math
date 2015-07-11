//  (C) Copyright John Maddock 2006.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_POWM1
#define BOOST_MATH_POWM1

#ifdef _MSC_VER
#pragma once
#endif

#include <boost/math/special_functions/math_fwd.hpp>
#include <boost/math/special_functions/log1p.hpp>
#include <boost/math/special_functions/expm1.hpp>
#include <boost/assert.hpp>

namespace boost{ namespace math{ namespace detail{

template <class T, class Policy>
inline T powm1_imp(const T a, const T z, const Policy& pol)
{
   BOOST_MATH_STD_USING

   //  Refer to the series expansion at z = 0 for (a-1)^z for the logic
   // in choosing this cutoff, see http://www.wolframalpha.com/input/?i=%28x%2B1%29^y
   if(fabs(a * z) < 1)
   {
      T p = fabs(a - 1) < 0.5f ? boost::math::log1p(a-1, pol) : log(a);
      p *= z;
      return boost::math::expm1(p, pol);
   }
   return pow(a, z) - 1;
}

} // detail

template <class T1, class T2>
inline typename tools::promote_args<T1, T2>::type 
   powm1(const T1 a, const T2 z)
{
   typedef typename tools::promote_args<T1, T2>::type result_type;
   return detail::powm1_imp(static_cast<result_type>(a), static_cast<result_type>(z), policies::policy<>());
}

template <class T1, class T2, class Policy>
inline typename tools::promote_args<T1, T2>::type 
   powm1(const T1 a, const T2 z, const Policy& pol)
{
   typedef typename tools::promote_args<T1, T2>::type result_type;
   return detail::powm1_imp(static_cast<result_type>(a), static_cast<result_type>(z), pol);
}

} // namespace math
} // namespace boost

#endif // BOOST_MATH_POWM1





