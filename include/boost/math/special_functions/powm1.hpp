//  (C) Copyright John Maddock 2006.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_POWM1
#define BOOST_MATH_POWM1

#include <boost/math/special_functions/log1p.hpp>
#include <boost/math/special_functions/expm1.hpp>
#include <boost/math/special_functions/math_fwd.hpp>
#include <boost/assert.hpp>

namespace boost{ namespace math{ namespace detail{

template <class T>
inline T powm1_imp(const T a, const T z)
{
   using namespace std;

   if((fabs(a) < 1) || (fabs(z) < 1))
   {
      T p = log(a) * z;
      if(fabs(p) < 2)
         return boost::math::expm1(p);
      // otherwise fall though:
   }
   return pow(a, z) - 1;
}

} // detail

template <class T1, class T2>
inline typename tools::promote_args<T1, T2>::type 
   powm1(const T1 a, const T2 z)
{
   typedef typename tools::promote_args<T1, T2>::type result_type;
   return detail::powm1_imp(static_cast<result_type>(a), static_cast<result_type>(z));
}

} // namespace math
} // namespace boost

#endif // BOOST_MATH_POWM1



