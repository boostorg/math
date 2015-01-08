//  Copyright (c) 2015 John Maddock
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
#ifndef BOOST_MATH_ELLINT_RG_HPP
#define BOOST_MATH_ELLINT_RG_HPP

#ifdef _MSC_VER
#pragma once
#endif

#include <boost/math/special_functions/math_fwd.hpp>
#include <boost/math/tools/config.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/math/policies/error_handling.hpp>
#include <boost/math/special_functions/ellint_rd.hpp>
#include <boost/math/special_functions/ellint_rf.hpp>

namespace boost { namespace math { namespace detail{

   template <typename T, typename Policy>
   T ellint_rg_imp(T x, T y, T z, const Policy& pol)
   {
      BOOST_MATH_STD_USING
      using namespace boost::math;

      static const char* function = "boost::math::ellint_rf<%1%>(%1%,%1%,%1%)";

      if(x < 0 || y < 0 || z < 0)
      {
         return policies::raise_domain_error<T>(function,
            "domain error, all arguments must be non-negative, "
            "only sensible result is %1%.",
            std::numeric_limits<T>::quiet_NaN(), pol);
      }
      //
      // Function is symmetric in x, y and z, but we require
      // (x - z)(y - z) >= 0 to avoid cancellation error in the result
      // which implies (for example) x >= z >= y
      //
      using std::swap;
      if(x < y)
         swap(x, y);
      if(x < z)
         swap(x, z);
      if(y > z)
         swap(y, z);
      
      BOOST_ASSERT(x >= z);
      BOOST_ASSERT(z >= y);
      //
      // Special cases from http://dlmf.nist.gov/19.20#ii
      //
      if(x == y)
      {
         if(y == z)
         {
            // x = y = z
            // This also works for x = y = z = 0 presumably.
            return sqrt(x);
         }
         else if(z == 0)
         {
            // x = y, z = 0
            return constants::pi<T>() * sqrt(x) / 4;
         }
         else
         {
            // x = y, z != 0
            swap(x, z);
            return (x == 0) ? sqrt(z) / 2 : (y * ellint_rc_imp(x, y, pol) + sqrt(x)) / 2;
         }
      }
      else if(y == z)
      {
         if(x == 0)
            return constants::pi<T>() * sqrt(y) / 4;
         else
            return (y == 0) ? sqrt(x) / 2 : (y * ellint_rc_imp(x, y, pol) + sqrt(x)) / 2;
      }

      return (z * ellint_rf_imp(x, y, z, pol)
         - (x - z) * (y - z) * ellint_rd_imp(x, y, z, pol) / 3
         + sqrt(x * y / z)) / 2;
   }

} // namespace detail

template <class T1, class T2, class T3, class Policy>
inline typename tools::promote_args<T1, T2, T3>::type 
   ellint_rg(T1 x, T2 y, T3 z, const Policy& pol)
{
   typedef typename tools::promote_args<T1, T2, T3>::type result_type;
   typedef typename policies::evaluation<result_type, Policy>::type value_type;
   return policies::checked_narrowing_cast<result_type, Policy>(
      detail::ellint_rg_imp(
         static_cast<value_type>(x),
         static_cast<value_type>(y),
         static_cast<value_type>(z), pol), "boost::math::ellint_rf<%1%>(%1%,%1%,%1%)");
}

template <class T1, class T2, class T3>
inline typename tools::promote_args<T1, T2, T3>::type 
   ellint_rg(T1 x, T2 y, T3 z)
{
   return ellint_rg(x, y, z, policies::policy<>());
}

}} // namespaces

#endif // BOOST_MATH_ELLINT_RG_HPP

