
///////////////////////////////////////////////////////////////////////////////
//  Copyright 2014 Anton Bikineev
//  Copyright 2014 Christopher Kormanyos
//  Copyright 2014 John Maddock
//  Copyright 2014 Paul Bristow
//  Distributed under the Boost
//  Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef _BOOST_HYPERGEOMETRIC_1F1_HPP_
  #define _BOOST_HYPERGEOMETRIC_1F1_HPP_

#include <boost/math/policies/policy.hpp>
#include <boost/math/policies/error_handling.hpp>
#include <boost/math/special_functions/detail/hypergeometric_series.hpp>
#include <boost/math/special_functions/detail/hypergeometric_asym.hpp>
#include <boost/math/special_functions/detail/hypergeometric_rational.hpp>
#include <boost/math/special_functions/detail/hypergeometric_1f1_recurrence.hpp>
#include <boost/math/special_functions/detail/hypergeometric_pade.hpp>
#include <boost/math/special_functions/detail/hypergeometric_1f1_bessel.hpp>

namespace boost { namespace math { namespace detail {

   // check when 1f1 series can't decay to polynom
   template <class T>
   inline bool check_hypergeometric_1f1_parameters(const T& a, const T& b)
   {
      BOOST_MATH_STD_USING

         if ((b <= 0) && (b == floor(b)))
         {
            if ((a >= 0) || (a < b) || (a != floor(a)))
               return false;
         }

      return true;
   }

   template <class T, class Policy>
   inline T hypergeometric_1f1_imp(const T& a, const T& b, const T& z, const Policy& pol)
   {
      BOOST_MATH_STD_USING // exp, fabs, sqrt

         static const char* const function = "boost::math::hypergeometric_1f1<%1%,%1%,%1%>(%1%,%1%,%1%)";

      if ((z == 0) || (a == 0))
         return T(1);

      // undefined result:
      if (!detail::check_hypergeometric_1f1_parameters(a, b))
         return policies::raise_domain_error<T>(
            function,
            "Function is indeterminate for negative integer b = %1%.",
            b,
            pol);

      // other checks:
      if (a == -1)
         return 1 - (z / b);

      const T b_minus_a = b - a;

      // 0f0 (exp) case;
      if (b_minus_a == 0)
         return exp(z);

      if ((b_minus_a == -1))
      {
         // for negative integer a and b is reasonable to use truncated series - polynomial
         if ((a < 0) && (a == ceil(a)))
            return detail::hypergeometric_1f1_generic_series(a, b, z, pol);

         return (1 + (z / b)) * exp(z);
      }

      if ((a == 1) && (b == 2))
         return (exp(z) - 1) / z;

      // asymptotic expansion
      // check region
      if (detail::hypergeometric_1f1_asym_region(a, b, z))
      {
         // check for poles in gamma for b
         if ((b > 0) || (b != floor(b)))
         {
            //check for poles in gamma for a
            if (((a > 0) || (a != floor(a))) && (z > 0))
               return detail::hypergeometric_1f1_asym_positive_series(a, b, z, pol);

            //check for poles in gamma for b
            if (((b_minus_a > 0) || (b_minus_a != floor(b_minus_a))) && (z < 0))
               return detail::hypergeometric_1f1_asym_negative_series(a, b, z, pol);
         }
      }

      if (fabs(b) >= fabs(100 * z)) // TODO: extend to multuiprecision
         return detail::hypergeometric_1f1_rational(a, b, z, pol);

      if (z < -1)
      {
         if (a == 1)
            return detail::hypergeometric_1f1_pade(b, z, pol);

         // Let's otherwise make z positive (almost always)
         // by Kummer's transformation
         // (we also don't transform if z belongs to [-1,0])
         return exp(z) * detail::hypergeometric_1f1_imp<T>(b_minus_a, b, -z, pol);
      }
      //
      // Check for initial divergence:
      //
      bool series_is_divergent = (a + 1) * z / (b + 1) < -1;
      //
      // If series starts off non-divergent, and becomes divergent later
      // then it's because both a and b are negative, so check for later
      // divergence as well:
      //
      if ((a < 0) && (b < 0) && (b > a))
      {
         T n = -floor(b);
         series_is_divergent = (a + n) * z / ((b + n) * n) < -1;
      }
      //
      // Test for alternating series due to negative a, 
      // in particular, see if the series is initially divergent
      // If so use the recurrence relation on a:
      //
      if (series_is_divergent)
      {
         if((a < 0) && (b > 0))
         return detail::hypergeometric_1f1_backward_recurrence_for_negative_a(a, b, z, pol);

         if (detail::hypergeometric_1f1_is_a_small_enough(a))
         {
            // TODO: this part has to be researched deeper
            const bool b_is_negative_and_greater_than_z = b < 0 ? (fabs(b) > fabs(z) ? 1 : 0) : 0;
            if ((a == ceil(a)) && !b_is_negative_and_greater_than_z)
               return detail::hypergeometric_1f1_backward_recurrence_for_negative_a(a, b, z, pol);
            else if ((2 * (z  * (b - (2 * a)))) > 0) // TODO: see when this methd is bad in opposite to usual taylor
               return detail::hypergeometric_1f1_13_3_7_series(a, b, z, pol);
            else if (b < a)
               return detail::hypergeometric_1f1_backward_recurrence_for_negative_b(a, b, z, pol);
         }
      }

      return detail::hypergeometric_1f1_generic_series(a, b, z, pol);
   }

} // namespace detail

template <class T1, class T2, class T3, class Policy>
inline typename tools::promote_args<T1, T2, T3>::type hypergeometric_1F1(T1 a, T2 b, T3 z, const Policy& /* pol */)
{
   BOOST_FPU_EXCEPTION_GUARD
      typedef typename tools::promote_args<T1, T2, T3>::type result_type;
   typedef typename policies::evaluation<result_type, Policy>::type value_type;
   typedef typename policies::normalise<
      Policy,
      policies::promote_float<false>,
      policies::promote_double<false>,
      policies::discrete_quantile<>,
      policies::assert_undefined<> >::type forwarding_policy;
   return policies::checked_narrowing_cast<result_type, Policy>(
      detail::hypergeometric_1f1_imp<value_type>(
         static_cast<value_type>(a),
         static_cast<value_type>(b),
         static_cast<value_type>(z),
         forwarding_policy()),
      "boost::math::hypergeometric_1f1<%1%>(%1%,%1%,%1%)");
}

template <class T1, class T2, class T3>
inline typename tools::promote_args<T1, T2, T3>::type hypergeometric_1F1(T1 a, T2 b, T3 z)
{
   return hypergeometric_1F1(a, b, z, policies::policy<>());
}


  } } // namespace boost::math

#endif // _BOOST_HYPERGEOMETRIC_2014_04_07_HPP_
