///////////////////////////////////////////////////////////////////////////////
//  Copyright 2014 Anton Bikineev
//  Copyright 2014 Christopher Kormanyos
//  Copyright 2014 John Maddock
//  Copyright 2014 Paul Bristow
//  Distributed under the Boost
//  Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_HYPERGEOMETRIC_1F1_HPP
#define BOOST_MATH_HYPERGEOMETRIC_1F1_HPP

#include <boost/math/policies/policy.hpp>
#include <boost/math/policies/error_handling.hpp>
#include <boost/math/special_functions/detail/hypergeometric_series.hpp>
#include <boost/math/special_functions/detail/hypergeometric_asym.hpp>
#include <boost/math/special_functions/detail/hypergeometric_rational.hpp>
#include <boost/math/special_functions/detail/hypergeometric_1F1_recurrence.hpp>
#include <boost/math/special_functions/detail/hypergeometric_pade.hpp>
#include <boost/math/special_functions/detail/hypergeometric_1F1_bessel.hpp>
#include <boost/math/special_functions/detail/hypergeometric_1F1_scaled_series.hpp>
#include <boost/math/special_functions/detail/hypergeometric_pFq_checked_series.hpp>
#include <boost/math/special_functions/detail/hypergeometric_1f1_addition_theorems_on_z.hpp>

namespace boost { namespace math { namespace detail {

   // check when 1F1 series can't decay to polynom
   template <class T>
   inline bool check_hypergeometric_1F1_parameters(const T& a, const T& b)
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
   T hypergeometric_1F1_imp(const T& a, const T& b, const T& z, const Policy& pol, int& log_scaling)
   {
      BOOST_MATH_STD_USING // exp, fabs, sqrt

         static const char* const function = "boost::math::hypergeometric_1F1<%1%,%1%,%1%>(%1%,%1%,%1%)";

      if ((z == 0) || (a == 0))
         return T(1);

      // undefined result:
      if (!detail::check_hypergeometric_1F1_parameters(a, b))
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
            return detail::hypergeometric_1F1_generic_series(a, b, z, pol, log_scaling, function);

         return (1 + (z / b)) * exp(z);
      }

      using std::expm1;
      if ((a == 1) && (b == 2))
         return expm1(z) / z;

      if (a == b)
         return exp(z);
      if ((b - a == b) && (fabs(z / b) < policies::get_epsilon<T, Policy>()))
         return 1;
      //
      // Asymptotic expansion for large z
      // TODO: check region for higher precision types.
      // Use recurrence relations to move to this region when a and b are also large.
      //
      if (detail::hypergeometric_1F1_asym_region(a, b, z, pol))
      {
         int saved_scale = log_scaling;
         T r = hypergeometric_1F1_asym_large_z_series(a, b, z, pol, log_scaling);
         //
         // Very occationally our convergence criteria don't quite go to full precision
         // and we end up with infinity:
         //
         if (boost::math::isfinite(r))
            return r;
         log_scaling = saved_scale;
      }

      if ((fabs(a * z / b) < 3.5) && (fabs(z * 100) < fabs(b)))
         return detail::hypergeometric_1F1_rational(a, b, z, pol);

      if (z < -1)
      {
         if (a == 1)
            return detail::hypergeometric_1F1_pade(b, z, pol);
         if ((fabs(z * a / b) < 2) && (fabs(z * (a + 10) / ((b + 10) * 3628800)) < 1))  // TODO check crossover for most accurate location
         {
            if ((boost::math::sign(b - a) == boost::math::sign(b)) && ((b > 0) || (b < -200)))
            {
               // Series is close enough to convergent that we should be OK,
               // In this domain b - a ~ b and since 1F1[a, a, z] = e^z 1F1[b-a, b, -z]
               // and 1F1[a, a, -z] = e^-z the result must necessarily be somewhere near unity.
               // We have to rule out b small and negative becuase if b crosses the origin early
               // in the series (before we're pretty much converged) then all bets are off.
               // Note that this can go badly wrong when b and z are both large and negative,
               // in that situation the series goes in waves of large and small values which
               // may or may not cancel out.  Likewise the initial part of the series may or may
               // not converge, and even if it does may or may not give a correct answer!
               // For example 1F1[-small, -1252.5, -1043.7] can loose up to ~800 digits due to
               // cancellation and is basically incalculable via this method.
               return hypergeometric_1F1_checked_series_impl(a, b, z, pol, log_scaling);
            }
         }
         if ((b < 4 * a) && (a < 0) && (-a < policies::get_max_series_iterations<Policy>()))  // TODO check crosover for best location
         {
            // Without this we get into an area where the series doesn't converge if b - a ~ b
            return hypergeometric_1F1_backward_recurrence_for_negative_a(a, b, z, pol, function);
         }
         if ((a > 0) && (b + 1 < a))
         {
            // Moving to a larger b value will allow us to apply Kummer's relation below:
            //return hypergeometric_1F1_fwd_on_b_imp(a, b, z, pol, log_scaling);
         }

         // Let's otherwise make z positive (almost always)
         // by Kummer's transformation
         // (we also don't transform if z belongs to [-1,0])
         int scaling = itrunc(z);
         T r = exp(z - scaling) * detail::hypergeometric_1F1_imp<T>(b_minus_a, b, -z, pol, log_scaling);
         log_scaling += scaling;
         return r;
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
      if (!series_is_divergent && (a < 0) && (b < 0) && (b > a))
      {
         //
         // We need to exclude situations where we're over the initial "hump"
         // in the series terms (ie series has already converged by the time
         // b crosses the origin:
         //
         //T fa = fabs(a);
         //T fb = fabs(b);
         T convergence_point = sqrt((a - 1) * (a - b)) - a;
         if (-b < convergence_point)
         {
            T n = -floor(b);
            series_is_divergent = (a + n) * z / ((b + n) * n) < -1;
         }
      }
      if (series_is_divergent && (b < -1) && (b > -5) && (a > b))
         series_is_divergent = false;  // don't bother with divergence, series will be OK

      //
      // Test for alternating series due to negative a,
      // in particular, see if the series is initially divergent
      // If so use the recurrence relation on a:
      //
      if (series_is_divergent)
      {
         if((a < 0) && (floor(a) == a))
            // This works amazingly well for negative integer a:
            return hypergeometric_1F1_backward_recurrence_for_negative_a(a, b, z, pol, function);
         //
         // In what follows we have to set limits on how large z can be otherwise
         // the Bessel series become large and divergent and all the digits cancel out.
         // The criteria are distinctly empiracle rather than based on a firm analysis
         // of the terms in the series.  In some cases we can use recurrence on z to
         // stretch out to larger z values, but these series tend to be subject to wild 
         // cancellation too.
         //
         if (b > 0)
         {
            T z_limit = fabs((2 * a - b) / (sqrt(fabs(a))));
            if (z < z_limit)
               return detail::hypergeometric_1F1_AS_13_3_7_tricomi(a, b, z, pol, log_scaling);

            int k = 1 + boost::math::itrunc(z - z_limit);
            // If k is too large we destroy all the digits in the result:
            if (k < 50)
            {
               return boost::math::detail::hypergeometric_1f1_recurrence_on_z_minus_zero(a, b, z - k, k, pol);
            }
         }
         else  // b < 0
         {
            if (a < 0)
            {
               T z_limit = fabs((2 * a - b) / (sqrt(fabs(a))));
               if (z < z_limit)
                  return detail::hypergeometric_1F1_AS_13_3_7_tricomi(a, b, z, pol, log_scaling);
            }
            else if(z < fabs((2 * a - b) / (sqrt(fabs(a * b)))))
               return detail::hypergeometric_1F1_AS_13_3_7_tricomi(a, b, z, pol, log_scaling);
         }

         // If we get here, then we've run out of methods to try, use the checked series which will
         // raise an error if the result is garbage:
         return hypergeometric_1F1_checked_series_impl(a, b, z, pol, log_scaling);
      }
      
      return detail::hypergeometric_1F1_generic_series(a, b, z, pol, log_scaling, function);
   }

   template <class T, class Policy>
   inline T hypergeometric_1F1_imp(const T& a, const T& b, const T& z, const Policy& pol)
   {
      BOOST_MATH_STD_USING // exp, fabs, sqrt
      int log_scaling = 0;
      T result = hypergeometric_1F1_imp(a, b, z, pol, log_scaling);
      //
      // Actual result will be result * e^log_scaling.
      //
      int max_scaling = boost::math::itrunc(boost::math::tools::log_max_value<T>()) - 2;
      T max_scale_factor = exp(T(max_scaling));

      while (log_scaling > max_scaling)
      {
         result *= max_scale_factor;
         log_scaling -= max_scaling;
      }
      while (log_scaling < -max_scaling)
      {
         result /= max_scale_factor;
         log_scaling += max_scaling;
      }
      if (log_scaling)
         result *= exp(T(log_scaling));
      return result;
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
      detail::hypergeometric_1F1_imp<value_type>(
         static_cast<value_type>(a),
         static_cast<value_type>(b),
         static_cast<value_type>(z),
         forwarding_policy()),
      "boost::math::hypergeometric_1F1<%1%>(%1%,%1%,%1%)");
}

template <class T1, class T2, class T3>
inline typename tools::promote_args<T1, T2, T3>::type hypergeometric_1F1(T1 a, T2 b, T3 z)
{
   return hypergeometric_1F1(a, b, z, policies::policy<>());
}


  } } // namespace boost::math

#endif // BOOST_MATH_HYPERGEOMETRIC_HPP
