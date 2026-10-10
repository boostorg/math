///////////////////////////////////////////////////////////////////////////////
//  Copyright 2014 Anton Bikineev
//  Copyright 2014 Christopher Kormanyos
//  Copyright 2014 John Maddock
//  Copyright 2014 Paul Bristow
//  Distributed under the Boost
//  Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_HYPERGEOMETRIC_2F0_HPP
#define BOOST_MATH_HYPERGEOMETRIC_2F0_HPP

#include <boost/math/policies/policy.hpp>
#include <boost/math/policies/error_handling.hpp>
#include <boost/math/special_functions/detail/hypergeometric_series.hpp>
#include <boost/math/special_functions/laguerre.hpp>
#include <boost/math/special_functions/hermite.hpp>
#include <boost/math/special_functions/fpclassify.hpp>
#include <boost/math/tools/fraction.hpp>

#ifndef BOOST_MATH_BUILD_MODULE
#include <cmath>
#endif

BOOST_MATH_NAMESPACE_BEGIN namespace detail {

   template <class T>
   struct hypergeometric_2F0_cf
   {
      //
      // We start this continued fraction at b on index -1
      // and treat the -1 and 0 cases as special cases.
      // We do this to avoid adding the continued fraction result
      // to 1 so that we can accurately evaluate for small results
      // as well as large ones.  See  http://functions.wolfram.com/07.31.10.0002.01
      //
      T a1, a2, z;
      int k;
      hypergeometric_2F0_cf(T a1_, T a2_, T z_) : a1(a1_), a2(a2_), z(z_), k(-2) {}
      typedef std::pair<T, T> result_type;

      result_type operator()()
      {
         ++k;
         if (k <= 0)
            return std::make_pair(z * a1 * a2, 1);
         return std::make_pair(-z * (a1 + k) * (a2 + k) / (k + 1), 1 + z * (a1 + k) * (a2 + k) / (k + 1));
      }
   };

   template <class T, class Policy>
   T hypergeometric_2F0_cf_imp(T a1, T a2, T z, const Policy& pol, const char* function)
   {
      using namespace BOOST_MATH_NAMESPACE;
      hypergeometric_2F0_cf<T> evaluator(a1, a2, z);
      std::uintmax_t max_iter = policies::get_max_series_iterations<Policy>();
      T cf = tools::continued_fraction_b(evaluator, policies::get_epsilon<T, Policy>(), max_iter);
      policies::check_series_iterations<T>(function, max_iter, pol);
      return cf;
   }


   //
   // 2F0(-n, b; z) for b > 0 and z > 0, where the series alternates and cancels badly.  Instead use the recurrence
   // g[k+1] = (1 - (k + b) z) g[k] + k z g[k-1], with g[0] = 1 and g[1] = 1 - b z.  Where g is close to the
   // minimal solution of the recurrence its rounding errors grow, so as in laguerre.hpp we find the rounding error
   // of each step almost exactly, and carry it forward with the same recurrence.  The result is then correctly
   // rounded unless the uncorrected recurrence was out by more than a few percent, which only happens for n > 50
   // and n z between about 1 and 20, and there we raise an evaluation_error rather than return a wrong value.
   // The rounded products come from fma(a, b, 0), which the compiler can't contract into a later addition:
   //
   template <class T, class Policy>
   T hypergeometric_2F0_neg_int_recurrence(unsigned n, const T& b, const T& z, const Policy& pol, const char* function)
   {
      BOOST_MATH_STD_USING
      using std::fma;
      T t;
      // g1 + e1 = 1 - b z exactly:
      T bz = fma(b, z, T(0));
      T rbz = fma(b, z, -bz);
      T g1 = 1 - bz;
      t = g1 - 1;
      T e1 = ((1 - (g1 - t)) + (-bz - t)) - rbz;
      T g0 = 1;
      T e0 = 0;
      for (unsigned k = 1; k < n; ++k)
      {
         // a + ea = 1 - (k + b) z and c + ec = k z exactly:
         T u = k + b;
         t = u - k;
         T eu = (k - (u - t)) + (b - t);
         T v = fma(u, z, T(0));
         T rv = fma(u, z, -v) + eu * z;
         T a = 1 - v;
         t = a - 1;
         T ea = ((1 - (a - t)) + (-v - t)) - rv;
         T c = fma(T(k), z, T(0));
         T ec = fma(T(k), z, -c);
         // s1 + r1 = a g1, s2 + r2 = c g0 and d + ed = s1 + s2 exactly:
         T s1 = fma(a, g1, T(0));
         T r1 = fma(a, g1, -s1);
         T s2 = fma(c, g0, T(0));
         T r2 = fma(c, g0, -s2);
         T d = s1 + s2;
         t = d - s1;
         T ed = (s1 - (d - t)) + (s2 - t);
         T e2 = a * e1 + c * e0 + ed + r1 + r2 + ea * g1 + ec * g0;
         g0 = g1;
         g1 = d;
         e0 = e1;
         e1 = e2;
      }
      if (!(BOOST_MATH_NAMESPACE::isfinite)(g1))
         return BOOST_MATH_NAMESPACE::policies::raise_overflow_error<T>(function, nullptr, pol);
      if (!(fabs(e1) < fabs(g1) / 32))
         return BOOST_MATH_NAMESPACE::policies::raise_evaluation_error<T>(function, "The recurrence loses too many digits to give an accurate result, last value was %1%", T(g1 + e1), pol);
      return g1 + e1;
   }

   template <class T, class Policy>
   inline T hypergeometric_2F0_imp(T a1, T a2, const T& z, const Policy& pol, bool asymptotic = false)
   {
      //
      // The terms in this series go to infinity unless one of a1 and a2 is a negative integer.
      //
      using std::swap;
      BOOST_MATH_STD_USING

      static const char* const function = "boost::math::hypergeometric_2F0<%1%,%1%,%1%>(%1%,%1%,%1%)";

      if (z == 0)
         return 1;

      bool is_a1_integer = (a1 == floor(a1));
      bool is_a2_integer = (a2 == floor(a2));

      if (!asymptotic && !is_a1_integer && !is_a2_integer)
         return BOOST_MATH_NAMESPACE::policies::raise_overflow_error<T>(function, nullptr, pol);
      if (!is_a1_integer || (a1 > 0))
      {
         swap(a1, a2);
         swap(is_a1_integer, is_a2_integer);
      }
      //
      // At this point a1 must be a negative integer:
      //
      if(!asymptotic && (!is_a1_integer || (a1 > 0)))
         return BOOST_MATH_NAMESPACE::policies::raise_overflow_error<T>(function, nullptr, pol);
      //
      // Special cases first:
      //
      if (a1 == 0)
         return 1;
      if ((a1 == a2 - 0.5f) && (z < 0))
      {
         // http://functions.wolfram.com/07.31.03.0083.01
         int n = static_cast<int>(static_cast<std::uintmax_t>(BOOST_MATH_NAMESPACE::lltrunc(-2 * a1)));
         T smz = sqrt(-z);
         return static_cast<T>(pow(2 / smz, T(-n)) * BOOST_MATH_NAMESPACE::hermite(n, 1 / smz, pol));  // Warning suppression: integer power returns at least a double
      }

      if (is_a1_integer && is_a2_integer)
      {
         if ((a1 < 1) && (a2 <= a1))
         {
            const unsigned int n = static_cast<unsigned int>(static_cast<std::uintmax_t>(BOOST_MATH_NAMESPACE::lltrunc(-a1)));
            const unsigned int m = static_cast<unsigned int>(static_cast<std::uintmax_t>(BOOST_MATH_NAMESPACE::lltrunc(-a2 - n)));

            return (pow(z, T(n)) * BOOST_MATH_NAMESPACE::factorial<T>(n, pol)) *
               BOOST_MATH_NAMESPACE::laguerre(n, m, -(1 / z), pol);
         }
         else if ((a2 < 1) && (a1 <= a2))
         {
            // function is symmetric for a1 and a2
            const unsigned int n = static_cast<unsigned int>(static_cast<std::uintmax_t>(BOOST_MATH_NAMESPACE::lltrunc(-a2)));
            const unsigned int m = static_cast<unsigned int>(static_cast<std::uintmax_t>(BOOST_MATH_NAMESPACE::lltrunc(-a1 - n)));

            return (pow(z, T(n)) * BOOST_MATH_NAMESPACE::factorial<T>(n, pol)) *
               BOOST_MATH_NAMESPACE::laguerre(n, m, -(1 / z), pol);
         }
      }

      if ((a2 > 0) && (z > 0))
         return hypergeometric_2F0_neg_int_recurrence(static_cast<unsigned>(BOOST_MATH_NAMESPACE::lltrunc(-a1)), a2, z, pol, function);

      if ((a1 * a2 * z < 0) && (a2 < -5) && (fabs(a1 * a2 * z) > 0.5))
      {
         // Series is alternating and maybe divergent at least for the first few terms
         // (until a2 goes positive), try the continued fraction:
         return hypergeometric_2F0_cf_imp(a1, a2, z, pol, function);
      }

      return detail::hypergeometric_2F0_generic_series(a1, a2, z, pol);
   }

} // namespace detail

BOOST_MATH_EXPORT template <class T1, class T2, class T3, class Policy>
inline typename tools::promote_args<T1, T2, T3>::type hypergeometric_2F0(T1 a1, T2 a2, T3 z, const Policy& /* pol */)
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
      detail::hypergeometric_2F0_imp<value_type>(
         static_cast<value_type>(a1),
         static_cast<value_type>(a2),
         static_cast<value_type>(z),
         forwarding_policy()),
      "boost::math::hypergeometric_2F0<%1%>(%1%,%1%,%1%)");
}

BOOST_MATH_EXPORT template <class T1, class T2, class T3>
inline typename tools::promote_args<T1, T2, T3>::type hypergeometric_2F0(T1 a1, T2 a2, T3 z)
{
   return hypergeometric_2F0(a1, a2, z, policies::policy<>());
}


  BOOST_MATH_NAMESPACE_END

#endif // BOOST_MATH_HYPERGEOMETRIC_HPP
