//  Copyright Nick Thompson, 2026.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_DISTRIBUTIONS_DETAIL_SADDLE_POINT_HPP
#define BOOST_MATH_DISTRIBUTIONS_DETAIL_SADDLE_POINT_HPP

//
// The building blocks of Loader's saddle-point evaluation of discrete densities:
// C. Loader, "Fast and Accurate Computation of Binomial Probabilities", 2000.
//

#include <boost/math/tools/config.hpp>

#ifndef BOOST_MATH_HAS_GPU_SUPPORT

#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/special_functions/log1p.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/math/policies/error_handling.hpp>

namespace boost { namespace math { namespace detail {

// ln(n!) minus Stirling's approximation (n + 1/2) ln n - n + ln(2 pi)/2, for n above the point
// where its Bernoulli series converges; that is minimum_argument_for_bernoulli_recursion, with a margin
// because just above it the series can still fail to converge (at n = 6.05 for double).
template <class RealType, class Policy>
inline bool stirlerr_series_converges(const RealType& n)
{
   return n > 2 * minimum_argument_for_bernoulli_recursion<RealType>();
}

template <class RealType, class Policy>
inline RealType stirlerr_series(const RealType& n, const Policy& pol)
{
   return bernoulli_stirling_series<RealType>(n, pol);
}

// The same for any n > 0. Below the series' range, step down from above it with
//   stirlerr(z) = stirlerr(z + 1) + (z + 1/2) log1pmx(1/z) + 1/(2z),
// whose terms are about 1/(12 z^2) and are each formed with an absolute error of about eps/z.
// (Computing ln(n!) - [(n + 1/2) ln n - n + ln(2 pi)/2] directly would cancel, with an absolute
// error of about ln(n!) eps, which grows with the series' threshold and so with the precision.)
template <class RealType, class Policy>
inline RealType stirlerr(const RealType& n, const Policy& pol)
{
   BOOST_MATH_STD_USING
   if (stirlerr_series_converges<RealType, Policy>(n))
   {
      return stirlerr_series(n, pol);
   }
   const RealType target = 2 * minimum_argument_for_bernoulli_recursion<RealType>();
   int steps = 1;
   while (!(n + steps > target))
   {
      ++steps;
   }
   RealType sum = stirlerr_series(RealType(n + steps), pol);
   for (int j = steps - 1; j >= 0; --j)
   {
      const RealType z = n + j;
      sum += (z + RealType(0.5)) * boost::math::log1pmx(RealType(1 / z), pol) + 1 / (2 * z);
   }
   return sum;
}

// The deviance term k ln(k/mean) + mean - k, for k > 0.
// For |v| < 1/2, with v = (k - mean)/(k + mean), sum its series in v^2, which has no cancellation;
// beyond that, the direct formula loses at most a factor of about 2.5.
template <class RealType, class Policy>
inline RealType bd0(const RealType& mean, const RealType& k, const Policy& pol)
{
   BOOST_MATH_STD_USING
   if (abs(k - mean) < (k + mean) / 2)
   {
      const RealType v = (k - mean) / (k + mean);
      const RealType v2 = v * v;
      RealType sum = (k - mean) * v;
      RealType term = 2 * k * v;
      const boost::math::uintmax_t max_iterations = policies::get_max_series_iterations<Policy>();
      for (boost::math::uintmax_t i = 1; i <= max_iterations; ++i)
      {
         term *= v2;
         const RealType next = term / (2 * i + 1);
         sum += next;
         if (abs(next) <= abs(sum) * tools::epsilon<RealType>())
         {
            return sum;
         }
      }
      return policies::raise_evaluation_error<RealType>("boost::math::detail::bd0<%1%>(%1%, %1%)",
          "Series did not converge, best approximation was %1%", sum, pol);
   }
   return k * log(k / mean) + mean - k;
}

// bd0 for a mean given as the ratio a / b, as when it is a product of counts over a count:
// v = (k b - a) / (k b + a) and a - k b are formed without first rounding the mean, which
// bd0's sensitivity of |k - mean| to its mean would otherwise amplify. This is exact while
// k b and a are exactly representable, e.g. below 2^53 for double.
template <class RealType, class Policy>
inline RealType bd0_ratio(const RealType& k, const RealType& a, const RealType& b, const Policy& pol)
{
   BOOST_MATH_STD_USING
   const RealType kb = k * b;
   const RealType diff = kb - a;
   if (abs(diff) < (kb + a) / 2)
   {
      const RealType v = diff / (kb + a);
      const RealType v2 = v * v;
      RealType sum = diff / b * v;
      RealType term = 2 * k * v;
      const boost::math::uintmax_t max_iterations = policies::get_max_series_iterations<Policy>();
      for (boost::math::uintmax_t i = 1; i <= max_iterations; ++i)
      {
         term *= v2;
         const RealType next = term / (2 * i + 1);
         sum += next;
         if (abs(next) <= abs(sum) * tools::epsilon<RealType>())
         {
            return sum;
         }
      }
      return policies::raise_evaluation_error<RealType>("boost::math::detail::bd0_ratio<%1%>(%1%, %1%, %1%)",
          "Series did not converge, best approximation was %1%", sum, pol);
   }
   return k * log(kb / a) - diff / b;
}

// The logarithm of the binomial density C(n, x) p^x q^(n - x), q = 1 - p, for 0 <= x <= n
// (Loader's dbinom_raw, in logarithms so that products of these can be exponentiated once).
// Its exponent is small near the mode, so near the mode it keeps its accuracy as n grows.
template <class RealType, class Policy>
inline RealType log_binomial_pdf_saddle_point(const RealType& x, const RealType& n, const RealType& p, const RealType& q, const Policy& pol)
{
   BOOST_MATH_STD_USING
   // The degenerate cases, as in R, where the general formulas would divide zero by zero.
   if ((p == 0) || (n == 0))
   {
      return x == 0 ? RealType(0) : RealType(-tools::max_value<RealType>());
   }
   if (q == 0)
   {
      return x == n ? RealType(0) : RealType(-tools::max_value<RealType>());
   }
   if (x == 0)
   {
      return p < RealType(0.1) ? RealType(-bd0(n * q, n, pol) - n * p) : RealType(n * log(q));
   }
   if (x == n)
   {
      return q < RealType(0.1) ? RealType(-bd0(n * p, n, pol) - n * q) : RealType(n * log(p));
   }
   const RealType lc = stirlerr(n, pol) - stirlerr(x, pol) - stirlerr(RealType(n - x), pol)
                     - bd0(RealType(n * p), x, pol) - bd0(RealType(n * q), RealType(n - x), pol);
   const RealType lf = log(constants::two_pi<RealType>()) + log(x) + boost::math::log1p(-x / n, pol);
   return lc - lf / 2;
}

}}} // namespaces

#endif // BOOST_MATH_HAS_GPU_SUPPORT

#endif // BOOST_MATH_DISTRIBUTIONS_DETAIL_SADDLE_POINT_HPP
