// Copyright 2008 John Maddock
//
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_DISTRIBUTIONS_DETAIL_HG_QUANTILE_HPP
#define BOOST_MATH_DISTRIBUTIONS_DETAIL_HG_QUANTILE_HPP

#include <boost/math/policies/error_handling.hpp>
#include <boost/math/distributions/detail/hypergeometric_pdf.hpp>

BOOST_MATH_NAMESPACE_BEGIN namespace detail{

//
// The result of summing the pdf up to p (or down to q = 1 - p):
//
//   up:        the smallest x with cdf(x) >= p, by the running sum,
//   nearest:   up or up - 1, whichever p is nearer to in a linear
//              interpolation of the cdf between them,
//   ambiguous: the running sum at up or up - 1 is within rounding error of
//              p, so the comparison with p must be made with cdf() itself.
//
// c says whether the caller's probability is q (the complement) or p.
//
struct hypergeometric_quantile_result
{
   std::uint64_t up;
   std::uint64_t nearest;
   bool ambiguous;
};

template <class T, class Policy>
hypergeometric_quantile_result hypergeometric_quantile_imp(T p, T q, bool c, T caller_epsilon, std::uint64_t r, std::uint64_t n, std::uint64_t N, const Policy& pol)
{
#ifdef _MSC_VER
#  pragma warning(push)
#  pragma warning(disable:4267)
#endif
   BOOST_MATH_STD_USING
   BOOST_FPU_EXCEPTION_GUARD
   // How far (relatively) the running sum may be from cdf(): the error of
   // each, doubled.  cdf() is rounded to the caller's type, which may be
   // narrower than T (float is evaluated in double):
   const T eps = (std::max)(tools::epsilon<T>(), caller_epsilon);
   const T fudge = eps * ((N <= BOOST_MATH_NAMESPACE::prime(BOOST_MATH_NAMESPACE::max_prime - 1)) ? 100 : 4 * N);
   // The caller compares with p if c is false, else with q.  When we search
   // with the other one, 1 - p or 1 - q and the cdf near 1 it is compared
   // with are only known to within a few epsilon absolutely:
   const T converted = 4 * eps;
   const std::uint64_t base = static_cast<std::uint64_t>((std::max)(0, static_cast<int>(n + r) - static_cast<int>(N)));
   const std::uint64_t lim = (std::min)(r, n);
   hypergeometric_quantile_result res;

   BOOST_MATH_INSTRUMENT_VARIABLE(p);
   BOOST_MATH_INSTRUMENT_VARIABLE(q);
   BOOST_MATH_INSTRUMENT_VARIABLE(r);
   BOOST_MATH_INSTRUMENT_VARIABLE(n);
   BOOST_MATH_INSTRUMENT_VARIABLE(N);
   BOOST_MATH_INSTRUMENT_VARIABLE(base);
   BOOST_MATH_INSTRUMENT_VARIABLE(lim);

   if(p <= 0.5)
   {
      //
      // Sum the cdf upwards: cum = cdf(x), prev = cdf(x - 1).
      //
      std::uint64_t x = base;
      T diff = hypergeometric_pdf<T>(x, r, n, N, pol);
      if (diff == 0)
      {
         ++x;
         // We want to skip through x values as fast as we can until we start getting non-zero values,
         // otherwise we're just making lots of expensive PDF calls:
         T log_pdf = BOOST_MATH_NAMESPACE::lgamma(static_cast<T>(n + 1), pol)
            + BOOST_MATH_NAMESPACE::lgamma(static_cast<T>(r + 1), pol)
            + BOOST_MATH_NAMESPACE::lgamma(static_cast<T>(N - n + 1), pol)
            + BOOST_MATH_NAMESPACE::lgamma(static_cast<T>(N - r + 1), pol)
            - BOOST_MATH_NAMESPACE::lgamma(static_cast<T>(N + 1), pol)
            - BOOST_MATH_NAMESPACE::lgamma(static_cast<T>(x + 1), pol)
            - BOOST_MATH_NAMESPACE::lgamma(static_cast<T>(n - x + 1), pol)
            - BOOST_MATH_NAMESPACE::lgamma(static_cast<T>(r - x + 1), pol)
            - BOOST_MATH_NAMESPACE::lgamma(static_cast<T>(N - n - r + x + 1), pol);
         while (log_pdf < tools::log_min_value<T>())
         {
            log_pdf += -log(static_cast<T>(x + 1)) + log(static_cast<T>(n - x)) + log(static_cast<T>(r - x)) - log(static_cast<T>(N - n - r + x + 1));
            ++x;
         }
         // By the time we get here, log_pdf may be fairly inaccurate due to
         // roundoff errors, get a fresh PDF calculation before proceeding:
         diff = hypergeometric_pdf<T>(x, r, n, N, pol);
      }
      T cum = diff;
      T prev = 0;
      while((cum < p) && (x < lim))
      {
         diff = (diff > tools::min_value<T>() * 8)
            ? T(n - x) * T(r - x) * diff / (T(x + 1) * T(N + x + 1 - n - r))
            : hypergeometric_pdf<T>(x + 1, r, n, N, pol);
         ++x;
         prev = cum;
         cum += diff;
         BOOST_MATH_INSTRUMENT_VARIABLE(x);
         BOOST_MATH_INSTRUMENT_VARIABLE(cum);
      }
      res.up = x;
      res.nearest = ((x != base) && (p < (prev + cum) / 2)) ? x - 1 : x;
      const T slack = p * fudge + (c ? converted : T(0));
      res.ambiguous = (cum < p) || (cum - p <= slack) || ((x != base) && (p - prev <= slack));
   }
   else
   {
      //
      // Sum the complement downwards: cum = cdf(complement(x)),
      // next = cdf(complement(x - 1)), and up is the smallest x with
      // cum <= q.
      //
      std::uint64_t x = lim;
      T diff = hypergeometric_pdf<T>(x, r, n, N, pol);
      if (diff == 0)
      {
         // We want to skip through x values as fast as we can until we start getting non-zero values,
         // otherwise we're just making lots of expensive PDF calls:
         --x;
         T log_pdf = BOOST_MATH_NAMESPACE::lgamma(static_cast<T>(n + 1), pol)
            + BOOST_MATH_NAMESPACE::lgamma(static_cast<T>(r + 1), pol)
            + BOOST_MATH_NAMESPACE::lgamma(static_cast<T>(N - n + 1), pol)
            + BOOST_MATH_NAMESPACE::lgamma(static_cast<T>(N - r + 1), pol)
            - BOOST_MATH_NAMESPACE::lgamma(static_cast<T>(N + 1), pol)
            - BOOST_MATH_NAMESPACE::lgamma(static_cast<T>(x + 1), pol)
            - BOOST_MATH_NAMESPACE::lgamma(static_cast<T>(n - x + 1), pol)
            - BOOST_MATH_NAMESPACE::lgamma(static_cast<T>(r - x + 1), pol)
            - BOOST_MATH_NAMESPACE::lgamma(static_cast<T>(N - n - r + x + 1), pol);
         while (log_pdf < tools::log_min_value<T>())
         {
            log_pdf += log(static_cast<T>(x)) - log(static_cast<T>(n - x + 1)) - log(static_cast<T>(r - x + 1)) + log(static_cast<T>(N - n - r + x));
            --x;
         }
         // By the time we get here, log_pdf may be fairly inaccurate due to
         // roundoff errors, get a fresh PDF calculation before proceeding:
         diff = hypergeometric_pdf<T>(x, r, n, N, pol);
      }
      // diff is pdf(x), so cdf(complement(x - 1)) = cum + diff:
      T cum = 0;
      while((cum + diff <= q) && (x > base))
      {
         cum += diff;
         diff = (diff > tools::min_value<T>() * 8)
            ? x * T(N + x - n - r) * diff / (T(1 + n - x) * T(1 + r - x))
            : hypergeometric_pdf<T>(x - 1, r, n, N, pol);
         --x;
         BOOST_MATH_INSTRUMENT_VARIABLE(x);
         BOOST_MATH_INSTRUMENT_VARIABLE(cum);
      }
      const T next = cum + diff;
      res.up = x;
      // In terms of p, x - 1 is nearer when p < cdf(x - 1) + pdf(x) / 2:
      res.nearest = ((x != base) && (q > cum + diff / 2)) ? x - 1 : x;
      const T slack = q * fudge + (c ? T(0) : converted);
      res.ambiguous = (q - cum <= slack) || ((x != base) && (next - q <= slack));
   }
   return res;
#ifdef _MSC_VER
#  pragma warning(pop)
#endif
}

template <class T, class Policy>
inline hypergeometric_quantile_result hypergeometric_quantile(T p, T q, bool c, std::uint64_t r, std::uint64_t n, std::uint64_t N, const Policy&)
{
   BOOST_FPU_EXCEPTION_GUARD
   typedef typename tools::promote_args<T>::type result_type;
   typedef typename policies::evaluation<result_type, Policy>::type value_type;
   typedef typename policies::normalise<
      Policy,
      policies::promote_float<false>,
      policies::promote_double<false>,
      policies::assert_undefined<> >::type forwarding_policy;

   return detail::hypergeometric_quantile_imp<value_type>(p, q, c, value_type(tools::epsilon<result_type>()), r, n, N, forwarding_policy());
}

} BOOST_MATH_NAMESPACE_END

#endif

