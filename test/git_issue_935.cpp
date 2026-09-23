// Copyright Nick Thompson, 2026
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0. (See accompanying file
// LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// See: https://github.com/boostorg/math/issues/935
// See: https://github.com/scipy/scipy/pull/17432
//
// The integer quantiles of the discrete distributions are found by solving
// the continuous cdf(x) = p and rounding.  When p is within a few ulps of
// cdf(k) the continuous root is only known to within noise, and the rounded
// result used to be off by one.  Check that quantile(cdf(k)) round-trips,
// with p nudged a few ulps either side of cdf(k), for every integer rounding
// policy, both tails, and each distribution that shares that code.

#ifndef BOOST_MATH_BUILD_MODULE
#include <boost/math/distributions/negative_binomial.hpp>
#include <boost/math/distributions/binomial.hpp>
#include <boost/math/distributions/poisson.hpp>
#else
import boost.math;
#endif

#ifndef BOOST_MATH_BUILD_MODULE
#include <cmath>
#include <vector>
#endif
#include "math_unit_test.hpp"

using namespace boost::math::policies;

enum class rounding { up, down, outwards, inwards };

// The answer each policy should give, from G[k] = cdf(k) - p (or
// p - cdf(complement(d, k)) for the complement), which increases with k.
// Returns -1 if the answer may lie beyond the tabulated k.
long expected_quantile(const std::vector<double>& G, rounding r, double lower_p)
{
   const bool up = (r == rounding::up) || (r == rounding::outwards && lower_p >= 0.5) || (r == rounding::inwards && lower_p < 0.5);
   const long n = static_cast<long>(G.size());
   if (G[0] >= 0)
      return 0;
   if (up)
   {
      // Smallest k with cdf(k) >= p, or the last of a run with cdf(k) == p:
      long k = 0;
      while ((k < n) && (G[k] < 0))
         ++k;
      if (k == n)
         return -1;
      if (G[k] == 0)
      {
         while ((k + 1 < n) && (G[k + 1] <= 0))
            ++k;
         if (k == n - 1)
            return -1;
      }
      return k;
   }
   // Largest k with cdf(k) <= p, or the first of a run with cdf(k) == p:
   long k = n - 1;
   while ((k >= 0) && (G[k] > 0))
      --k;
   if (k == n - 1)
      return -1;
   if (G[k] == 0)
   {
      while ((k > 0) && (G[k - 1] >= 0))
         --k;
   }
   return k;
}

template <class Dist>
void check_round_trips(const Dist& d, rounding r)
{
   std::vector<double> lower, upper;
   for (int k = 0; (k < 400) && (k <= support(d).second); ++k)
   {
      lower.push_back(cdf(d, static_cast<double>(k)));
      upper.push_back(cdf(complement(d, static_cast<double>(k))));
      if (lower.back() > 0.9999)
         break;
   }
   std::vector<double> G(lower.size());
   for (int comp = 0; comp < 2; ++comp)
   {
      for (std::size_t k = 0; k < lower.size(); ++k)
      {
         for (int ulps : { 0, -1, 1, -10, 10 })
         {
            double p = comp ? upper[k] : lower[k];
            for (int i = 0; i < std::abs(ulps); ++i)
               p = std::nextafter(p, ulps < 0 ? 0.0 : 1.0);
            if ((p <= 0) || (p >= 1))
               continue;
            for (std::size_t i = 0; i < lower.size(); ++i)
               G[i] = comp ? p - upper[i] : lower[i] - p;
            const long expected = expected_quantile(G, r, comp ? 1 - p : p);
            if (expected < 0)
               continue;
            const double q = comp ? quantile(complement(d, p)) : quantile(d, p);
            CHECK_EQUAL(q, static_cast<double>(expected));
         }
      }
   }
}

template <class Policy>
void test_policy(rounding r)
{
   for (double successes : { 0.5, 3.0, 5.0, 17.5 })
      for (double p : { 0.05, 0.5, 0.9 })
         check_round_trips(boost::math::negative_binomial_distribution<double, Policy>(successes, p), r);
   for (double trials : { 5.0, 100.0 })
      for (double p : { 0.1, 0.5, 0.9 })
         check_round_trips(boost::math::binomial_distribution<double, Policy>(trials, p), r);
   for (double mean : { 0.1, 4.5, 30.0 })
      check_round_trips(boost::math::poisson_distribution<double, Policy>(mean), r);
}

int main()
{
   // The case reported from SciPy: cdf(6) nudged down by 10 ulps must give 6.
   using scipy_policy = policy<discrete_quantile<integer_round_up>>;
   boost::math::negative_binomial_distribution<double, scipy_policy> nb(5.0, 0.5);
   double c = cdf(nb, 6.0);
   c -= 10 * (std::nextafter(c, 1.0) - c);
   CHECK_EQUAL(quantile(nb, c), 6.0);

   test_policy<policy<discrete_quantile<integer_round_up>>>(rounding::up);
   test_policy<policy<discrete_quantile<integer_round_down>>>(rounding::down);
   test_policy<policy<discrete_quantile<integer_round_outwards>>>(rounding::outwards);
   test_policy<policy<discrete_quantile<integer_round_inwards>>>(rounding::inwards);

   return boost::math::test::report_errors();
}
