// Copyright Nick Thompson, 2026
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0. (See accompanying file
// LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// See: https://github.com/boostorg/math/issues/1203
//
// By default the quantile of a discrete distribution is the smallest k with
// cdf(k) >= p, and quantile(complement(d, q)) the smallest k with
// cdf(complement(d, k)) <= q, as in R, SciPy, Julia and Mathematica.  The
// expected values below are from R 4.6.1 (qpois, qbinom, qnbinom, qgeom and
// qhyper); SciPy 1.18.1 and Distributions.jl 0.25.131 give the same values.

#ifndef BOOST_MATH_BUILD_MODULE
#include <boost/math/distributions/poisson.hpp>
#include <boost/math/distributions/binomial.hpp>
#include <boost/math/distributions/negative_binomial.hpp>
#include <boost/math/distributions/geometric.hpp>
#include <boost/math/distributions/hypergeometric.hpp>
#include <boost/math/distributions/bernoulli.hpp>
#else
import boost.math;
#endif

#ifndef BOOST_MATH_BUILD_MODULE
#include <cmath>
#endif
#include "math_unit_test.hpp"

template <class Dist>
void check(const Dist& d, const double (&expected)[9])
{
   const double p[] = { 0.01, 0.1, 0.25, 0.49, 0.5, 0.51, 0.75, 0.9, 0.99 };
   for (int i = 0; i < 9; ++i)
   {
      CHECK_EQUAL(quantile(d, p[i]), expected[i]);
      // The upper quantiles are the lower ones in reverse:
      CHECK_EQUAL(quantile(complement(d, p[i])), expected[8 - i]);
   }
}

int main()
{
   // The issue's example, which used to give one less for every p < 0.5:
   {
      const double expected[] = {
         2599, 2613, 2621, 2628, 2634, 2638, 2642, 2646, 2649, 2652, 2655, 2658, 2660, 2663, 2665, 2667, 2669, 2671, 2673, 2675,
         2677, 2679, 2681, 2682, 2684, 2685, 2687, 2689, 2690, 2692, 2693, 2695, 2696, 2697, 2699, 2700, 2702, 2703, 2704, 2706,
         2707, 2708, 2710, 2711, 2712, 2714, 2715, 2716, 2718, 2719, 2720, 2722, 2723, 2724, 2726, 2727, 2728, 2729, 2731, 2732,
         2734, 2735, 2736, 2738, 2739, 2740, 2742, 2743, 2745, 2746, 2748, 2749, 2751, 2753, 2754, 2756, 2758, 2759, 2761, 2763,
         2765, 2767, 2769, 2771, 2773, 2775, 2778, 2780, 2783, 2786, 2789, 2793, 2796, 2800, 2805, 2811, 2818, 2827, 2841 };
      boost::math::poisson_distribution<double> d(2719.13);
      for (int i = 1; i <= 99; ++i)
      {
         CHECK_EQUAL(quantile(d, i / 100.0), expected[i - 1]);
         CHECK_EQUAL(quantile(complement(d, i / 100.0)), expected[99 - i]);
      }
   }

   {
      const double expected[] = { 8, 11, 13, 15, 15, 15, 17, 19, 23 };
      check(boost::math::binomial_distribution<double>(50, 0.3), expected);
   }
   {
      const double expected[] = { 0, 3, 4, 7, 7, 7, 10, 13, 20 };
      check(boost::math::negative_binomial_distribution<double>(5, 0.4), expected);
   }
   {
      // The geometric used to return a real value whatever the policy:
      const double expected[] = { 0, 0, 1, 3, 3, 3, 6, 10, 20 };
      check(boost::math::geometric_distribution<double>(0.2), expected);
   }
   {
      // 25 defective items in 60, sample 20:
      const double expected[] = { 4, 6, 7, 8, 8, 8, 10, 11, 13 };
      check(boost::math::hypergeometric_distribution<double>(25, 20, 60), expected);
   }

   // p exactly on a step: quantile(cdf(k)) == k, and one ulp above gives k + 1.
   {
      boost::math::poisson_distribution<double> d(2719.13);
      boost::math::geometric_distribution<double> g(0.2);
      boost::math::hypergeometric_distribution<double> h(25, 20, 60);
      for (double k : { 2600.0, 2718.0, 2719.0, 2840.0 })
      {
         CHECK_EQUAL(quantile(d, cdf(d, k)), k);
         CHECK_EQUAL(quantile(d, std::nextafter(cdf(d, k), 1.0)), k + 1);
         CHECK_EQUAL(quantile(complement(d, cdf(complement(d, k)))), k);
         CHECK_EQUAL(quantile(complement(d, std::nextafter(cdf(complement(d, k)), 0.0))), k + 1);
      }
      for (double k : { 1.0, 4.0, 13.0 })
      {
         CHECK_EQUAL(quantile(g, cdf(g, k)), k);
         CHECK_EQUAL(quantile(complement(g, cdf(complement(g, k)))), k);
      }
      for (double k : { 4.0, 8.0, 12.0 })
      {
         CHECK_EQUAL(quantile(h, cdf(h, k)), k);
         CHECK_EQUAL(quantile(h, std::nextafter(cdf(h, k), 1.0)), k + 1);
         CHECK_EQUAL(quantile(complement(h, cdf(complement(h, k)))), k);
         CHECK_EQUAL(quantile(complement(h, std::nextafter(cdf(complement(h, k)), 0.0))), k + 1);
      }
   }

   // Near p = 1 the cdf can round to the same value at consecutive k.  The
   // smallest k with cdf(k) >= p is then the first of that run.
   {
      boost::math::poisson_distribution<double> d(4.5);
      for (double k = 20; k < 60; ++k)
      {
         const double p = cdf(d, k);
         if ((p < 1) && (cdf(d, k - 1) < p))
            CHECK_EQUAL(quantile(d, p), k);
      }
   }

   // A long run: with a small success fraction the cdf rounds to the same
   // value just below 1 for a thousand or so consecutive k.
   {
      using namespace boost::math::policies;
      boost::math::geometric_distribution<double> up(0.001);
      boost::math::geometric_distribution<double, policy<discrete_quantile<integer_round_down>>> down(0.001);
      const double p = std::nextafter(1.0, 0.0);
      const double ku = quantile(up, p);
      CHECK_EQUAL(cdf(up, ku) >= p, true);
      CHECK_EQUAL(cdf(up, ku - 1) < p, true);
      const double kd = quantile(down, p);
      CHECK_EQUAL(cdf(down, kd) <= p, true);
      CHECK_EQUAL(cdf(down, kd + 1) > p, true);
   }

   // q one ulp below cdf(complement(d, 0)): the answer is 1, not 0.
   {
      boost::math::geometric_distribution<double> d(0.3);
      const double q = std::nextafter(cdf(complement(d, 0.0)), 0.0);
      CHECK_EQUAL(quantile(complement(d, q)), 1.0);
      const double p = std::nextafter(cdf(d, 0.0), 1.0);
      CHECK_EQUAL(quantile(d, p), 1.0);
   }

   // The Bernoulli has its own quantile, with the same definition:
   for (double sf : { 0.25, 0.5, 0.75 })
   {
      boost::math::bernoulli_distribution<double> d(sf);
      for (double k : { 0.0, 1.0 })
      {
         CHECK_EQUAL(quantile(d, cdf(d, k)), k);
         CHECK_EQUAL(quantile(complement(d, cdf(complement(d, k)))), k);
      }
      // P(X > 0) == sf, so q >= sf gives 0 and q < sf gives 1:
      CHECK_EQUAL(quantile(complement(d, std::nextafter(sf, 1.0))), 0.0);
      CHECK_EQUAL(quantile(complement(d, std::nextafter(sf, 0.0))), 1.0);
      CHECK_EQUAL(quantile(d, std::nextafter(1 - sf, 1.0)), 1.0);
   }

#ifndef BOOST_NO_EXCEPTIONS
   // Running out of iterations while settling the result is an error, not a
   // silently returned guess:
   {
      using namespace boost::math::policies;
      boost::math::geometric_distribution<double, policy<max_root_iterations<0>>> d(0.5);
      bool thrown = false;
      try
      {
         (void)quantile(d, 0.75);  // Has to look at k = 0 as well as k = 1.
      }
      catch (const boost::math::evaluation_error&)
      {
         thrown = true;
      }
      CHECK_EQUAL(thrown, true);
   }
#endif

   return boost::math::test::report_errors();
}
