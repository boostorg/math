//  Copyright 2026 Nick Thompson
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// Building blocks for entropy(dist): Stirling remainders, which remove the cancellation between
// the x log(x) parts of lgamma and digamma at large parameters, and the Shannon entropy of a
// discrete distribution from the log ratios of its probabilities.

#ifndef BOOST_MATH_DISTRIBUTIONS_DETAIL_ENTROPY_HPP
#define BOOST_MATH_DISTRIBUTIONS_DETAIL_ENTROPY_HPP

#include <boost/math/tools/config.hpp>
#include <boost/math/tools/precision.hpp>
#include <boost/math/special_functions/bernoulli.hpp>
#include <boost/math/special_functions/digamma.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/special_functions/log1p.hpp>
#include <boost/math/special_functions/detail/two_product.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/math/policies/error_handling.hpp>
#ifndef BOOST_MATH_BUILD_MODULE
#include <algorithm>
#include <cmath>
#include <cstdint>
#endif

BOOST_MATH_NAMESPACE_BEGIN
namespace detail {

// Below this, entropies use lgamma and digamma directly; at or above it, the Stirling series
// below have converged to full precision before their terms start to grow.
template <class RealType, class Policy>
inline RealType entropy_series_threshold()
{
   return (std::max)(RealType(10), RealType(policies::digits<RealType, Policy>()) * RealType(0.14));
}

// Asymptotic series for large x of nu(x) = log(x) - 1/(2x) - digamma(x) and of
// r(x) = mu(x) + (x - 1) nu(x), where mu(x) = lgamma(x) - (x - 1/2) log(x) + x - log(2 pi)/2.
template <class RealType, class Policy>
inline void entropy_stirling_series(RealType x, RealType& r, RealType& nu, const Policy& pol)
{
   BOOST_MATH_STD_USING
   const RealType eps = policies::get_epsilon<RealType, Policy>();
   r = 0;
   nu = 0;
   RealType x_power = x; // x^(2k - 1)
   for (int k = 1; k < 100; ++k)
   {
      RealType b2k = bernoulli_b2n<RealType>(k, pol);
      RealType nu_term = b2k / (2 * k * x_power * x);
      RealType r_term = b2k / ((2 * k - 1) * x_power) - nu_term;
      r += r_term;
      nu += nu_term;
      if ((fabs(r_term) <= eps * fabs(r)) && (fabs(nu_term) <= eps * fabs(nu)))
         break;
      x_power *= x * x;
   }
}

// mu(x) = lgamma(x) - (x - 1/2) log(x) + x - log(2 pi)/2.
template <class RealType, class Policy>
inline RealType entropy_stirling_mu(RealType x, const Policy& pol)
{
   BOOST_MATH_STD_USING
   if (x < entropy_series_threshold<RealType, Policy>())
   {
      return BOOST_MATH_NAMESPACE::lgamma(x, pol) - (x - constants::half<RealType>()) * log(x) + x - constants::log_root_two_pi<RealType>();
   }
   const RealType eps = policies::get_epsilon<RealType, Policy>();
   RealType mu = 0;
   RealType x_power = x;
   for (int k = 1; k < 100; ++k)
   {
      RealType term = bernoulli_b2n<RealType>(k, pol) / (2 * k * (2 * k - 1) * x_power);
      mu += term;
      if (fabs(term) <= eps * fabs(mu))
         break;
      x_power *= x * x;
   }
   return mu;
}

// lgamma(z + t) - lgamma(z) - t log(z), for z and z + t > 0. With u = t/z it is
// z phi(u) - log1p(u)/2 + mu(z + t) - mu(z), where phi(u) = (1 + u) log1p(u) - u = u log1p(u) + log1pmx(u)
// has no cancellation to speak of, unlike the difference of the lgammas.
template <class RealType, class Policy>
inline RealType entropy_lgamma_delta(RealType z, RealType t, const Policy& pol)
{
   BOOST_MATH_STD_USING
   if (t == 0)
      return 0;
   RealType u = t / z;
   RealType l = BOOST_MATH_NAMESPACE::log1p(u, pol);
   return z * (u * l + BOOST_MATH_NAMESPACE::log1pmx(u, pol)) - l / 2 + entropy_stirling_mu(RealType(z + t), pol) - entropy_stirling_mu(z, pol);
}

// The entropy of gamma(a, 1): a + lgamma(a) + (1 - a) digamma(a).
// For large a, log(2 pi e a)/2 - 1/(2a) + r(a), with r from entropy_stirling_series.
template <class RealType, class Policy>
inline RealType gamma_entropy_unit_scale(RealType a, const Policy& pol)
{
   BOOST_MATH_STD_USING
   if (a < entropy_series_threshold<RealType, Policy>())
   {
      return a + BOOST_MATH_NAMESPACE::lgamma(a, pol) + (1 - a) * BOOST_MATH_NAMESPACE::digamma(a, pol);
   }
   RealType r, nu;
   entropy_stirling_series(a, r, nu, pol);
   return log(a) / 2 + constants::log_root_two_pi<RealType>() + constants::half<RealType>() - 1 / (2 * a) + r;
}

// The entropy of inverse_gamma(a, 1): a + lgamma(a) - (1 + a) digamma(a).
// For large a, log(2 pi e)/2 - 3 log(a)/2 + 1/(2a) + r(a) + 2 nu(a).
template <class RealType, class Policy>
inline RealType inverse_gamma_entropy_unit_scale(RealType a, const Policy& pol)
{
   BOOST_MATH_STD_USING
   if (a < entropy_series_threshold<RealType, Policy>())
   {
      return a + BOOST_MATH_NAMESPACE::lgamma(a, pol) - (1 + a) * BOOST_MATH_NAMESPACE::digamma(a, pol);
   }
   RealType r, nu;
   entropy_stirling_series(a, r, nu, pol);
   return constants::log_root_two_pi<RealType>() + constants::half<RealType>() - 3 * log(a) / 2 + 1 / (2 * a) + r + 2 * nu;
}

// a * b - c * d for a * b close to c * d: with p1 + e1 = a * b and p2 + e2 = c * d exactly (two_product,
// which also works for types whose fma isn't fused), p1 - p2 is exact by Sterbenz's lemma.
template <class RealType>
inline RealType entropy_difference_of_products(RealType a, RealType b, RealType c, RealType d)
{
   two_product<RealType> product;
   RealType p1, e1, p2, e2;
   product(a, b, p1, e1);
   product(c, d, p2, e2);
   return (p1 - p2) + (e1 - e2);
}

// Neumaier's compensated sum: the direct sums can run to ~1e5 terms (negative binomials with small p).
template <class RealType>
struct entropy_compensated_sum
{
   RealType sum = 0;
   RealType c = 0;
   void add(RealType x)
   {
      BOOST_MATH_STD_USING
      RealType t = sum + x;
      c += (fabs(sum) >= fabs(x)) ? RealType((sum - t) + x) : RealType((x - t) + sum);
      sum = t;
   }
   RealType value() const { return sum + c; }
};

// Shannon entropy -sum p log(p) of a discrete distribution, from the log ratios S = log(p / p_0) of its
// probabilities at points spaced h apart, relative to a point near the mode (S = 0 there):
//    H = log(h Z) - sum(w S) / Z,  where w = exp(S) and Z = sum(w).
// Normalising by Z means that p_0 itself is never needed. With h = 1 at the integers this is exact.
// For a smooth pmf with standard deviation sigma >> 1, the sum over the integers equals the integral of
// the pmf's continuation (through lgamma) up to about exp(-2 pi^2 sigma^2), by the Poisson summation
// formula, and the trapezoidal rule with step h = sigma / c has the same integral up to exp(-2 pi^2 c^2).
// So a few dozen points suffice for any sigma.
//
// right(S, ratio) and left(S, ratio) set S at the next point to the right (left) of the last one and
// return false when there is none (the end of the support). On the integers they also set ratio to
// p(next) / p(last), and w is then the running product of the ratios wherever that has the smaller error
// bound (j steps against |S|; it is far better for tiny Poisson means, say); elsewhere ratio = 0. rho_floor is a lower bound for the ratio of successive
// weights far out in the right tail (needed only when that ratio increases, as for negative binomials
// with r < 1); the tail is bounded as a geometric series.
template <class RealType, class Right, class Left, class Policy>
RealType discrete_entropy_sum(Right& right, Left& left, RealType h, RealType rho_floor, const char* function, const Policy& pol)
{
   BOOST_MATH_STD_USING
   const RealType eps = policies::get_epsilon<RealType, Policy>();
   const std::uintmax_t max_iter = policies::get_max_series_iterations<Policy>();
   entropy_compensated_sum<RealType> z;  // sum of w, without the w = 1 at S = 0.
   entropy_compensated_sum<RealType> zs; // sum of w S.
   std::uintmax_t iter = 0;
   for (int side = 0; side < 2; ++side)
   {
      RealType s_prev = 0;
      RealType w_product = 1;
      RealType steps = 0;
      RealType rho_min = side == 0 ? rho_floor : RealType(0);
      for (;;)
      {
         RealType s;
         RealType ratio = 0;
         if (!(side == 0 ? right(s, ratio) : left(s, ratio)))
            break;
         w_product *= ratio;
         steps += 1;
         if (++iter > max_iter)
         {
            return policies::raise_evaluation_error<RealType>(function, "Series did not converge, the entropy is approximately %1%.",
               RealType(log(h) + BOOST_MATH_NAMESPACE::log1p(z.value(), pol) - zs.value() / (1 + z.value())), pol);
         }
         RealType w = ((ratio > 0) && (steps < fabs(s))) ? w_product : RealType(exp(s));
         z.add(w);
         zs.add(w * s);
         RealType rho = (std::max)(RealType(exp(s - s_prev)), rho_min);
         s_prev = s;
         if (rho < 1)
         {
            // Bound on the remaining sum of w (1 + |S|), S growing at most linearly in the tail:
            RealType tail = w * (1 + fabs(s)) * rho / ((1 - rho) * (1 - rho));
            if (tail <= eps * (z.value() + fabs(zs.value())) / 16)
               break;
         }
      }
   }
   RealType zv = z.value();
   return log(h) + BOOST_MATH_NAMESPACE::log1p(zv, pol) - zs.value() / (1 + zv);
}

// The number of standard deviations beyond which a near-normal pmf is below eps / e^4,
// and c with exp(-2 pi^2 c^2) below eps / e^4: the trapezoidal step is sigma / c.
template <class RealType, class Policy>
inline RealType discrete_entropy_tail_sigmas()
{
   BOOST_MATH_STD_USING
   return sqrt(2 * (RealType(policies::digits<RealType, Policy>()) * constants::ln_two<RealType>() + 4));
}

template <class RealType, class Policy>
inline RealType discrete_entropy_step_ratio()
{
   BOOST_MATH_STD_USING
   return sqrt((RealType(policies::digits<RealType, Policy>()) * constants::ln_two<RealType>() + 4) / (2 * constants::pi_sqr<RealType>())) + constants::half<RealType>();
}

// log(num/den), with d = num - den computed accurately by the caller, for successive pmf ratios;
// sets ratio = num/den, or den/num going left.
template <class RealType, class Policy>
inline RealType entropy_log_ratio(RealType num, RealType den, RealType d, RealType& ratio, bool left, const Policy& pol)
{
   BOOST_MATH_STD_USING
   ratio = left ? RealType(den / num) : RealType(num / den);
   return (fabs(d) < den / 2) ? BOOST_MATH_NAMESPACE::log1p(d / den, pol) : RealType(log(num / den));
}

} // namespace detail
BOOST_MATH_NAMESPACE_END

#endif // BOOST_MATH_DISTRIBUTIONS_DETAIL_ENTROPY_HPP
