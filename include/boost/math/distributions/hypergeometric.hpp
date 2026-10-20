// Copyright 2008 Gautam Sewani
// Copyright 2008 John Maddock
// Copyright 2021 Paul A. Bristow
//
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_DISTRIBUTIONS_HYPERGEOMETRIC_HPP
#define BOOST_MATH_DISTRIBUTIONS_HYPERGEOMETRIC_HPP

#include <boost/math/distributions/detail/common_error_handling.hpp>
#include <boost/math/distributions/complement.hpp>
#include <boost/math/distributions/detail/hypergeometric_pdf.hpp>
#include <boost/math/distributions/detail/hypergeometric_cdf.hpp>
#include <boost/math/distributions/detail/hypergeometric_quantile.hpp>
#include <boost/math/distributions/detail/inv_discrete_quantile.hpp>
#include <boost/math/special_functions/fpclassify.hpp>
#ifndef BOOST_MATH_BUILD_MODULE
#include <cstdint>
#endif

#include <boost/math/distributions/detail/entropy.hpp>

BOOST_MATH_NAMESPACE_BEGIN

   BOOST_MATH_EXPORT template <class RealType = double, class Policy = policies::policy<> >
   class hypergeometric_distribution
   {
   public:
      typedef RealType value_type;
      typedef Policy policy_type;

      hypergeometric_distribution(std::uint64_t r, std::uint64_t n, std::uint64_t N) // Constructor. r=defective/failures/success, n=trials/draws, N=total population.
         : m_n(n), m_N(N), m_r(r)
      {
         static const char* function = "boost::math::hypergeometric_distribution<%1%>::hypergeometric_distribution";
         RealType ret;
         check_params(function, &ret);
      }
      // Accessor functions.
      std::uint64_t total() const
      {
         return m_N;
      }

      std::uint64_t defective() const // successes/failures/events
      {
         return m_r;
      }

      std::uint64_t sample_count()const
      {
         return m_n;
      }

      bool check_params(const char* function, RealType* result)const
      {
         if(m_r > m_N)
         {
            *result = BOOST_MATH_NAMESPACE::policies::raise_domain_error<RealType>(
               function, "Parameter r out of range: must be <= N but got %1%", static_cast<RealType>(m_r), Policy());
            return false;
         }
         if(m_n > m_N)
         {
            *result = BOOST_MATH_NAMESPACE::policies::raise_domain_error<RealType>(
               function, "Parameter n out of range: must be <= N but got %1%", static_cast<RealType>(m_n), Policy());
            return false;
         }
         return true;
      }
      bool check_x(std::uint64_t x, const char* function, RealType* result)const
      {
         if(x < static_cast<std::uint64_t>((std::max)(INT64_C(0), static_cast<std::int64_t>(m_n + m_r) - static_cast<std::int64_t>(m_N))))
         {
            *result = BOOST_MATH_NAMESPACE::policies::raise_domain_error<RealType>(
               function, "Random variable out of range: must be > 0 and > m + r - N but got %1%", static_cast<RealType>(x), Policy());
            return false;
         }
         if(x > (std::min)(m_r, m_n))
         {
            *result = BOOST_MATH_NAMESPACE::policies::raise_domain_error<RealType>(
               function, "Random variable out of range: must be less than both n and r but got %1%", static_cast<RealType>(x), Policy());
            return false;
         }
         return true;
      }

   private:
      // Data members:
      std::uint64_t m_n;  // number of items picked or drawn.
      std::uint64_t m_N; // number of "total" items.
      std::uint64_t m_r; // number of "defective/successes/failures/events items.

   }; // class hypergeometric_distribution

   BOOST_MATH_EXPORT typedef hypergeometric_distribution<double> hypergeometric;

   BOOST_MATH_EXPORT template <class RealType, class Policy>
   inline const std::pair<std::uint64_t, std::uint64_t> range(const hypergeometric_distribution<RealType, Policy>& dist)
   { // Range of permissible values for random variable x.
#ifdef _MSC_VER
#  pragma warning(push)
#  pragma warning(disable:4267)
#endif
      const auto r = dist.defective();
      const auto n = dist.sample_count();
      const auto N = dist.total();
      const auto l = static_cast<std::uint64_t>((std::max)(INT64_C(0), static_cast<std::int64_t>(n + r) - static_cast<std::int64_t>(N)));
      const auto u = (std::min)(r, n);
      return std::make_pair(l, u);
#ifdef _MSC_VER
#  pragma warning(pop)
#endif
   }

   BOOST_MATH_EXPORT template <class RealType, class Policy>
   inline const std::pair<std::uint64_t, std::uint64_t> support(const hypergeometric_distribution<RealType, Policy>& d)
   {
      return range(d);
   }

   BOOST_MATH_EXPORT template <class RealType, class Policy>
   inline RealType pdf(const hypergeometric_distribution<RealType, Policy>& dist, const std::uint64_t& x)
   {
      static const char* function = "boost::math::pdf(const hypergeometric_distribution<%1%>&, const %1%&)";
      RealType result = 0;
      if(!dist.check_params(function, &result))
         return result;
      if(!dist.check_x(x, function, &result))
         return result;

      return BOOST_MATH_NAMESPACE::detail::hypergeometric_pdf<RealType>(
         x, dist.defective(), dist.sample_count(), dist.total(), Policy());
   }

   BOOST_MATH_EXPORT template <class RealType, class Policy, class U>
   inline RealType pdf(const hypergeometric_distribution<RealType, Policy>& dist, const U& x)
   {
      BOOST_MATH_STD_USING
      static const char* function = "boost::math::pdf(const hypergeometric_distribution<%1%>&, const %1%&)";
      RealType r = static_cast<RealType>(x);
      auto u = static_cast<std::uint64_t>(lltrunc(r, typename policies::normalise<Policy, policies::rounding_error<policies::ignore_error> >::type()));
      if(u == r)
      {
         return pdf(dist, u);
      }
      // Non-integer x: the analytic continuation of C(r, x) C(N - r, n - x) / C(N, n).
      RealType result = 0;
      if(!dist.check_params(function, &result))
         return result;
      const RealType lower = static_cast<RealType>((std::max)(INT64_C(0), static_cast<std::int64_t>(dist.sample_count() + dist.defective()) - static_cast<std::int64_t>(dist.total())));
      const RealType upper = static_cast<RealType>((std::min)(dist.defective(), dist.sample_count()));
      if(!(BOOST_MATH_NAMESPACE::isfinite)(r) || !(r > lower) || !(r < upper))
      {
         return BOOST_MATH_NAMESPACE::policies::raise_domain_error<RealType>(
            function, "Random variable out of range: must be between max(0, n + r - N) and min(n, r) but got %1%", r, Policy());
      }
      // Here x > 0, so its truncation u is floor(x).
      return BOOST_MATH_NAMESPACE::detail::hypergeometric_pdf_noninteger(
         r, u, dist.defective(), dist.sample_count(), dist.total(), Policy());
   }

   BOOST_MATH_EXPORT template <class RealType, class Policy>
   inline RealType cdf(const hypergeometric_distribution<RealType, Policy>& dist, const std::uint64_t& x)
   {
      static const char* function = "boost::math::cdf(const hypergeometric_distribution<%1%>&, const %1%&)";
      RealType result = 0;
      if(!dist.check_params(function, &result))
         return result;
      if(!dist.check_x(x, function, &result))
         return result;

      return BOOST_MATH_NAMESPACE::detail::hypergeometric_cdf<RealType>(
         x, dist.defective(), dist.sample_count(), dist.total(), false, Policy());
   }

   BOOST_MATH_EXPORT template <class RealType, class Policy, class U>
   inline RealType cdf(const hypergeometric_distribution<RealType, Policy>& dist, const U& x)
   {
      BOOST_MATH_STD_USING
      static const char* function = "boost::math::cdf(const hypergeometric_distribution<%1%>&, const %1%&)";
      RealType r = static_cast<RealType>(x);
      auto u = static_cast<std::uint64_t>(lltrunc(r, typename policies::normalise<Policy, policies::rounding_error<policies::ignore_error> >::type()));
      if(u != r)
      {
         return BOOST_MATH_NAMESPACE::policies::raise_domain_error<RealType>(
            function, "Random variable out of range: must be an integer but got %1%", r, Policy());
      }
      return cdf(dist, u);
   }

   BOOST_MATH_EXPORT template <class RealType, class Policy>
   inline RealType cdf(const complemented2_type<hypergeometric_distribution<RealType, Policy>, std::uint64_t>& c)
   {
      static const char* function = "boost::math::cdf(const hypergeometric_distribution<%1%>&, const %1%&)";
      RealType result = 0;
      if(!c.dist.check_params(function, &result))
         return result;
      if(!c.dist.check_x(c.param, function, &result))
         return result;

      return BOOST_MATH_NAMESPACE::detail::hypergeometric_cdf<RealType>(
         c.param, c.dist.defective(), c.dist.sample_count(), c.dist.total(), true, Policy());
   }

   BOOST_MATH_EXPORT template <class RealType, class Policy, class U>
   inline RealType cdf(const complemented2_type<hypergeometric_distribution<RealType, Policy>, U>& c)
   {
      BOOST_MATH_STD_USING
      static const char* function = "boost::math::cdf(const hypergeometric_distribution<%1%>&, const %1%&)";
      RealType r = static_cast<RealType>(c.param);
      auto u = static_cast<std::uint64_t>(lltrunc(r, typename policies::normalise<Policy, policies::rounding_error<policies::ignore_error> >::type()));
      if(u != r)
      {
         return BOOST_MATH_NAMESPACE::policies::raise_domain_error<RealType>(
            function, "Random variable out of range: must be an integer but got %1%", r, Policy());
      }
      return cdf(complement(c.dist, u));
   }

   namespace detail {
   //
   // hypergeometric_quantile sums the pdf to find the answer.  Away from a
   // step of the cdf that sum decides it; within rounding error of one the
   // comparison with p is made using cdf() itself, as for the other discrete
   // distributions, so that quantile(cdf(k)) == k.
   //
   template <class Dist>
   inline typename Dist::value_type hypergeometric_round_quantile(const Dist& d, const hypergeometric_quantile_result& s, typename Dist::value_type p, bool c, const policies::discrete_quantile<policies::integer_round_up>& tag)
   {
      typedef typename Dist::value_type value_type;
      return s.ambiguous ? round_discrete_quantile(d, value_type(s.up), p, c, tag) : value_type(s.up);
   }
   template <class Dist>
   inline typename Dist::value_type hypergeometric_round_quantile(const Dist& d, const hypergeometric_quantile_result& s, typename Dist::value_type p, bool c, const policies::discrete_quantile<policies::integer_round_down>& tag)
   {
      typedef typename Dist::value_type value_type;
      if (s.ambiguous)
         return round_discrete_quantile(d, value_type(value_type(s.up) - 1), p, c, tag);
      // cdf(up) > p clearly, so the answer is the integer below it, if any:
      return s.up > support(d).first ? value_type(s.up - 1) : value_type(s.up);
   }
   template <class Dist>
   inline typename Dist::value_type hypergeometric_round_quantile(const Dist& d, const hypergeometric_quantile_result& s, typename Dist::value_type p, bool c, const policies::discrete_quantile<policies::integer_round_outwards>&)
   {
      return is_lower_quantile(p, c)
         ? hypergeometric_round_quantile(d, s, p, c, policies::discrete_quantile<policies::integer_round_down>())
         : hypergeometric_round_quantile(d, s, p, c, policies::discrete_quantile<policies::integer_round_up>());
   }
   template <class Dist>
   inline typename Dist::value_type hypergeometric_round_quantile(const Dist& d, const hypergeometric_quantile_result& s, typename Dist::value_type p, bool c, const policies::discrete_quantile<policies::integer_round_inwards>&)
   {
      return is_lower_quantile(p, c)
         ? hypergeometric_round_quantile(d, s, p, c, policies::discrete_quantile<policies::integer_round_up>())
         : hypergeometric_round_quantile(d, s, p, c, policies::discrete_quantile<policies::integer_round_down>());
   }
   template <class Dist>
   inline typename Dist::value_type hypergeometric_round_quantile(const Dist& d, const hypergeometric_quantile_result& s, typename Dist::value_type p, bool c, const policies::discrete_quantile<policies::integer_round_nearest>& tag)
   {
      // The nearer of up - 1 and up, as for the other distributions (which
      // also never return more than up):
      typedef typename Dist::value_type value_type;
      return s.ambiguous ? round_discrete_quantile(d, value_type(s.nearest), p, c, tag) : value_type(s.nearest);
   }
   } // namespace detail

   BOOST_MATH_EXPORT template <class RealType, class Policy>
   inline RealType quantile(const hypergeometric_distribution<RealType, Policy>& dist, const RealType& p)
   {
      BOOST_MATH_STD_USING // for ADL of std functions

      // Checking function argument
      RealType result = 0;
      const char* function = "boost::math::quantile(const hypergeometric_distribution<%1%>&, %1%)";
      if (false == dist.check_params(function, &result)) 
         return result;

      if(false == detail::check_probability(function, p, &result, Policy())) 
         return result;

      // Only the ends of the support have cdf == 0 or 1 exactly, but the
      // computed cdf can round to 0 or 1 further in, so handle them here:
      if (p == 0)
         return static_cast<RealType>(support(dist).first);
      if (p == 1)
         return static_cast<RealType>(support(dist).second);
      return detail::hypergeometric_round_quantile(dist, detail::hypergeometric_quantile(p, RealType(1 - p), false, dist.defective(), dist.sample_count(), dist.total(), Policy()), p, false, typename Policy::discrete_quantile_type());
   } // quantile

   BOOST_MATH_EXPORT template <class RealType, class Policy>
   inline RealType quantile(const complemented2_type<hypergeometric_distribution<RealType, Policy>, RealType>& c)
   {
      BOOST_MATH_STD_USING // for ADL of std functions

      // Checking function argument
      RealType result = 0;
      const char* function = "quantile(const complemented2_type<hypergeometric_distribution<%1%>, %1%>&)";
      if (false == c.dist.check_params(function, &result)) 
         return result;
      if (false == detail::check_probability(function, c.param, &result, Policy()))
         return result;

      // As above:
      if (c.param == 1)
         return static_cast<RealType>(support(c.dist).first);
      if (c.param == 0)
         return static_cast<RealType>(support(c.dist).second);
      return detail::hypergeometric_round_quantile(c.dist, detail::hypergeometric_quantile(RealType(1 - c.param), c.param, true, c.dist.defective(), c.dist.sample_count(), c.dist.total(), Policy()), c.param, true, typename Policy::discrete_quantile_type());
   } // quantile

   // https://www.wolframalpha.com/input/?i=kurtosis+hypergeometric+distribution

   BOOST_MATH_EXPORT template <class RealType, class Policy>
   inline RealType mean(const hypergeometric_distribution<RealType, Policy>& dist)
   {
      return static_cast<RealType>(dist.defective() * dist.sample_count()) / dist.total();
   } // RealType mean(const hypergeometric_distribution<RealType, Policy>& dist)

   BOOST_MATH_EXPORT template <class RealType, class Policy>
   inline RealType variance(const hypergeometric_distribution<RealType, Policy>& dist)
   {
      RealType r = static_cast<RealType>(dist.defective());
      RealType n = static_cast<RealType>(dist.sample_count());
      RealType N = static_cast<RealType>(dist.total());
      return n * r  * (N - r) * (N - n) / (N * N * (N - 1));
   } // RealType variance(const hypergeometric_distribution<RealType, Policy>& dist)

   BOOST_MATH_EXPORT template <class RealType, class Policy>
   inline RealType mode(const hypergeometric_distribution<RealType, Policy>& dist)
   {
      BOOST_MATH_STD_USING
      RealType r = static_cast<RealType>(dist.defective());
      RealType n = static_cast<RealType>(dist.sample_count());
      RealType N = static_cast<RealType>(dist.total());
      return floor((r + 1) * (n + 1) / (N + 2));
   }

   BOOST_MATH_EXPORT template <class RealType, class Policy>
   inline RealType skewness(const hypergeometric_distribution<RealType, Policy>& dist)
   {
      BOOST_MATH_STD_USING
      RealType r = static_cast<RealType>(dist.defective());
      RealType n = static_cast<RealType>(dist.sample_count());
      RealType N = static_cast<RealType>(dist.total());
      return (N - 2 * r) * sqrt(N - 1) * (N - 2 * n) / (sqrt(n * r * (N - r) * (N - n)) * (N - 2));
   } // RealType skewness(const hypergeometric_distribution<RealType, Policy>& dist)

   BOOST_MATH_EXPORT template <class RealType, class Policy>
   inline RealType kurtosis_excess(const hypergeometric_distribution<RealType, Policy>& dist)
   {
      // https://www.wolframalpha.com/input/?i=kurtosis+hypergeometric+distribution shown as plain text:
      //  mean | (m n)/N
      //  standard deviation | sqrt((m n(N - m) (N - n))/(N - 1))/N
      //  variance | (m n(1 - m/N) (N - n))/((N - 1) N)
      //  skewness | (sqrt(N - 1) (N - 2 m) (N - 2 n))/((N - 2) sqrt(m n(N - m) (N - n)))
      //  kurtosis | ((N - 1) N^2 ((3 m(N - m) (n^2 (-N) + (n - 2) N^2 + 6 n(N - n)))/N^2 - 6 n(N - n) + N(N + 1)))/(m n(N - 3) (N - 2) (N - m) (N - n))
     // Kurtosis[HypergeometricDistribution[n, m, N]]
      RealType m = static_cast<RealType>(dist.defective()); // Failures or success events. (Also symbols K or M are used).
      RealType n = static_cast<RealType>(dist.sample_count()); // draws or trials.
      RealType n2 = n * n; // n^2
      RealType N = static_cast<RealType>(dist.total()); // Total population from which n draws or trials are made.
      RealType N2 = N * N; // N^2
      // result = ((N - 1) N^2 ((3 m(N - m) (n^2 (-N) + (n - 2) N^2 + 6 n(N - n)))/N^2 - 6 n(N - n) + N(N + 1)))/(m n(N - 3) (N - 2) (N - m) (N - n));
      RealType result = ((N-1)*N2*((3*m*(N-m)*(n2*(-N)+(n-2)*N2+6*n*(N-n)))/N2-6*n*(N-n)+N*(N+1)))/(m*n*(N-3)*(N-2)*(N-m)*(N-n));
      // Agrees with kurtosis hypergeometric distribution(50,200,500) kurtosis = 2.96917
      // N[kurtosis[hypergeometricdistribution(50,200,500)], 55]  2.969174035736058474901169623721804275002985337280263464
      return result;
   } // RealType kurtosis_excess(const hypergeometric_distribution<RealType, Policy>& dist)

   BOOST_MATH_EXPORT template <class RealType, class Policy>
   inline RealType kurtosis(const hypergeometric_distribution<RealType, Policy>& dist)
   {
      return kurtosis_excess(dist) + 3;
   } // RealType kurtosis_excess(const hypergeometric_distribution<RealType, Policy>& dist)
   // Not GPU enabled: the large-variance path uses Bernoulli numbers.
   BOOST_MATH_EXPORT template <class RealType, class Policy>
   inline RealType entropy(const hypergeometric_distribution<RealType, Policy>& dist)
   {
      BOOST_MATH_STD_USING
      // No closed form: -sum p log(p) from the log ratios of the probabilities (detail::discrete_entropy_sum).
      constexpr auto function = "boost::math::entropy(const hypergeometric_distribution<%1%>&)";
      RealType result = 0;
      if (!dist.check_params(function, &result))
         return result;
      const std::uint64_t r = dist.defective();
      const std::uint64_t n = dist.sample_count();
      const std::uint64_t N = dist.total();
      const std::uint64_t k_min = (n + r > N) ? n + r - N : 0;
      const std::uint64_t k_max = (std::min)(n, r);
      if (k_min == k_max)
         return 0;
      const RealType rr = static_cast<RealType>(r);
      const RealType nn = static_cast<RealType>(n);
      const RealType NN = static_cast<RealType>(N);
      // The mode, floor((r + 1)(n + 1) / (N + 2)); only the starting point of the sums:
      std::uint64_t mode = (r + 1 <= UINT64_MAX / (n + 1))
         ? ((r + 1) * (n + 1)) / (N + 2)
         : static_cast<std::uint64_t>(std::floor(static_cast<double>(r + 1) * static_cast<double>(n + 1) / static_cast<double>(N + 2)));
      mode = (std::max)(k_min, (std::min)(k_max, mode));
      const RealType sigmas = detail::discrete_entropy_tail_sigmas<RealType, Policy>();
      const RealType threshold = detail::entropy_series_threshold<RealType, Policy>();
      RealType sigma = sqrt(nn * (rr / NN) * (1 - rr / NN) * (NN - nn) / (NN - 1));
      // log p(k) = -lgamma(k + 1) - lgamma(r - k + 1) - lgamma(n - k + 1) - lgamma(N - r - n + k + 1) + const:
      const RealType z1 = static_cast<RealType>(mode + 1);
      const RealType z2 = static_cast<RealType>(r - mode + 1);
      const RealType z3 = static_cast<RealType>(n - mode + 1);
      const RealType z4 = static_cast<RealType>(N + mode + 1 - r - n);
      RealType margin = sigmas * sigma + threshold + 2;
      if ((sigma >= 2 * sigmas) && ((std::min)((std::min)(z1, z2), (std::min)(z3, z4)) >= margin))
      {
         // Trapezoidal rule on the continuation, with linear coefficient alpha = log(z2 z3 / (z1 z4)):
         RealType alpha = BOOST_MATH_NAMESPACE::log1p(RealType(detail::difference_of_products(z2, z3, z1, z4) / (z1 * z4)), Policy());
         RealType h = sigma / detail::discrete_entropy_step_ratio<RealType, Policy>();
         auto s_of = [&](RealType t)
         {
            return t * alpha - detail::entropy_lgamma_delta(z1, t, Policy()) - detail::entropy_lgamma_delta(z2, RealType(-t), Policy())
               - detail::entropy_lgamma_delta(z3, RealType(-t), Policy()) - detail::entropy_lgamma_delta(z4, t, Policy());
         };
         int j_right = 0;
         int j_left = 0;
         auto right = [&](RealType& s, RealType&)
         {
            RealType t = (++j_right) * h;
            if ((std::min)(z2, z3) - t < threshold)
               return false;
            s = s_of(t);
            return true;
         };
         auto left = [&](RealType& s, RealType&)
         {
            RealType t = -(++j_left) * h;
            if ((std::min)(z1, z4) + t < threshold)
               return false;
            s = s_of(t);
            return true;
         };
         return detail::discrete_entropy_sum(right, left, h, RealType(0), function, Policy());
      }
      // Direct sum from the mode, with p(k + 1) / p(k) = (r - k)(n - k) / ((k + 1)(N - r - n + k + 1)),
      // whose numerator minus denominator is (r + 1)(n + 1) - (k + 1)(N + 2):
      std::uint64_t k_right = mode;
      std::uint64_t k_left = mode;
      detail::compensated_sum<RealType> s_right;
      detail::compensated_sum<RealType> s_left;
      auto log_ratio = [&](std::uint64_t k, RealType& ratio, bool left)
      {
         RealType kk = static_cast<RealType>(k);
         RealType num = (rr - kk) * (nn - kk);
         RealType den = (kk + 1) * static_cast<RealType>(N + k + 1 - r - n);
         RealType d = detail::difference_of_products(RealType(rr + 1), RealType(nn + 1), RealType(kk + 1), RealType(NN + 2));
         return detail::entropy_log_ratio(num, den, d, ratio, left, Policy());
      };
      auto right = [&](RealType& s, RealType& ratio)
      {
         if (k_right == k_max)
            return false;
         s_right.add(log_ratio(k_right, ratio, false));
         k_right += 1;
         s = s_right.value();
         return true;
      };
      auto left = [&](RealType& s, RealType& ratio)
      {
         if (k_left == k_min)
            return false;
         k_left -= 1;
         s_left.add(-log_ratio(k_left, ratio, true));
         s = s_left.value();
         return true;
      };
      return detail::discrete_entropy_sum(right, left, RealType(1), RealType(0), function, Policy());
   }

BOOST_MATH_NAMESPACE_END

// This include must be at the end, *after* the accessors
// for this distribution have been defined, in order to
// keep compilers that support two-phase lookup happy.
#include <boost/math/distributions/detail/derived_accessors.hpp>

#endif // include guard
