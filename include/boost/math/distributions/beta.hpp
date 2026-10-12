// boost\math\distributions\beta.hpp

// Copyright John Maddock 2006.
// Copyright Paul A. Bristow 2006.
// Copyright Matt Borland 2023.

// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

// http://en.wikipedia.org/wiki/Beta_distribution
// http://www.itl.nist.gov/div898/handbook/eda/section3/eda366h.htm
// http://mathworld.wolfram.com/BetaDistribution.html

// The Beta Distribution is a continuous probability distribution.
// The beta distribution is used to model events which are constrained to take place
// within an interval defined by maxima and minima,
// so is used extensively in PERT and other project management systems
// to describe the time to completion.
// The cdf of the beta distribution is used as a convenient way
// of obtaining the sum over a set of binomial outcomes.
// The beta distribution is also used in Bayesian statistics.

#ifndef BOOST_MATH_DIST_BETA_HPP
#define BOOST_MATH_DIST_BETA_HPP

#include <boost/math/tools/config.hpp>
#include <boost/math/tools/tuple.hpp>
#include <boost/math/distributions/fwd.hpp>
#include <boost/math/special_functions/beta.hpp> // for beta.
#include <boost/math/special_functions/digamma.hpp> // for entropy.
#include <boost/math/special_functions/bernoulli.hpp> // for entropy.
#include <boost/math/special_functions/zeta.hpp> // for entropy.
#include <boost/math/special_functions/log1p.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/math/distributions/complement.hpp> // complements.
#include <boost/math/distributions/detail/common_error_handling.hpp> // error checks
#include <boost/math/special_functions/fpclassify.hpp> // isnan.
#include <boost/math/tools/roots.hpp> // for root finding.
#include <boost/math/policies/error_handling.hpp>

#if defined (BOOST_MSVC)
#  pragma warning(push)
#  pragma warning(disable: 4702) // unreachable code
// in domain_error_imp in error_handling
#endif

BOOST_MATH_NAMESPACE_BEGIN
    namespace beta_detail
    {
      // Common error checking routines for beta distribution functions:
      template <class RealType, class Policy>
      BOOST_MATH_GPU_ENABLED inline bool check_alpha(const char* function, const RealType& alpha, RealType* result, const Policy& pol)
      {
        if(!(BOOST_MATH_NAMESPACE::isfinite)(alpha) || (alpha <= 0))
        {
          *result = policies::raise_domain_error<RealType>(
            function,
            "Alpha argument is %1%, but must be > 0 !", alpha, pol);
          return false;
        }
        return true;
      } // bool check_alpha

      template <class RealType, class Policy>
      BOOST_MATH_GPU_ENABLED inline bool check_beta(const char* function, const RealType& beta, RealType* result, const Policy& pol)
      {
        if(!(BOOST_MATH_NAMESPACE::isfinite)(beta) || (beta <= 0))
        {
          *result = policies::raise_domain_error<RealType>(
            function,
            "Beta argument is %1%, but must be > 0 !", beta, pol);
          return false;
        }
        return true;
      } // bool check_beta

      template <class RealType, class Policy>
      BOOST_MATH_GPU_ENABLED inline bool check_prob(const char* function, const RealType& p, RealType* result, const Policy& pol)
      {
        if((p < 0) || (p > 1) || !(BOOST_MATH_NAMESPACE::isfinite)(p))
        {
          *result = policies::raise_domain_error<RealType>(
            function,
            "Probability argument is %1%, but must be >= 0 and <= 1 !", p, pol);
          return false;
        }
        return true;
      } // bool check_prob

      template <class RealType, class Policy>
      BOOST_MATH_GPU_ENABLED inline bool check_x(const char* function, const RealType& x, RealType* result, const Policy& pol)
      {
        if(!(BOOST_MATH_NAMESPACE::isfinite)(x) || (x < 0) || (x > 1))
        {
          *result = policies::raise_domain_error<RealType>(
            function,
            "x argument is %1%, but must be >= 0 and <= 1 !", x, pol);
          return false;
        }
        return true;
      } // bool check_x

      template <class RealType, class Policy>
      BOOST_MATH_GPU_ENABLED inline bool check_dist(const char* function, const RealType& alpha, const RealType& beta, RealType* result, const Policy& pol)
      { // Check both alpha and beta.
        return check_alpha(function, alpha, result, pol)
          && check_beta(function, beta, result, pol);
      } // bool check_dist

      template <class RealType, class Policy>
      BOOST_MATH_GPU_ENABLED inline bool check_dist_and_x(const char* function, const RealType& alpha, const RealType& beta, RealType x, RealType* result, const Policy& pol)
      {
        return check_dist(function, alpha, beta, result, pol)
          && beta_detail::check_x(function, x, result, pol);
      } // bool check_dist_and_x

      template <class RealType, class Policy>
      BOOST_MATH_GPU_ENABLED inline bool check_dist_and_prob(const char* function, const RealType& alpha, const RealType& beta, RealType p, RealType* result, const Policy& pol)
      {
        return check_dist(function, alpha, beta, result, pol)
          && check_prob(function, p, result, pol);
      } // bool check_dist_and_prob

      template <class RealType, class Policy>
      BOOST_MATH_GPU_ENABLED inline bool check_mean(const char* function, const RealType& mean, RealType* result, const Policy& pol)
      {
        if(!(BOOST_MATH_NAMESPACE::isfinite)(mean) || (mean <= 0))
        {
          *result = policies::raise_domain_error<RealType>(
            function,
            "mean argument is %1%, but must be > 0 !", mean, pol);
          return false;
        }
        return true;
      } // bool check_mean
      template <class RealType, class Policy>
      BOOST_MATH_GPU_ENABLED inline bool check_variance(const char* function, const RealType& variance, RealType* result, const Policy& pol)
      {
        if(!(BOOST_MATH_NAMESPACE::isfinite)(variance) || (variance <= 0))
        {
          *result = policies::raise_domain_error<RealType>(
            function,
            "variance argument is %1%, but must be > 0 !", variance, pol);
          return false;
        }
        return true;
      } // bool check_variance

      // Below this, entropy uses lgamma and digamma directly; at or above it, the Stirling series
      // below has converged to full precision before its terms start to grow.
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

      // Entropy for a and b near 1, where H vanishes to second order and the usual form cancels to
      // first order. Here H = g(a) + g(b) + k(a + b), with g(x) = lgamma(x) - (x - 1) digamma(x) and
      // k(s) = (s - 2) digamma(s) - lgamma(s), which have zero value and slope at x = 1 and s = 2.
      // Integrating their derivatives, -(x - 1) trigamma(x) and (s - 2) trigamma(s), term by term
      // with the series for trigamma(1 + t) and trigamma(2 + t) gives, for |t| < 1,
      //   g(1 + t) = -sum_{n >= 0} (-1)^n (n + 1)/(n + 2) zeta(n + 2) t^(n + 2),
      //   k(2 + t) =  sum_{n >= 0} (-1)^n (n + 1)/(n + 2) (zeta(n + 2) - 1) t^(n + 2).
      // zeta(n + 2) - 1 is only accurate to absolute eps, so k uses its series for |t| <= 1/2 only;
      // further out it has no cancellation to speak of.
      template <class RealType, class Policy>
      inline RealType entropy_zeta_series(RealType t, bool minus_one, const Policy& pol)
      {
        BOOST_MATH_STD_USING
        const RealType eps = policies::get_epsilon<RealType, Policy>();
        RealType t_power = t * t;
        RealType sum = 0;
        for (int n = 0; n < 500; ++n)
        {
          RealType z = zeta(RealType(n + 2), pol);
          RealType term = RealType(n + 1) / (n + 2) * (minus_one ? z - 1 : z) * t_power;
          sum += (n & 1) ? -term : term;
          if (fabs(term) <= eps * fabs(sum))
            break;
          t_power *= t;
        }
        return sum;
      }

      template <class RealType, class Policy>
      inline RealType entropy_near_uniform(RealType a, RealType b, const Policy& pol)
      {
        BOOST_MATH_STD_USING
        RealType s = a + b;
        RealType k = (fabs(s - 2) <= RealType(0.5))
          ? entropy_zeta_series(RealType(s - 2), true, pol)
          : (s - 2) * digamma(s, pol) - lgamma(s, pol);
        auto g = [&pol](RealType x)
        {
          return (fabs(x - 1) <= RealType(0.5))
            ? RealType(-entropy_zeta_series(RealType(x - 1), false, pol))
            : RealType(lgamma(x, pol) - (x - 1) * digamma(x, pol));
        };
        return k + g(a) + g(b);
      }
    } // namespace beta_detail

    // typedef beta_distribution<double> beta;
    // is deliberately NOT included to avoid a name clash with the beta function.
    // Use beta_distribution<> mybeta(...) to construct type double.

    BOOST_MATH_EXPORT template <class RealType = double, class Policy = policies::policy<> >
    class beta_distribution
    {
    public:
      typedef RealType value_type;
      typedef Policy policy_type;

      BOOST_MATH_GPU_ENABLED beta_distribution(RealType l_alpha = 1, RealType l_beta = 1) : m_alpha(l_alpha), m_beta(l_beta)
      {
        RealType result;
        beta_detail::check_dist(
           "boost::math::beta_distribution<%1%>::beta_distribution",
          m_alpha,
          m_beta,
          &result, Policy());
      } // beta_distribution constructor.
      // Accessor functions:
      BOOST_MATH_GPU_ENABLED RealType alpha() const
      {
        return m_alpha;
      }
      BOOST_MATH_GPU_ENABLED RealType beta() const
      { // .
        return m_beta;
      }

      // Estimation of the alpha & beta parameters.
      // http://en.wikipedia.org/wiki/Beta_distribution
      // gives formulae in section on parameter estimation.
      // Also NIST EDA page 3 & 4 give the same.
      // http://www.itl.nist.gov/div898/handbook/eda/section3/eda366h.htm
      // http://www.epi.ucdavis.edu/diagnostictests/betabuster.html

      BOOST_MATH_GPU_ENABLED static RealType find_alpha(
        RealType mean, // Expected value of mean.
        RealType variance) // Expected value of variance.
      {
        constexpr auto function = "boost::math::beta_distribution<%1%>::find_alpha";
        RealType result = 0; // of error checks.
        if(false ==
            (
              beta_detail::check_mean(function, mean, &result, Policy())
              && beta_detail::check_variance(function, variance, &result, Policy())
            )
          )
        {
          return result;
        }
        return mean * (( (mean * (1 - mean)) / variance)- 1);
      } // RealType find_alpha

      BOOST_MATH_GPU_ENABLED static RealType find_beta(
        RealType mean, // Expected value of mean.
        RealType variance) // Expected value of variance.
      {
        constexpr auto function = "boost::math::beta_distribution<%1%>::find_beta";
        RealType result = 0; // of error checks.
        if(false ==
            (
              beta_detail::check_mean(function, mean, &result, Policy())
              &&
              beta_detail::check_variance(function, variance, &result, Policy())
            )
          )
        {
          return result;
        }
        return (1 - mean) * (((mean * (1 - mean)) /variance)-1);
      } //  RealType find_beta

      // Estimate alpha & beta from either alpha or beta, and x and probability.
      // Uses for these parameter estimators are unclear.

      BOOST_MATH_GPU_ENABLED static RealType find_alpha(
        RealType beta, // from beta.
        RealType x, //  x.
        RealType probability) // cdf
      {
        constexpr auto function = "boost::math::beta_distribution<%1%>::find_alpha";
        RealType result = 0; // of error checks.
        if(false ==
            (
             beta_detail::check_prob(function, probability, &result, Policy())
             &&
             beta_detail::check_beta(function, beta, &result, Policy())
             &&
             beta_detail::check_x(function, x, &result, Policy())
            )
          )
        {
          return result;
        }
        return static_cast<RealType>(ibeta_inva(beta, x, probability, Policy()));
      } // RealType find_alpha(beta, a, probability)

      BOOST_MATH_GPU_ENABLED static RealType find_beta(
        // ibeta_invb(T b, T x, T p); (alpha, x, cdf,)
        RealType alpha, // alpha.
        RealType x, // probability x.
        RealType probability) // probability cdf.
      {
        constexpr auto function = "boost::math::beta_distribution<%1%>::find_beta";
        RealType result = 0; // of error checks.
        if(false ==
            (
              beta_detail::check_prob(function, probability, &result, Policy())
              &&
              beta_detail::check_alpha(function, alpha, &result, Policy())
              &&
              beta_detail::check_x(function, x, &result, Policy())
            )
          )
        {
          return result;
        }
        return static_cast<RealType>(ibeta_invb(alpha, x, probability, Policy()));
      } //  RealType find_beta(alpha, x, probability)

    private:
      RealType m_alpha; // Two parameters of the beta distribution.
      RealType m_beta;
    }; // template <class RealType, class Policy> class beta_distribution

    #ifdef __cpp_deduction_guides
    BOOST_MATH_EXPORT template <class RealType>
    beta_distribution(RealType)->beta_distribution<typename BOOST_MATH_NAMESPACE::tools::promote_args<RealType>::type>;
    BOOST_MATH_EXPORT template <class RealType>
    beta_distribution(RealType, RealType)->beta_distribution<typename BOOST_MATH_NAMESPACE::tools::promote_args<RealType>::type>;
    #endif

    BOOST_MATH_EXPORT template <class RealType, class Policy>
    BOOST_MATH_GPU_ENABLED inline const BOOST_MATH_NAMESPACE::pair<RealType, RealType> range(const beta_distribution<RealType, Policy>& /* dist */)
    { // Range of permissible values for random variable x.
      using BOOST_MATH_NAMESPACE::tools::max_value;
      return BOOST_MATH_NAMESPACE::pair<RealType, RealType>(static_cast<RealType>(0), static_cast<RealType>(1));
    }

    BOOST_MATH_EXPORT template <class RealType, class Policy>
    BOOST_MATH_GPU_ENABLED inline const BOOST_MATH_NAMESPACE::pair<RealType, RealType> support(const beta_distribution<RealType, Policy>&  /* dist */)
    { // Range of supported values for random variable x.
      // This is range where cdf rises from 0 to 1, and outside it, the pdf is zero.
      return BOOST_MATH_NAMESPACE::pair<RealType, RealType>(static_cast<RealType>(0), static_cast<RealType>(1));
    }

    BOOST_MATH_EXPORT template <class RealType, class Policy>
    BOOST_MATH_GPU_ENABLED inline RealType mean(const beta_distribution<RealType, Policy>& dist)
    { // Mean of beta distribution = np.
      return  dist.alpha() / (dist.alpha() + dist.beta());
    } // mean

    BOOST_MATH_EXPORT template <class RealType, class Policy>
    BOOST_MATH_GPU_ENABLED inline RealType variance(const beta_distribution<RealType, Policy>& dist)
    { // Variance of beta distribution = np(1-p).
      RealType a = dist.alpha();
      RealType b = dist.beta();
      return  (a * b) / ((a + b ) * (a + b) * (a + b + 1));
    } // variance

    BOOST_MATH_EXPORT template <class RealType, class Policy>
    BOOST_MATH_GPU_ENABLED inline RealType mode(const beta_distribution<RealType, Policy>& dist)
    {
      constexpr auto function = "boost::math::mode(beta_distribution<%1%> const&)";

      RealType result;
      if ((dist.alpha() <= 1))
      {
        result = policies::raise_domain_error<RealType>(
          function,
          "mode undefined for alpha = %1%, must be > 1!", dist.alpha(), Policy());
        return result;
      }

      if ((dist.beta() <= 1))
      {
        result = policies::raise_domain_error<RealType>(
          function,
          "mode undefined for beta = %1%, must be > 1!", dist.beta(), Policy());
        return result;
      }
      RealType a = dist.alpha();
      RealType b = dist.beta();
      return (a-1) / (a + b - 2);
    } // mode

    //template <class RealType, class Policy>
    //inline RealType median(const beta_distribution<RealType, Policy>& dist)
    //{ // Median of beta distribution is not defined.
    //  return tools::domain_error<RealType>(function, "Median is not implemented, result is %1%!", std::numeric_limits<RealType>::quiet_NaN());
    //} // median

    //But WILL be provided by the derived accessor as quantile(0.5).

    BOOST_MATH_EXPORT template <class RealType, class Policy>
    BOOST_MATH_GPU_ENABLED inline RealType skewness(const beta_distribution<RealType, Policy>& dist)
    {
      BOOST_MATH_STD_USING // ADL of std functions.
      RealType a = dist.alpha();
      RealType b = dist.beta();
      return (2 * (b-a) * sqrt(a + b + 1)) / ((a + b + 2) * sqrt(a * b));
    } // skewness

    BOOST_MATH_EXPORT template <class RealType, class Policy>
    BOOST_MATH_GPU_ENABLED inline RealType kurtosis_excess(const beta_distribution<RealType, Policy>& dist)
    {
      RealType a = dist.alpha();
      RealType b = dist.beta();
      RealType a_2 = a * a;
      RealType n = 6 * (a_2 * a - a_2 * (2 * b - 1) + b * b * (b + 1) - 2 * a * b * (b + 2));
      RealType d = a * b * (a + b + 2) * (a + b + 3);
      return  n / d;
    } // kurtosis_excess

    BOOST_MATH_EXPORT template <class RealType, class Policy>
    BOOST_MATH_GPU_ENABLED inline RealType kurtosis(const beta_distribution<RealType, Policy>& dist)
    {
      return 3 + kurtosis_excess(dist);
    } // kurtosis

    // Not GPU enabled: the series use Bernoulli numbers and zeta.
    BOOST_MATH_EXPORT template <class RealType, class Policy>
    inline RealType entropy(const beta_distribution<RealType, Policy>& dist)
    {
      BOOST_MATH_STD_USING
      // H = log(B(a, b)) - (a - 1) digamma(a) - (b - 1) digamma(b) + (a + b - 2) digamma(a + b).
      // For large a or b those terms are ~a log(a) while H ~ log(a), so with s = a + b, the
      // Stirling forms of lgamma and digamma cancel the large parts exactly:
      //   H = log(2 pi e ab / s^3) / 2 - 1/(2a) - 1/(2b) + 1/s + r(a) + r(b) - r(s) + nu(s),
      // where nu and r are the remainders computed by entropy_stirling_series. A small parameter
      // b keeps its lgamma and digamma: log(2 pi e b)/2 - 1/(2b) + r(b) = lgamma(b) + b - (b - 1) digamma(b).
      RealType a = dist.alpha();
      RealType b = dist.beta();
      const RealType threshold = beta_detail::entropy_series_threshold<RealType, Policy>();
      if ((a < threshold) && (b < threshold))
      {
        return beta_detail::entropy_near_uniform(a, b, Policy());
      }
      if (a < b)
      {
        BOOST_MATH_GPU_SAFE_SWAP(a, b);
      }
      RealType s = a + b;
      RealType r_a, nu_a, r_s, nu_s;
      beta_detail::entropy_stirling_series(a, r_a, nu_a, Policy());
      beta_detail::entropy_stirling_series(s, r_s, nu_s, Policy());
      // The terms that depend on a and s only; log(a/s) = -log1p(b/a):
      RealType result = -BOOST_MATH_NAMESPACE::log1p(b / a, Policy()) / 2 - log(s) - 1 / (2 * a) + 1 / s + r_a - r_s + nu_s;
      if (b < threshold)
      {
        return result + lgamma(b, Policy()) + b - (b - 1) * digamma(b, Policy());
      }
      RealType r_b, nu_b;
      beta_detail::entropy_stirling_series(b, r_b, nu_b, Policy());
      // log(b/s) / 2 + log(s) / 2 + the b terms, given result's -log(s):
      return result + constants::log_root_two_pi<RealType>() + constants::half<RealType>()
        - BOOST_MATH_NAMESPACE::log1p(a / b, Policy()) / 2 + log(s) / 2 - 1 / (2 * b) + r_b;
    } // entropy

    BOOST_MATH_EXPORT template <class RealType, class Policy>
    BOOST_MATH_GPU_ENABLED inline RealType pdf(const beta_distribution<RealType, Policy>& dist, const RealType& x)
    { // Probability Density/Mass Function.
      BOOST_FPU_EXCEPTION_GUARD

      constexpr auto function = "boost::math::pdf(beta_distribution<%1%> const&, %1%)";

      BOOST_MATH_STD_USING // for ADL of std functions

      RealType a = dist.alpha();
      RealType b = dist.beta();

      // Argument checks:
      RealType result = 0;
      if(false == beta_detail::check_dist_and_x(
        function,
        a, b, x,
        &result, Policy()))
      {
        return result;
      }
      using BOOST_MATH_NAMESPACE::beta;

      // Corner cases: check_x ensures x element of [0, 1], but PDF is 0 for x = 0 and x = 1. PDF EQN:
      // https://wikimedia.org/api/rest_v1/media/math/render/svg/125fdaa41844a8703d1a8610ac00fbf3edacc8e7
      if(x == 0)
      {
        if (a == 1)
        {
          return static_cast<RealType>(1 / beta(a, b));
        }
        else if (a < 1)
        {
          policies::raise_overflow_error<RealType>(function, nullptr, Policy());
        }
        else
        {
          return RealType(0);
        }
      }
      else if (x == 1)
      {
        if (b == 1)
        {
          return static_cast<RealType>(1 / beta(a, b));
        }
        else if (b < 1)
        {
          policies::raise_overflow_error<RealType>(function, nullptr, Policy());
        }
        else
        {
          return RealType(0);
        }
      }
      
      return static_cast<RealType>(ibeta_derivative(a, b, x, Policy()));
    } // pdf

    BOOST_MATH_EXPORT template <class RealType, class Policy>
    BOOST_MATH_GPU_ENABLED inline RealType cdf(const beta_distribution<RealType, Policy>& dist, const RealType& x)
    { // Cumulative Distribution Function beta.
      BOOST_MATH_STD_USING // for ADL of std functions

      constexpr auto function = "boost::math::cdf(beta_distribution<%1%> const&, %1%)";

      RealType a = dist.alpha();
      RealType b = dist.beta();

      // Argument checks:
      RealType result = 0;
      if(false == beta_detail::check_dist_and_x(
        function,
        a, b, x,
        &result, Policy()))
      {
        return result;
      }
      // Special cases:
      if (x == 0)
      {
        return 0;
      }
      else if (x == 1)
      {
        return 1;
      }
      return static_cast<RealType>(ibeta(a, b, x, Policy()));
    } // beta cdf

    BOOST_MATH_EXPORT template <class RealType, class Policy>
    BOOST_MATH_GPU_ENABLED inline RealType cdf(const complemented2_type<beta_distribution<RealType, Policy>, RealType>& c)
    { // Complemented Cumulative Distribution Function beta.

      BOOST_MATH_STD_USING // for ADL of std functions

      constexpr auto function = "boost::math::cdf(beta_distribution<%1%> const&, %1%)";

      RealType const& x = c.param;
      beta_distribution<RealType, Policy> const& dist = c.dist;
      RealType a = dist.alpha();
      RealType b = dist.beta();

      // Argument checks:
      RealType result = 0;
      if(false == beta_detail::check_dist_and_x(
        function,
        a, b, x,
        &result, Policy()))
      {
        return result;
      }
      if (x == 0)
      {
        return RealType(1);
      }
      else if (x == 1)
      {
        return RealType(0);
      }
      // Calculate cdf beta using the incomplete beta function.
      // Use of ibeta here prevents cancellation errors in calculating
      // 1 - x if x is very small, perhaps smaller than machine epsilon.
      return static_cast<RealType>(ibetac(a, b, x, Policy()));
    } // beta cdf

    BOOST_MATH_EXPORT template <class RealType, class Policy>
    BOOST_MATH_GPU_ENABLED inline RealType quantile(const beta_distribution<RealType, Policy>& dist, const RealType& p)
    { // Quantile or Percent Point beta function or
      // Inverse Cumulative probability distribution function CDF.
      // Return x (0 <= x <= 1),
      // for a given probability p (0 <= p <= 1).
      // These functions take a probability as an argument
      // and return a value such that the probability that a random variable x
      // will be less than or equal to that value
      // is whatever probability you supplied as an argument.

      constexpr auto function = "boost::math::quantile(beta_distribution<%1%> const&, %1%)";

      RealType result = 0; // of argument checks:
      RealType a = dist.alpha();
      RealType b = dist.beta();
      if(false == beta_detail::check_dist_and_prob(
        function,
        a, b, p,
        &result, Policy()))
      {
        return result;
      }
      // Special cases:
      if (p == 0)
      {
        return RealType(0);
      }
      if (p == 1)
      {
        return RealType(1);
      }
      return static_cast<RealType>(ibeta_inv(a, b, p, static_cast<RealType*>(nullptr), Policy()));
    } // quantile

    BOOST_MATH_EXPORT template <class RealType, class Policy>
    BOOST_MATH_GPU_ENABLED inline RealType quantile(const complemented2_type<beta_distribution<RealType, Policy>, RealType>& c)
    { // Complement Quantile or Percent Point beta function .
      // Return the number of expected x for a given
      // complement of the probability q.

      constexpr auto function = "boost::math::quantile(beta_distribution<%1%> const&, %1%)";

      //
      // Error checks:
      RealType q = c.param;
      const beta_distribution<RealType, Policy>& dist = c.dist;
      RealType result = 0;
      RealType a = dist.alpha();
      RealType b = dist.beta();
      if(false == beta_detail::check_dist_and_prob(
        function,
        a,
        b,
        q,
        &result, Policy()))
      {
        return result;
      }
      // Special cases:
      if(q == 1)
      {
        return RealType(0);
      }
      if(q == 0)
      {
        return RealType(1);
      }

      return static_cast<RealType>(ibetac_inv(a, b, q, static_cast<RealType*>(nullptr), Policy()));
    } // Quantile Complement

  BOOST_MATH_NAMESPACE_END

// This include must be at the end, *after* the accessors
// for this distribution have been defined, in order to
// keep compilers that support two-phase lookup happy.
#include <boost/math/distributions/detail/derived_accessors.hpp>

#if defined (BOOST_MSVC)
# pragma warning(pop)
#endif

#endif // BOOST_MATH_DIST_BETA_HPP


