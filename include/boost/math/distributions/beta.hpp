// boost\math\distributions\beta.hpp

// Copyright John Maddock 2006.
// Copyright Paul A. Bristow 2006.

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

#ifndef BOOST_MATH_DIST_BETA_HPP
#define BOOST_MATH_DIST_BETA_HPP

#include <boost/math/special_functions/beta.hpp> // for beta.
#include <boost/math/distributions/complement.hpp> // complements.
#include <boost/math/distributions/detail/common_error_handling.hpp> // error checks
#include <boost/math/special_functions/fpclassify.hpp> // isnan.
#include <boost/math/tools/roots.hpp> // for root finding.

#if defined (BOOST_MSVC) && defined(BOOST_MATH_THROW_ON_DOMAIN_ERROR)
#  pragma warning(push)
#  pragma warning(disable: 4702) // unreachable code
// in domain_error_imp in error_handling
#endif

namespace boost
{
  namespace math
  {
    namespace beta_detail
    {
      // Common error checking routines for beta distribution functions:
      template <class RealType>
      inline bool check_alpha(const char* function, const RealType& alpha, RealType* result)
      {
        if(!(boost::math::isfinite)(alpha) || (alpha < 0))
        {
          *result = tools::domain_error<RealType>(
            function,
            "Alpha argument is %1%, but must be >= 0 !", alpha);
          return false;
        }
        return true;
      } // bool check_alpha
      template <class RealType>
      inline bool check_beta(const char* function, const RealType& beta, RealType* result)
      {
        if(!(boost::math::isfinite)(beta) || (beta < 0))
        {
          *result = tools::domain_error<RealType>(
            function,
            "Beta argument is %1%, but must be >= 0 !", beta);
          return false;
        }
        return true;
      } // bool check_beta
      template <class RealType>
      inline bool check_prob(const char* function, const RealType& p, RealType* result)
      {
        if((p < 0) || (p > 1) || !(boost::math::isfinite)(p))
        {
          *result = tools::domain_error<RealType>(
            function,
            "Probability argument is %1%, but must be >= 0 and <= 1 !", p);
          return false;
        }
        return true;
      } // bool check_prob
      template <class RealType>
      inline bool check_x(const char* function, const RealType& x, RealType* result)
      {
        if(!(boost::math::isfinite)(x) || (x < 0) || (x > 1))
        {
          *result = tools::domain_error<RealType>(
            function,
            "x argument is %1%, but must be >= 0 and <= 1 !", x);
          return false;
        }
        return true;
      } // bool check_x
      template <class RealType>
      inline bool check_dist(const char* function, const RealType& alpha, const RealType& beta, RealType* result)
      { // Check both alpha and beta.
        return check_alpha(function, alpha, result)
          && check_beta(function, beta, result);
      } // bool check_dist
      template <class RealType>
      inline bool check_dist_and_x(const char* function, const RealType& alpha, const RealType& beta, RealType x, RealType* result)
      {
        return check_dist(function, alpha, beta, result) 
          && check_x(function, x, result);
      } // bool check_dist_and_x
      template <class RealType>
      inline bool check_dist_and_prob(const char* function, const RealType& alpha, const RealType& beta, RealType p, RealType* result)
      {
        return check_dist(function, alpha, beta, result) 
          && check_prob(function, p, result);
      } // bool check_dist_and_prob

      template <class RealType>
      inline bool check_mean(const char* function, const RealType& mean, RealType* result)
      {
        if(!(boost::math::isfinite)(mean) || (mean < 0))
        {
          *result = tools::domain_error<RealType>(
            function,
            "mean argument is %1%, but must be >= 0 !", mean);
          return false;
        }
        return true;
      } // bool check_mean
      template <class RealType>
      inline bool check_variance(const char* function, const RealType& variance, RealType* result)
      {
        if(!(boost::math::isfinite)(variance) || (variance < 0))
        {
          *result = tools::domain_error<RealType>(
            function,
            "variance argument is %1%, but must be >= 0 !", variance);
          return false;
        }
        return true;
      } // bool check_variance
  
    }

    // typedef beta_distribution<double> beta;
    // is deliberately NOT included to avoid a name clash with beta function.
    // Use beta_distribution<> mybeta(...) to construct type double. 

    template <class RealType = double>
    class beta_distribution
    {
    public:
      typedef RealType value_type;

      beta_distribution(RealType alpha, RealType beta) : m_alpha(alpha), m_beta(beta)
      {
        RealType result;
        beta_detail::check_dist(
          BOOST_CURRENT_FUNCTION,
          m_alpha,
          m_beta,
          &result);
      } // beta_distribution constructor.
      // Accessor functions:
      RealType alpha() const
      { 
        return m_alpha;
      }
      RealType beta() const
      { // .
        return m_beta;
      }

      // Estimation of the alpha & beta parameters.
      // http://en.wikipedia.org/wiki/Beta_distribution
      // gives formulae in section on parameter estimation.
      static RealType estimate_alpha(
        RealType mean, // Expected value of mean for probability x.
        RealType variance,// Expected value of variance for probability x.
        RealType probability) // x
      {
        RealType result; // of error checks.
        if(false == 
          beta_detail::check_prob(
          BOOST_CURRENT_FUNCTION, probability, &result)
          &&
          beta_detail::check_mean(
          BOOST_CURRENT_FUNCTION, mean, &result)
          &&
          beta_detail::check_variance(
          BOOST_CURRENT_FUNCTION, variance, &result)
          )
        {
          return result;
        }
        return mean * (( (mean * (1 - mean)) / variance)- 1);
      } // RealType estimate_alpha

      static RealType estimate_beta(
        RealType mean, // Expected value of mean for probability x.
        RealType variance,// Expected value of variance for probability x.
        RealType probability) // x
      {
        RealType result; // of error checks.
        if(false ==
          beta_detail::check_prob(
          BOOST_CURRENT_FUNCTION, probability, &result)
          &&
          beta_detail::check_mean(
          BOOST_CURRENT_FUNCTION, mean, &result)
          &&
          beta_detail::check_variance(
          BOOST_CURRENT_FUNCTION, variance, &result)
          )
        {
          return result;
        }
        return (1 - mean) * (((mean * (1 - mean)) /variance)-1);
      } //  RealType estimate_beta

    private:
      RealType m_alpha; // Two parameters of the beta distribution.
      RealType m_beta;
    }; // template <class RealType> class beta_distribution


    template <class RealType>
    inline RealType mean(const beta_distribution<RealType>& dist)
    { // Mean of beta distribution = np.
      return  dist.alpha() / (dist.alpha() + dist.beta());
    } // mean

    template <class RealType>
    inline RealType variance(const beta_distribution<RealType>& dist)
    { // Variance of beta distribution = np(1-p).
      RealType a = dist.alpha();
      RealType b = dist.beta();
      return  (a * b) / ((a + b ) * (a + b) * (a + b + 1));
    } // variance

    template <class RealType>
    RealType pdf(const beta_distribution<RealType>& dist, const RealType x)
    { // Probability Density/Mass Function.
      BOOST_FPU_EXCEPTION_GUARD

      using boost::math::tools::domain_error;
      using namespace std; // for ADL of std functions

      RealType a = dist.alpha();
      RealType b = dist.beta();

      // Argument checks:
      RealType result;
      if(false == beta_detail::check_dist_and_x(
        BOOST_CURRENT_FUNCTION,
        a, b, x,
        &result))
      {
        return result;
      }
      using boost::math::beta;
      return pow(x, (a-1)) * pow((1 - x), (b-1))/ beta(a, b);
    } // pdf

    template <class RealType>
    RealType cdf(const beta_distribution<RealType>& dist, const RealType x)
    { // Cumulative Distribution Function beta.
      using boost::math::tools::domain_error;
      using namespace std; // for ADL of std functions

      RealType a = dist.alpha();
      RealType b = dist.beta();

      // Argument checks:
      RealType result;
      if(false == beta_detail::check_dist_and_x(
        BOOST_CURRENT_FUNCTION,
        a, b, x,
        &result))
      {
        return result;
      }
      return ibeta(a, b, x);
    } // beta cdf

    template <class RealType>
    RealType cdf(const complemented2_type<beta_distribution<RealType>, RealType>& c)
    { // Complemented Cumulative Distribution Function beta.

      using boost::math::tools::domain_error;
      using namespace std; // for ADL of std functions

      RealType const& x = c.param;
      beta_distribution<RealType> const& dist = c.dist;
      RealType a = dist.alpha();
      RealType b = dist.beta();

      // Argument checks:
      RealType result;
      if(false == beta_detail::check_dist_and_x(
        BOOST_CURRENT_FUNCTION,
        a, b, x,
        &result))
      {
        return result;
      }
      // Calculate cdf beta using the incomplete beta function.
      // Use of ibeta here prevents cancellation errors in calculating
      // 1 - x if x is very small, perhaps smaller than machine epsilon.
      return 1 - ibeta(a, b, x);
    } // beta cdf

    namespace detail
    {
      template <class RealType>
      struct beta_functor
      {
        beta_functor(const beta_distribution<RealType>& d, const RealType& target, bool c = false)
          : dist(d), t(target), complement(c) {}

        RealType operator()(const RealType x)
        {
          if(x >= dist.alpha())
            return 1; // Any positive value will do.
          return complement ? t - cdf(boost::math::complement(dist, x)) : cdf(dist, x) - t;
        }
      private:
        const beta_distribution<RealType>& dist;
        RealType t;
        bool complement;
      }; // struct beta_functor
    } // namespace detail

    template <class RealType>
    RealType quantile(const beta_distribution<RealType>& dist, const RealType& p)
    { // Quantile or Percent Point beta function or .
      // Inverse cumulative probability distribution.
      // Return x (0 <= x <= 1),
      // for a given probability p (0 <= p <= 1).
      // These functions take a probability as an argument
      // and return a value such that the probability that a random variable x
      // will be less than or equal to that value
      // is whatever probability you supplied as an argument. 

      // Argument checks:
      RealType result;
      RealType a = dist.alpha();
      RealType b = dist.beta();
      if(false == beta_detail::check_dist_and_prob(
        BOOST_CURRENT_FUNCTION,
        a, b, p,
        &result))
      {
        return result;
      }
      // Special cases:
      if (p == 0)
      {
        return 0;
      }
      if (p == 1)
      {
        return 1;
      }
      // Solve for quantile numerically:
      detail::beta_functor<RealType> f(dist, p);
      tools::eps_tolerance<RealType> tol(tools::digits<RealType>());
      boost::uintmax_t max_iter = 1000;
      std::pair<RealType, RealType> r = tools::bracket_and_solve_root(
        f,
        dist.alpha() / 2, // (dist.alpha + dist.beta()) / 4
        static_cast<RealType>(8),
        true,
        tol,
        max_iter);
      // Special cases?
      if(max_iter >= 1000)
      {
        tools::logic_error<RealType>(BOOST_CURRENT_FUNCTION, "Unable to locate the root within a reasonable number of iterations, closest approximation so far was %1%", r.first);
      }
      // return centre point of range found:
      return r.first + (r.second - r.first) / 2;
    } // quantile

    template <class RealType>
    RealType quantile(const complemented2_type<beta_distribution<RealType>, RealType>& c)
    { // Complement Quantile or Percent Point beta function .
      // Return the number of expected x for a given
      // complement of the probability q.
      //
      // Error checks:
      RealType q = c.param;
      const beta_distribution<RealType>& dist = c.dist;
      RealType result;
      if(false == beta_detail::check_dist_and_prob(
        BOOST_CURRENT_FUNCTION,
        dist.alpha(),
        dist.beta(),
        q,
        &result))
      {
        return result;
      }
      // Special cases:
      if(q == 1)
      {  
        return 0;
      }
      if(q == 0)
      { 
        return 1;
      }

      // Solve for quantile numerically:
      detail::beta_functor<RealType> f(dist, q, true);
      tools::eps_tolerance<RealType> tol(tools::digits<RealType>());
      boost::uintmax_t max_iter = 1000;
      std::pair<RealType, RealType> r = tools::bracket_and_solve_root(
        f,
        dist.alpha() / 2, // (dist.alpha + dist.beta()) / 4
        static_cast<RealType>(8),
        true,
        tol,
        max_iter);
      if(max_iter >= 1000)
      {
        tools::logic_error<RealType>(BOOST_CURRENT_FUNCTION, "Unable to locate the root within a reasonable number of iterations, closest approximation so far was %1%", r.first);
      }
      // return centre point of range found:
      return r.first + (r.second - r.first) / 2;
    } // Quantile Complement

    template <class RealType>
    inline RealType mode(const beta_distribution<RealType>& dist)
    {
      RealType result;
      if ((dist.alpha() <= 1))
      {
        result = tools::domain_error<RealType>(
          BOOST_CURRENT_FUNCTION,
          "mode undefined for alpha = %1%, must be > 1!", dist.alpha());
        return result;
      }

      if ((dist.beta() <= 1))
      {
        result = tools::domain_error<RealType>(
          BOOST_CURRENT_FUNCTION,
          "mode undefined for beta = %1%, must be > 1!", dist.beta());
        return result;
      }
      RealType a = dist.alpha();
      RealType b = dist.beta();
      return (a-1) / (a + b - 2);
    } // mode

    template <class RealType>
    inline RealType skewness(const beta_distribution<RealType>& dist)
    {
      using namespace std; // ADL of std functions.
      RealType a = dist.alpha();
      RealType b = dist.beta();
      return (2 * (b-a) * sqrt(a + B + 1)) / ((a + b + 2) * sqrt(a * b));
    } // skewness

    template <class RealType>
    inline RealType kurtosis(const beta_distribution<RealType>& dist)
    {
      RealType a = dist.alpha();
      RealType b = dist.beta();
      return 3 + 6 * a * a * a - a * a * (2 * b -1) + b * b * (b + 1) - 2 * a * b * (b + 2) /
        a * b * (a + b + 2) * (a + b + 3);
    } // kurtosis

    template <class RealType>
    inline RealType kurtosis_excess(const beta_distribution<RealType>& dist)
    {
      RealType a = dist.alpha();
      RealType b = dist.beta();
      return 6 * a * a * a - a * a * (2 * b -1) + b * b * (b + 1) - 2 * a * b * (b + 2) /
        a * b * (a + b + 2) * (a + b + 3);
    } // kurtosis_excess

  } // namespace math
} // namespace boost

// This include must be at the end, *after* the accessors
// for this distribution have been defined, in order to
// keep compilers that support two-phase lookup happy.
#include <boost/math/distributions/detail/derived_accessors.hpp>

#endif // BOOST_MATH_DIST_BETA_HPP

