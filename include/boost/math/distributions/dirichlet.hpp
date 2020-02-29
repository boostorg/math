// boost/math/distributions/dirichlet.hpp

// Copyright Mrityunjay Tripathi 2020.

// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

// https://en.wikipedia.org/wiki/Dirichlet_distribution
// https://mast.queensu.ca/~communications/Papers/msc-jiayu-lin.pdf

// The Dirichlet distribution is a family of continuous multivariate probability
// distributions parameterized by a vector 'alpha' of positive reals.
// It is a multivariate generalization of the beta distribution, hence its
// alternative name of multivariate beta distribution (MBD).
// Dirichlet distributions are commonly used as prior distributions in
// Bayesian statistics, and in fact the Dirichlet distribution is the
// conjugate prior of the categorical distribution and multinomial distribution.

#ifndef BOOST_MATH_DISTRIBUTIONS_DIRICHLET_HPP
#define BOOST_MATH_DISTRIBUTIONS_DIRICHLET_HPP

#include <boost/math/distributions/fwd.hpp>
#include <boost/multiprecision/float128.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/distributions/complement.hpp>
#include <boost/math/special_functions/digamma.hpp>
#include <boost/math/special_functions/fpclassify.hpp>
#include <boost/math/distributions/detail/common_error_handling.hpp>

#if defined(BOOST_MSVC)
#pragma warning(push)
#pragma warning(disable : 4702) // unreachable code
// in domain_error_imp in error_handling
#endif

#include <utility>

namespace boost
{
namespace math
{
namespace dirichlet_detail
{
// Common error checking routines for dirichlet distribution function:
template <class RandomAccessContainer, class Policy>
inline bool check_alpha(const char *function,
                        const RandomAccessContainer &alpha,
                        typename RandomAccessContainer::value_type *result,
                        const Policy &pol)
{
  using RealType = RandomAccessContainer::value_type;
  for (decltype(alpha.size()) i = 0; i < alpha.size(); ++i)
  {
    if (!(boost::math::isfinite)(alpha[i]) || (alpha[i] <= 0))
    {
      *result = policies::raise_domain_error<RealType>(
          function,
          "alpha Parameter is %1%, but must be > 0 !", alpha[i], pol);
      return false;
    }
  }
  return true;
} // bool check_alpha

template <class RandomAccessContainer, class Policy>
inline bool check_prob(const char *function,
                       const RandomAccessContainer &p,
                       typename RandomAccessContainer::value_type *result,
                       const Policy &pol)
{
  using RealType = RandomAccessContainer::value_type;
  for (decltype(p.size()) i = 0; i < p.size(); ++i)
  {
    if ((p[i] < 0) || (p[i] > 1) || !(boost::math::isfinite)(p[i]))
    {
      *result = policies::raise_domain_error<RealType>(
          function,
          "Probability argument is %1%, but must be >= 0 and <= 1 !", p[i], pol);
      return false;
    }
  }
  return true;
} // bool check_prob

template <class RandomAccessContainer, class Policy>
inline bool check_x(const char *function,
                    const RandomAccessContainer &x,
                    typename RandomAccessContainer::value_type *result,
                    const Policy &pol)
{
  using RealType = typename RandomAccessContainer::value_type;
  for (decltype(x.size()) i = 0; i < x.size(); ++i)
  {
    if (!(boost::math::isfinite)(x[i]) || (x[i] < 0) || (x[i] > 1))
    {
      *result = policies::raise_domain_error<RealType>(
          function,
          "x argument is %1%, but must be >= 0 and <= 1 !", x[i], pol);
      return false;
    }
  }
  return true;
} // bool check_x

template <class RandomAccessContainer, class Policy>
inline bool check_dist(const char *function,
                       const RandomAccessContainer &alpha,
                       typename RandomAccessContainer::value_type *result,
                       const Policy &pol)
{
  return check_alpha(function, alpha, result, pol);
} // bool check_dist

template <class RandomAccessContainer, class Policy>
inline bool check_dist_and_x(const char *function,
                             const RandomAccessContainer &alpha,
                             const RandomAccessContainer &x,
                             typename RandomAccessContainer::value_type *result,
                             const Policy &pol)
{
  return check_dist(function, alpha, result, pol) && check_x(function, x, result, pol);
} // bool check_dist_and_x

template <class RandomAccessContainer, class Policy>
inline bool check_dist_and_prob(const char *function,
                                const RandomAccessContainer &alpha,
                                const RandomAccessContainer &p,
                                typename RandomAccessContainer::value_type *result,
                                const Policy &pol)
{
  return check_dist(function, alpha, result, pol) && check_prob(function, p, result, pol);
} // bool check_dist_and_prob

template <class RandomAccessContainer, class Policy>
inline bool check_mean(const char *function,
                       const RandomAccessContainer &mean,
                       typename RandomAccessContainer::value_type *result,
                       const Policy &pol)
{
  using RealType = typename RandomAccessContainer::value_type;
  for (decltype(mean.size()) i = 0; i < mean.size(); ++i)
  {
    if (!(boost::math::isfinite)(mean[i]) || (mean[i] <= 0))
    {
      *result = policies::raise_domain_error<RealType>(
          function,
          "mean argument is %1%, but must be > 0 !", mean[i], pol);
      return false;
    }
  }
  return true;
} // bool check_mean

template <class RandomAccessContainer, class Policy>
inline bool check_variance(const char *function,
                           const RandomAccessContainer &variance,
                           typename RandomAccessContainer::value_type *result,
                           const Policy &pol)
{
  using RealType = typename RandomAccessContainer::value_type;
  for (decltype(variance.size()) i = 0; i < variance.size(); ++i)
  {
    if (!(boost::math::isfinite)(variance[i]) || (variance[i] <= 0))
    {
      *result = policies::raise_domain_error<RealType>(
          function,
          "variance argument is %1%, but must be > 0 !", variance[i], pol);
      return false;
    }
  }
  return true;
} // bool check_variance

template <class RandomAccessContainer, class Policy>
inline bool check_mean_and_variance(const char *function,
                                    const RandomAccessContainer &mean,
                                    const RandomAccessContainer &variance,
                                    typename RandomAccessContainer::value_type *result,
                                    const Policy &pol)
{
  return check_mean(function, mean, result, pol) && check_variance(function, variance, result, pol);
} // bool check_mean_and_variance

template <class RandomAccessContainer>
inline typename RandomAccessContainer::value_type mvar_beta(
    const RandomAccessContainer &alpha,
    const typename RandomAccessContainer::value_type &b)
{
  // B(a1,a2,...ak) = tgamma(a1+a2+...+ak)/(tgamma(a1)*tgamma(a2)...*tgamma(ak)
  using RealType = typename RandomAccessContainer::value_type;
  RealType mb;
  RealType alpha_sum = accumulate(alpha.begin(), alpha.end(), b * alpha.size());
  for (decltype(alpha.size()) i = 0; i < alpha.size(); ++i)
  {
    mb *= tgamma(alpha[i] + b);
  }
  mb /= tgamma(alpha_sum);
  return mb;
} // mvar_beta

template <class RandomAccessContainer>
inline typename RandomAccessContainer::value_type alpha0(const RandomAccessContainer &alpha)
{
  return accumulate(alpha.begin(), alpha.end(), 0);
}
} // namespace dirichlet_detail

template <class RandomAccessContainer = std::vector<double>, class Policy = policies::policy<>>
class dirichlet_distribution
{
  using RealType = typename RandomAccessContainer::value_type;

public:
  dirichlet_distribution(RandomAccessContainer &&alpha) : m_alpha(alpha)
  {
    RealType result;
    dirichlet_detail::check_dist(
        "boost::math::dirichlet_distribution<%1%>::dirichlet_distribution",
        alpha,
        &result, Policy());
  } // dirichlet_distribution constructor.

  // Get the concentration parameters.
  RandomAccessContainer &alpha() const { return m_alpha; }
  // Set the concentration parameters.
  RandomAccessContainer &alpha() { return m_alpha; }

  // Get the order of concentration parameters.
  decltype(m_alpha.size()) &order() const { return m_alpha.size(); }

  // Get alpha from mean and variance.
  static void find_alpha(
      RandomAccessContainer &&mean,     // Expected value of mean.
      RandomAccessContainer &&variance) // Expected value of variance.
  {
    assert(("Dimensions of mean and variance must be same!", mean.size() == variance.size()));
    static const char *function = "boost::math::dirichlet_distribution<%1%>::alpha_from_mean_and_variance";
    RealType result = 0; // of error checks.
    if (!(dirichlet_detail::check_mean(function, mean, &result, Policy()) && dirichlet_detail::check_variance(function, variance, &result, Policy())))
    {
      return result;
    }
    for (decltype(mean.size()) i = 0; i < mean.size(); ++i)
    {
      m_alpha[i] = mean[i] * (((mean[i] * (1 - mean[i])) / variance[i]) - 1);
    }
  } // void find_alpha

private:
  RandomAccessContainer m_alpha; // https://en.wikipedia.org/wiki/Concentration_parameter.
}; // template <class RealType, class Policy> class dirichlet_distribution

template <class RandomAccessContainer, class Policy>
inline const std::pair<
    typename RandomAccessContainer::value_type,
    typename RandomAccessContainer::value_type>
range(const dirichlet_distribution<RandomAccessContainer, Policy> & /* dist */)
{ // Range of permissible values for random variable x.
  using boost::math::tools::max_value;
  using RealType = typename RandomAccessContainer::value_type;
  return std::pair<RealType, RealType>(static_cast<RealType>(0), static_cast<RealType>(1));
}

template <class RandomAccessContainer, class Policy>
inline const std::pair<
    typename RandomAccessContainer::value_type,
    typename RandomAccessContainer::value_type>
support(const dirichlet_distribution<RandomAccessContainer, Policy> & /* dist */)
{ // Range of supported values for random variable x.
  // This is range where cdf rises from 0 to 1, and outside it, the pdf is zero.
  using RealType = typename RandomAccessContainer::value_type;
  return std::pair<RealType, RealType>(static_cast<RealType>(0), static_cast<RealType>(1));
}

template <class RandomAccessContainer, class Policy>
inline RandomAccessContainer mean(const dirichlet_distribution<RandomAccessContainer, Policy> &dist)
{ // Mean of dirichlet distribution = c[i]/sum(c).
  // using RealType = typename RandomAccessContainer::value_type;
  RandomAccessContainer m(dist.order());
  for (decltype(dist.order()) i = 0; i < dist.order(); ++i)
  {
    m[i] = dist.m_alpha[i] / dirichlet_detail::alpha0(dist.m_alpha);
  }
  return m;
} // mean

template <class RandomAccessContainer, class Policy>
inline RandomAccessContainer variance(const dirichlet_distribution<RandomAccessContainer, Policy> &dist)
{
  RandomAccessContainer v(dist.order());
  for (decltype(dist.order()) i = 0; i < dist.order(); ++i)
  {
    v[i] = dist.m_alpha[i] / dirichlet_detail::alpha0(dist.m_alpha) * (1 - dist.m_alpha[i] / dirichlet_detail::alpha0(dist.alpha)) / (1 + dirichlet_detail::alpha0(dist.alpha));
  }
  return v;
} // variance

template <class RandomAccessContainer, class Policy>
inline RandomAccessContainer standard_deviation(const dirichlet_distribution<RandomAccessContainer, Policy> &dist)
{
  RandomAccessContainer std = variance(dist);
  for (decltype(dist.order()) i = 0; i < s.size(); ++i)
  {
    std[i] = std::sqrt(std[i]);
  }
  return std;
} // standard_deviation

template <class RandomAccessContainer, class Policy>
inline RandomAccessContainer mode(const dirichlet_distribution<RandomAccessContainer, Policy> &dist)
{
  using RealType = typename RandomAccessContainer::value_type;
  static const char *function = "boost::math::mode(dirichlet_distribution<%1%> const&)";
  RandomAccessContainer m(dist.order());
  for (decltype(dist.order()) i = 0; i < m.size(); ++i)
  {
    if (dist.m_alpha[i] <= 1)
    {
      result = policies::raise_domain_error<RealType>(
          function,
          "mode undefined for alpha = %1%, must be > 1!", dist.m_alpha[i], Policy());
      return result;
    }
    else
    {
      m[i] = (dist.m_alpha[i] - 1) / (dirichlet_detail::alpha0(dist.m_alpha) - dist.order());
    }
  }
  return m;
} // mode

template <class RandomAccessContainer, class Policy>
inline typename RandomAccessContainer::value_type entropy(const dirichlet_distribution<RandomAccessContainer, Policy> &dist)
{
  using RealType = typename RandomAccessContainer::value_type;
  RealType ent = std::log(dirichlet_detail::mvar_beta(dist.m_alpha, 0)) + (dirichlet_detail::alpha0(dist.m_alpha) - dist.order()) * digamma(dirichlet_detail::alpha0(dist.m_alpha));
  for (decltype(dist.order()) i = 0; i < dist.order(); ++i)
  {
    ent += (dist.m_alpha[i] - 1) * digamma(dist.m_alpha[i]);
  }
  return ent;
}

template <class RandomAccessContainer, class Policy>
inline RandomAccessContainer skewness(const dirichlet_distribution<RandomAccessContainer, Policy> &dist)
{
  using RealType = typename RandomAccessContainer::value_type;
  RandomAccessContainer s(dist.order());
  RealType A = dirichlet_detail::alpha0(dist.m_alpha);
  RealType aj;
  for (decltype(dist.order()) i = 0; i < dist.order(); ++i)
  {
    aj = dist.m_alpha[i];
    s[i] = std::sqrt(aj * (A + 1) / (A - aj)) * ((aj + 2) * (aj + 1) * A * A / (aj * (A + 2) * (A - aj)) - 3 - aj * (A + 1) / (A - aj));
  }
  return s;
}

template <class RandomAccessContainer, class Policy>
inline RandomAccessContainer kurtosis(const dirichlet_distribution<RandomAccessContainer, Policy> &dist)
{
  using RealType = typename RandomAccessContainer::value_type;
  using std::pow;
  RandomAccessContainer k(dist.order());
  RealType A = dirichlet_detail::alpha0(dist.m_alpha);
  RealType aj;
  for (decltype(dist.order()) i = 0; i < dist.order(); ++i)
  {
    aj = dist.m_alpha[i];
    k[i] = ((aj + 2) * (aj + 1) * ((aj + 3) * A / (A + 3) / aj - 4) + 6 * (aj + 1) * aj / (A + 1) / A - 3 * pow(aj / A, 2)) / std::pow((A - aj) / A / (A + 1), 2);
  }
  return k;
}

template <class RandomAccessContainer, class Policy>
inline RandomAccessContainer kurtosis_excess(const dirichlet_distribution<RandomAccessContainer, Policy> &dist)
{
  RandomAccessContainer ke = kurtosis(dist);
  for (decltype(dist.order()) i = 0; i < dist.order(); ++i)
  {
    ke[i] = ke[i] - 3;
  }
  return ke;
}

template <class RandomAccessContainer, class Policy>
inline typename RandomAccessContainer::value_type pdf(
    const dirichlet_distribution<RandomAccessContainer, Policy> &dist,
    const RandomAccessContainer &x)
{ // Probability Density/Mass Function.
  using RealType = typename RandomAccessContainer::value_type;
  using std::pow;
  BOOST_FPU_EXCEPTION_GUARD
  static const char *function = "boost::math::pdf(dirichlet_distribution<%1%> const&, %1%)";
  BOOST_MATH_STD_USING // for ADL of std functions
  RealType result = 0;
  if (!dirichlet_detail::check_dist_and_x(function, x, &result, Policy()))
  {
    return result;
  }

  RealType f = 1;
  for (decltype(dist.order()) i = 0; i < dist.order(); ++i)
  {
    f *= pow(x[i], dist.m_alpha[i] - 1);
  }
  f /= dirichlet_detail::mvar_beta(dist.m_alpha, 0);
  return f;
} // pdf

template <class RandomAccessContainer, class Policy>
inline typename RandomAccessContainer::value_type cdf(
    const dirichlet_distribution<RandomAccessContainer, Policy> &dist,
    const RandomAccessContainer &x)
{ // Cumulative Distribution Function dirichlet.
  using RealType = typename RandomAccessContainer::value_type;
  using std::pow;
  BOOST_MATH_STD_USING // for ADL of std functions

  static const char *function = "boost::math::cdf(dirichlet_distribution<%1%> const&, %1%)";
  RealType result = 0; // Arguments check.
  if (!dirichlet_detail::check_dist_and_x(function, dist.alpha, x, &result, Policy()))
  {
    return result;
  }
  RealType c = 1;
  for (decltype(dist.order()) i = 0; i < dist.order(); ++i)
  {
    c *= pow(x[i], dist.m_alpha[i]) / tgamma(dist.m_alpha[i]) / dist.m_alpha[i];
  }
  c *= tgamma(dirichlet_detail::alpha0(dist.m_alpha));
  return c;
} // dirichlet cdf

template <class RandomAccessContainer, class Policy>
inline typename RandomAccessContainer::value_type cdf(
    const complemented2_type<dirichlet_distribution<RandomAccessContainer, Policy>,
    typename RandomAccessContainer::value_type> &c)
{ // Complemented Cumulative Distribution Function dirichlet.
  using RealType = typename RandomAccessContainer::value_type;
  using std::pow;
  BOOST_MATH_STD_USING // for ADL of std functions
  static const char *function = "boost::math::cdf(dirichlet_distribution<%1%> const&, %1%)";
  RealType const &x = c.param;
  dirichlet_distribution<RealType, Policy> const &dist = c.dist;
  RealType result = 0; // Argument checks.
  if (!dirichlet_detail::check_dist_and_x(function, x, &result, Policy()))
  {
    return result;
  }
  RealType cumm = 1;
  for (decltype(dist.order()) i = 0; i < dist.order(); ++i)
  {
    cumm *= pow(x[i], dist.m_alpha[i]) / tgamma(dist.m_alpha[i]) / dist.m_alpha[i];
  }
  cumm *= tgamma(dirichlet_detail::alpha0(dist.m_alpha));
  return cumm;
} // dirichlet cdf

} // namespace math
} // namespace boost

// This include must be at the end, *after* the accessors
// for this distribution have been defined, in order to
// keep compilers that support two-phase lookup happy.
#include <boost/math/distributions/detail/derived_accessors.hpp>

#if defined(BOOST_MSVC)
#pragma warning(pop)
#endif

#endif // BOOST_MATH_DIST_DIRICHLET_HPP
