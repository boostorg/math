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
// It is a multivariate generalization of the dirichlet distribution, hence its
// alternative name of multivariate dirichlet distribution (MBD).
// Dirichlet distributions are commonly used as prior distributions in
// Bayesian statistics, and in fact the Dirichlet distribution is the
// conjugate prior of the categorical distribution and multinomial distribution.

#ifndef BOOST_MATH_DIST_DIRICHLET_HPP
#define BOOST_MATH_DIST_DIRICHLET_HPP

#include <boost/math/distributions/fwd.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/special_functions/digamma.hpp>
#include <boost/math/distributions/complement.hpp>
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
template <class VectorType, class RealType, class Policy>
inline bool check_concentration(const char *function,
                                const VectorType &concentration,
                                RealType *result,
                                const Policy &pol)
{
  for (size_t i = 0; i < concentration.size(); ++i)
  {
    if (!(boost::math::isfinite)(i) || (i <= 0))
    {
      *result = policies::raise_domain_error<RealType>(
          function,
          "Concentration Parameter is %1%, but must be > 0 !", concentration, pol);
      return false;
    }
  }
  return true;
} // bool check_concentration

template <class VectorType, class RealType, class Policy>
inline bool check_prob(const char *function,
                       const VectorType &p,
                       RealType *result,
                       const Policy &pol)
{
  for (size_t i = 0; i < p.size(); ++i)
  {
    if ((i < 0) || (i > 1) || !(boost::math::isfinite)(i))
    {
      *result = policies::raise_domain_error<RealType>(
          function,
          "Probability argument is %1%, but must be >= 0 and <= 1 !", i, pol);
      return false;
    }
  }
  return true;
} // bool check_prob

template <class VectorType, class RealType, class Policy>
inline bool check_x(const char *function,
                    const VectorType &x,
                    RealType *result,
                    const Policy &pol)
{
  for (size_t i = 0; i < x.size(); ++i)
  {
    if (!(boost::math::isfinite)(x) || (x < 0) || (x > 1))
    {
      *result = policies::raise_domain_error<RealType>(
          function,
          "x argument is %1%, but must be >= 0 and <= 1 !", x, pol);
      return false;
    }
  }
  return true;
} // bool check_x

template <class VectorType, class RealType, class Policy>
inline bool check_dist(const char *function,
                       const VectorType &concentration,
                       RealType *result,
                       const Policy &pol)
{
  return check_concentration(function, concentration, result, pol);
} // bool check_dist

template <class VectorType, class RealType, class Policy>
inline bool check_dist_and_x(const char *function,
                             const VectorType &concentration,
                             const VectorType &x,
                             RealType *result,
                             const Policy &pol)
{
  return check_dist(function, concentration, result, pol) && check_x(function, x, result, pol);
} // bool check_dist_and_x

template <class VectorType, class RealType, class Policy>
inline bool check_dist_and_prob(const char *function,
                                const VectorType &concentration,
                                const VectorType &p,
                                RealType *result,
                                const Policy &pol)
{
  return check_dist(function, concentration, result, pol) && check_prob(function, p, result, pol);
} // bool check_dist_and_prob

template <class VectorType, class RealType, class Policy>
inline bool check_mean(const char *function,
                       const VectorType &mean,
                       RealType *result,
                       const Policy &pol)
{
  for (size_t i = 0; i < mean.size(); ++i)
  {
    if (!(boost::math::isfinite)(i) || (i <= 0))
    {
      *result = policies::raise_domain_error<RealType>(
          function,
          "mean argument is %1%, but must be > 0 !", i, pol);
      return false;
    }
  }
  return true;
} // bool check_mean

template <class VectorType, class RealType, class Policy>
inline bool check_variance(const char *function,
                           const VectorType &variance,
                           RealType *result,
                           const Policy &pol)
{
  for (size_t i = 0; i < variance.size(); ++i)
  {
    if (!(boost::math::isfinite)(i) || (i <= 0))
    {
      *result = policies::raise_domain_error<RealType>(
          function,
          "variance argument is %1%, but must be > 0 !", i, pol);
      return false;
    }
  }
  return true;
} // bool check_variance
} // namespace dirichlet_detail

template <class VectorType = std::vector<double>, class RealType = double, class Policy = policies::policy<>>
class dirichlet_distribution
{
public:
  dirichlet_distribution(VectorType concentration) : concentration(concentration)
  {
    RealType result;
    dirichlet_detail::check_dist(
        "boost::math::dirichlet_distribution<%1%>::dirichlet_distribution",
        concentration,
        &result, Policy());
    sum_concentration = accumulate(concentration.begin(), concentration.end(), 0);
  } // dirichlet_distribution constructor.

  // Accessor functions:
  VectorType Concentration() const
  {
    return concentration;
  }

  size_t Order() const
  {
    return concentration.size();
  }

  static VectorType find_concentration(
      VectorType mean,     // Expected value of mean.
      VectorType variance) // Expected value of variance.
  {
    assert(("Dimensions of mean and variance must be same!", mean.size() == variance.size()));
    static const char *function = "boost::math::dirichlet_distribution<%1%>::find_concentration";
    RealType result = 0; // of error checks.
    if (!(dirichlet_detail::check_mean(function, mean, &result, Policy()) && dirichlet_detail::check_variance(function, variance, &result, Policy())))
    {
      return result;
    }
    VectorType c;
    for (size_t i = 0; i < mean.size(); ++i)
    {
      c.push_back(mean[i] * (((mean[i] * (1 - mean[i])) / variance[i]) - 1));
    }
    return c;
  } // RealType find_concentration

  // TODO
  // static VectorType find_concentration(
  //     VectorType x,           //  x.
  //     VectorType probability) // cdf
  // {
  //   assert(("", x.size() == probability.size()));
  //   static const char *function = "boost::math::dirichlet_distribution<%1%>::find_conentration";
  //   RealType result = 0; // of error checks.
  //   if (!(dirichlet_detail::check_prob(function, probability, &result, Policy()) && dirichlet_detail::check_x(function, x, &result, Policy())))
  //   {
  //     return result;
  //   }
  //   return ;
  // } // RealType find_concentration(x, probability)

private:
  VectorType concentration; // https://en.wikipedia.org/wiki/Concentration_parameter.
  RealType sum_concentration;
}; // template <class RealType, class Policy> class dirichlet_distribution


template <class VectorType, class RealType, class Policy>
inline const std::pair<RealType, RealType> range(const dirichlet_distribution<VectorType, RealType, Policy> & /* dist */)
{ // Range of permissible values for random variable x.
  using boost::math::tools::max_value;
  return std::pair<RealType, RealType>(static_cast<RealType>(0), static_cast<RealType>(1));
}


template <class VectorType, class RealType, class Policy>
inline const std::pair<RealType, RealType> support(const dirichlet_distribution<VectorType, RealType, Policy> & /* dist */)
{ // Range of supported values for random variable x.
  // This is range where cdf rises from 0 to 1, and outside it, the pdf is zero.
  return std::pair<RealType, RealType>(static_cast<RealType>(0), static_cast<RealType>(1));
}


template <class VectorType, class RealType, class Policy>
inline VectorType mean(const dirichlet_distribution<VectorType, RealType, Policy> &dist)
{ // Mean of dirichlet distribution = c[i]/sum(c).
  VectorType m;
  for (size_t i = 0; i < dist.Order(); ++i)
  {
    m.push_back(dist.concentration[i] / dist.sum_concentration);
  }
  return m;

} // mean


template <class VectorType, class RealType, class Policy>
inline VectorType variance(const dirichlet_distribution<VectorType, RealType, Policy> &dist)
{
  VectorType v;
  for (size_t i = 0; i < dist.Order(); ++i)
  {
    v.push_back(dist.concentration[i] / dist.sum_concentration * (1 - dist.concentration[i] / dist.sum_concentration) / (1 + dist.sum_concentration));
  }
  return v;
} // variance


template <class VectorType, class RealType, class Policy>
inline VectorType mode(const dirichlet_distribution<VectorType, RealType, Policy> &dist)
{
  static const char *function = "boost::math::mode(dirichlet_distribution<%1%> const&)";
  VectorType m;
  for (size_t i = 0; i < dist.Order(); ++i)
  {
    if ((dist.concentration[i] <= 1))
    {
      result = policies::raise_domain_error<RealType>(
          function,
          "mode undefined for alpha = %1%, must be > 1!", dist.alpha(), Policy());
      return result;
    }
    else
    {
      m.push_back((dist.concentration[i] - 1) / (dist.sum_concentration - dist.Order()));
    }
  }
  return m;
} // mode


template <class VectorType, class RealType, class Policy>
inline RealType entropy(const dirichlet_distribution<VectorType, RealType, Policy> &dist)
{
  RealType t1 = 1;
  for (size_t i = 0; i < dist.Order(); ++i)
  {
    t1 *= tgamma(dist.concentration[i]);
  }
  t1 = std::log(t1 / tgamma(dist.sum_concentration));
  RealType t2 = (dist.sum_concentration - dist.Order()) * digamma(dist.sum_concentration);
  RealType t3 = 0;
  for (size_t i = 0; i < dist.Order(); ++i)
  {
    t3 += (dist.concentration[i] - 1) * digamma(dist.concentration[i]);
  }
  return t1 + t2 - t3;
}


template <class VectorType, class RealType, class Policy>
inline RealType pdf(const dirichlet_distribution<VectorType, RealType, Policy> &dist, const VectorType &x)
{ // Probability Density/Mass Function.
  BOOST_FPU_EXCEPTION_GUARD

  static const char *function = "boost::math::pdf(dirichlet_distribution<%1%> const&, %1%)";

  BOOST_MATH_STD_USING // for ADL of std functions

  // Argument checks:
  RealType result = 0;
  if (!dirichlet_detail::check_dist_and_x(function, x, &result, Policy()))
  {
    return result;
  }
  using boost::math::tgamma;
  RealType f = 1;
  for (size_t i = 0; i < dist.Order(); ++i)
  {
    f *= std::pow(x[i], dist.concentration[i] - 1);
  }
  f /= dist.normalizing_factor;
  return f;
} // pdf


template <class VectorType, class RealType, class Policy>
inline RealType cdf(const dirichlet_distribution<VectorType, RealType, Policy> &dist, const VectorType &x)
{ // Cumulative Distribution Function dirichlet.
  BOOST_MATH_STD_USING // for ADL of std functions

      static const char *function = "boost::math::cdf(dirichlet_distribution<%1%> const&, %1%)";

  // Argument checks:
  RealType result = 0;
  if (!dirichlet_detail::check_dist_and_x(function, dist.concentration, x, &result, Policy()))
  {
    return result;
  }
  RealType c = 1;
  for (size_t i = 0; i < dist.Order(); ++i)
  {
    c *= std::pow(x[i], dist.concentration[i]) / tgamma(dist.concentration[i]) / dist.concentration[i];
  }
  c *= tgamma(dist.sum_concentration);
  return c;
} // dirichlet cdf

template <class VectorType, class RealType, class Policy>
inline RealType cdf(const complemented2_type<dirichlet_distribution<VectorType, RealType, Policy>, RealType> &c)
{ // Complemented Cumulative Distribution Function dirichlet.

  BOOST_MATH_STD_USING // for ADL of std functions

  static const char *function = "boost::math::cdf(dirichlet_distribution<%1%> const&, %1%)";

  RealType const &x = c.param;
  dirichlet_distribution<RealType, Policy> const &dist = c.dist;

  // Argument checks:
  RealType result = 0;
  if (!dirichlet_detail::check_dist_and_x(function, x, &result, Policy()))
  {
    return result;
  }
  RealType cumm = 1;
  for (size_t i = 0; i < dist.Order(); ++i)
  {
    cumm *= std::pow(x[i], dist.concentration[i]) / tgamma(dist.concentration[i]) / dist.concentration[i];
  }
  cumm *= tgamma(dist.sum_concentration);
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

#endif // BOOST_MATH_DIST_dirichlet_HPP
