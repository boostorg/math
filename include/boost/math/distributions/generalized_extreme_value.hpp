//  Copyright Nick Thompson 2026.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_DISTRIBUTIONS_GENERALIZED_EXTREME_VALUE_HPP
#define BOOST_MATH_DISTRIBUTIONS_GENERALIZED_EXTREME_VALUE_HPP

// The generalized extreme value distribution, with location mu, scale sigma and shape xi as in
// https://en.wikipedia.org/wiki/Generalized_extreme_value_distribution:
//
//    cdf(x) = exp(-t(x)),   t(x) = (1 + xi (x - mu)/sigma)^(-1/xi),   or exp(-(x - mu)/sigma) for xi = 0.
//
// xi > 0 is the Frechet type, bounded below at mu - sigma/xi; xi < 0 the reversed Weibull type, bounded
// above there; xi = 0 the Gumbel distribution, extreme_value_distribution. (SciPy's genextreme takes c = -xi.)
// log t = -log1p(xi z)/xi and the quantile's expm1(-xi log y)/xi tend to their Gumbel limits as xi -> 0
// without cancellation.

#include <boost/math/tools/config.hpp>
#include <boost/math/tools/precision.hpp>
#include <boost/math/tools/tuple.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/math/special_functions/expm1.hpp>
#include <boost/math/special_functions/fpclassify.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/special_functions/log1p.hpp>
#include <boost/math/special_functions/zeta.hpp>
#include <boost/math/distributions/complement.hpp>
#include <boost/math/distributions/extreme_value.hpp>
#include <boost/math/distributions/fwd.hpp>
#include <boost/math/distributions/detail/common_error_handling.hpp>
#include <boost/math/policies/error_handling.hpp>
#include <boost/math/policies/policy.hpp>

#ifndef BOOST_MATH_BUILD_MODULE
#include <cmath>
#include <limits>
#include <utility>
#endif

BOOST_MATH_NAMESPACE_BEGIN

BOOST_MATH_EXPORT template <class RealType = double, class Policy = policies::policy<> >
class generalized_extreme_value_distribution
{
public:
   using value_type = RealType;
   using policy_type = Policy;

   explicit generalized_extreme_value_distribution(RealType location = 0, RealType scale = 1, RealType shape = 0)
      : m_location(location), m_scale(scale), m_shape(shape)
   {
      constexpr auto function = "boost::math::generalized_extreme_value_distribution<%1%>::generalized_extreme_value_distribution";
      RealType result;
      detail::check_location(function, location, &result, Policy());
      detail::check_scale(function, scale, &result, Policy());
      detail::check_finite(function, shape, &result, Policy());
   }

   RealType location() const { return m_location; }
   RealType scale() const { return m_scale; }
   RealType shape() const { return m_shape; }

private:
   RealType m_location;
   RealType m_scale;
   RealType m_shape;
};

BOOST_MATH_EXPORT using generalized_extreme_value = generalized_extreme_value_distribution<double>;

#ifdef __cpp_deduction_guides
BOOST_MATH_EXPORT template <class RealType>
generalized_extreme_value_distribution(RealType)->generalized_extreme_value_distribution<typename BOOST_MATH_NAMESPACE::tools::promote_args<RealType>::type>;
BOOST_MATH_EXPORT template <class RealType>
generalized_extreme_value_distribution(RealType, RealType)->generalized_extreme_value_distribution<typename BOOST_MATH_NAMESPACE::tools::promote_args<RealType>::type>;
BOOST_MATH_EXPORT template <class RealType>
generalized_extreme_value_distribution(RealType, RealType, RealType)->generalized_extreme_value_distribution<typename BOOST_MATH_NAMESPACE::tools::promote_args<RealType>::type>;
#endif

namespace detail {

template <class RealType, class Policy>
inline bool check_gev(const char* function, const generalized_extreme_value_distribution<RealType, Policy>& dist, RealType* result)
{
   return check_location(function, dist.location(), result, Policy())
      && check_scale(function, dist.scale(), result, Policy())
      && check_finite(function, dist.shape(), result, Policy());
}

// Where x lies relative to the support, and inside it L = -log t(x) = log1p(xi z)/xi, so that cdf(x) = exp(-exp(-L)).
enum class gev_region { below, inside, endpoint, above };

// The support is decided against the end point mu - sigma/xi that support() reports, before x is standardized,
// whose rounding could otherwise carry a point outside onto it. Inside, 1 + xi z can still round to 0 or below
// next to the end point, where the end point's limit is the answer to working precision.
template <class RealType>
gev_region gev_log_t(RealType x, RealType mu, RealType sigma, RealType xi, RealType* L)
{
   BOOST_MATH_STD_USING
   RealType z = (x - mu) / sigma;
   if (xi == 0)
   {
      *L = z;
      return gev_region::inside;
   }
   RealType endpoint = mu - sigma / xi;
   if (x == endpoint)
   {
      return gev_region::endpoint;
   }
   if ((xi > 0) == (x < endpoint))
   {
      // Beyond the lower end point for xi > 0, the upper for xi < 0:
      return xi > 0 ? gev_region::below : gev_region::above;
   }
   RealType s = xi * z;
   if (s <= -1)
   {
      return gev_region::endpoint;
   }
   // For tiny xi z, xi z may have lost bits to underflow; log1p(s)/xi = z (1 - s/2 + ...):
   *L = abs(s) < tools::epsilon<RealType>() ? z * (1 - s / 2) : RealType(BOOST_MATH_NAMESPACE::log1p(s) / xi);
   return gev_region::inside;
}

// (y^(-xi) - 1)/xi, the standardized quantile at t = y; -log y for xi = 0. expm1(-xi log y) carries
// log's absolute error times xi, so where |xi log y| > 1 pow, which then has nothing to cancel, is better.
template <class RealType, class Policy>
RealType gev_standard_quantile(RealType y, RealType log_y, RealType xi, const Policy& pol)
{
   BOOST_MATH_STD_USING
   RealType a = -xi * log_y;
   if (xi == 0 || abs(a) < tools::epsilon<RealType>())
   {
      return -log_y * (1 + a / 2);
   }
   if (abs(a) > 1)
   {
      return (pow(y, -xi) - 1) / xi;
   }
   return BOOST_MATH_NAMESPACE::expm1(a, pol) / xi;
}

// expm1(x) - x, without the cancellation for small x.
template <class RealType, class Policy>
RealType gev_expm1mx(RealType x, const Policy& pol)
{
   BOOST_MATH_STD_USING
   if (abs(x) > RealType(0.5))
   {
      return BOOST_MATH_NAMESPACE::expm1(x, pol) - x;
   }
   RealType term = x * x / 2;
   RealType sum = term;
   for (unsigned n = 3; abs(term) > tools::epsilon<RealType>() * abs(sum); ++n)
   {
      term *= x / n;
      sum += term;
   }
   return sum;
}

// sum_k c[k] A_k for k = 2, 3, 4, where A_k = lgamma(1 - k xi) - k lgamma(1 - xi) is the log of
// E[W^k]/E[W]^k for W = t^(-xi), the moments of which are Gamma(1 - k xi). The combinations the
// moments need cancel at leading order(s) in xi, so for small xi this sums the Taylor series
// lgamma(1 - z) = euler z + sum_{j >= 2} zeta(j) z^j / j term by term, in which they cancel exactly:
//    sum_k c[k] A_k = sum_{j >= 2} zeta(j)/j sum_k c[k] (k^j - k) xi^j.
// The series converges for |k xi| < 1 with ratio |k xi|; it is used to |k xi| = 3/4 for the largest k,
// beyond which the lgamma differences lose at most a digit or so.
template <class RealType, class Policy>
RealType gev_log_moment_combination(const int (&c)[3], RealType xi, const Policy& pol)
{
   BOOST_MATH_STD_USING
   int k_max = c[2] != 0 ? 4 : (c[1] != 0 ? 3 : 2);
   if (k_max * abs(xi) > RealType(0.75))
   {
      RealType l1 = BOOST_MATH_NAMESPACE::lgamma(1 - xi, pol);
      RealType sum = 0;
      for (int k = 2; k <= 4; ++k)
      {
         if (c[k - 2] != 0)
         {
            sum += c[k - 2] * (BOOST_MATH_NAMESPACE::lgamma(1 - k * xi, pol) - k * l1);
         }
      }
      return sum;
   }
   RealType power[3] = { 2, 3, 4 };   // k^j
   RealType xi_j = xi;                // xi^j
   RealType sum = 0;
   const std::uintmax_t max_terms = policies::get_max_series_iterations<Policy>();
   for (std::uintmax_t j = 2; ; ++j)
   {
      RealType weight = 0;
      for (int k = 2; k <= 4; ++k)
      {
         power[k - 2] *= k;
         weight += c[k - 2] * (power[k - 2] - k);
      }
      xi_j *= xi;
      RealType term = BOOST_MATH_NAMESPACE::zeta(RealType(j), pol) / j * weight * xi_j;
      sum += term;
      // The weights of the leading orders are zero, so only stop on a nonzero term:
      if (weight != 0 && abs(term) <= tools::epsilon<RealType>() * abs(sum))
      {
         return sum;
      }
      if (j > max_terms)
      {
         return policies::raise_evaluation_error<RealType>("boost::math::generalized_extreme_value moments", "Series did not converge for shape %1%.", xi, pol);
      }
   }
}

// The moments of xi this small are their Gumbel limits to working precision.
template <class RealType>
bool gev_is_gumbel(RealType xi)
{
   BOOST_MATH_STD_USING
   return abs(xi) < tools::epsilon<RealType>();
}

// The mean for xi >= 1 and the variance for xi >= 1/2, which diverge to +infinity (max_value for a type without one).
template <class RealType>
RealType gev_infinity()
{
   return std::numeric_limits<RealType>::has_infinity ? std::numeric_limits<RealType>::infinity() : tools::max_value<RealType>();
}

template <class RealType, class Policy>
RealType gev_undefined_moment(const char* function, const char* message, RealType xi)
{
   return policies::raise_domain_error<RealType>(function, message, xi, Policy());
}

} // namespace detail

BOOST_MATH_EXPORT template <class RealType, class Policy>
inline std::pair<RealType, RealType> range(const generalized_extreme_value_distribution<RealType, Policy>& /*dist*/)
{ // Range of permissible values for random variable x; the pdf is zero outside the support.
   using BOOST_MATH_NAMESPACE::tools::max_value;
   return std::pair<RealType, RealType>(
      std::numeric_limits<RealType>::has_infinity ? -std::numeric_limits<RealType>::infinity() : -max_value<RealType>(),
      std::numeric_limits<RealType>::has_infinity ? std::numeric_limits<RealType>::infinity() : max_value<RealType>());
}

BOOST_MATH_EXPORT template <class RealType, class Policy>
inline std::pair<RealType, RealType> support(const generalized_extreme_value_distribution<RealType, Policy>& dist)
{ // Range where the cdf rises from 0 to 1: bounded at mu - sigma/xi below for xi > 0, above for xi < 0.
   using BOOST_MATH_NAMESPACE::tools::max_value;
   RealType xi = dist.shape();
   RealType endpoint = xi == 0 ? RealType(0) : RealType(dist.location() - dist.scale() / xi);
   if (xi > 0)
   {
      return std::pair<RealType, RealType>(endpoint, max_value<RealType>());
   }
   if (xi < 0)
   {
      return std::pair<RealType, RealType>(-max_value<RealType>(), endpoint);
   }
   return std::pair<RealType, RealType>(-max_value<RealType>(), max_value<RealType>());
}

BOOST_MATH_EXPORT template <class RealType, class Policy>
inline RealType logpdf(const generalized_extreme_value_distribution<RealType, Policy>& dist, const RealType& x)
{
   BOOST_MATH_STD_USING
   constexpr auto function = "boost::math::logpdf(const generalized_extreme_value_distribution<%1%>&, %1%)";
   RealType result = 0;
   if (!detail::check_gev(function, dist, &result) || !detail::check_x_not_NaN(function, x, &result, Policy()))
   {
      return result;
   }
   RealType xi = dist.shape();
   if ((BOOST_MATH_NAMESPACE::isinf)(x))
   {
      return -std::numeric_limits<RealType>::infinity();
   }
   RealType L;
   switch (detail::gev_log_t(x, dist.location(), dist.scale(), xi, &L))
   {
   case detail::gev_region::inside:
      break;
   case detail::gev_region::endpoint:
      // The bounded end, where t = 0 for xi < 0 and the density is t^(1 + xi)/sigma; t = infinity for xi > 0.
      if (xi == -1)
      {
         return -log(dist.scale());
      }
      if (xi < -1)
      {
         return policies::raise_overflow_error<RealType>(function, nullptr, Policy());
      }
      return -std::numeric_limits<RealType>::infinity();
   default:
      return -std::numeric_limits<RealType>::infinity();
   }
   // log pdf = -log sigma + (1 + xi) log t - t, with log t = -L:
   RealType t = exp(-L);
   if ((BOOST_MATH_NAMESPACE::isinf)(t))
   {
      return -std::numeric_limits<RealType>::infinity();
   }
   return -log(dist.scale()) - (1 + xi) * L - t;
}

BOOST_MATH_EXPORT template <class RealType, class Policy>
inline RealType pdf(const generalized_extreme_value_distribution<RealType, Policy>& dist, const RealType& x)
{
   BOOST_MATH_STD_USING
   constexpr auto function = "boost::math::pdf(const generalized_extreme_value_distribution<%1%>&, %1%)";
   RealType result = 0;
   if (!detail::check_gev(function, dist, &result) || !detail::check_x_not_NaN(function, x, &result, Policy()))
   {
      return result;
   }
   RealType xi = dist.shape();
   if ((BOOST_MATH_NAMESPACE::isinf)(x))
   {
      return 0;
   }
   RealType L;
   switch (detail::gev_log_t(x, dist.location(), dist.scale(), xi, &L))
   {
   case detail::gev_region::inside:
      break;
   case detail::gev_region::endpoint:
      if (xi == -1)
      {
         return 1 / dist.scale();
      }
      if (xi < -1)
      {
         return policies::raise_overflow_error<RealType>(function, nullptr, Policy());
      }
      return 0;
   default:
      return 0;
   }
   RealType t = exp(-L);
   if ((BOOST_MATH_NAMESPACE::isinf)(t))
   {
      return 0;
   }
   // t^(1 + xi) exp(-t) / sigma:
   return exp(-(1 + xi) * L - t) / dist.scale();
}

namespace detail {

// cdf, its complement, and their logs; region and t decide them all.
template <class RealType, class Policy>
RealType gev_cdf_imp(const generalized_extreme_value_distribution<RealType, Policy>& dist, const RealType& x, bool complement, bool take_log, const char* function)
{
   BOOST_MATH_STD_USING
   RealType result = 0;
   if (!check_gev(function, dist, &result) || !check_x_not_NaN(function, x, &result, Policy()))
   {
      return result;
   }
   const RealType zero = take_log ? RealType(-std::numeric_limits<RealType>::infinity()) : RealType(0);
   const RealType one = take_log ? RealType(0) : RealType(1);
   if ((BOOST_MATH_NAMESPACE::isinf)(x))
   {
      return (x < 0) != complement ? zero : one;
   }
   RealType xi = dist.shape();
   RealType L;
   switch (gev_log_t(x, dist.location(), dist.scale(), xi, &L))
   {
   case gev_region::inside:
      break;
   case gev_region::endpoint:
      return (xi > 0) != complement ? zero : one;
   case gev_region::below:
      return complement ? one : zero;
   default:
      return complement ? zero : one;
   }
   RealType t = exp(-L);
   if (!complement)
   {
      return take_log ? RealType(-t) : RealType(exp(-t));
   }
   RealType q = -BOOST_MATH_NAMESPACE::expm1(-t, Policy());
   if (!take_log)
   {
      return q;
   }
   // log(1 - exp(-t)): log t - t/2 = -L - t/2 once t is below epsilon (t may have underflowed), through q
   // while it is small, else through log1p:
   if (t < tools::epsilon<RealType>())
   {
      return -L - t / 2;
   }
   return t <constants::ln_two<RealType>() ? RealType(log(q)) : RealType(BOOST_MATH_NAMESPACE::log1p(-exp(-t), Policy()));
}

template <class RealType, class Policy>
RealType gev_quantile_imp(const generalized_extreme_value_distribution<RealType, Policy>& dist, const RealType& p, bool complement, const char* function)
{
   BOOST_MATH_STD_USING
   RealType result = 0;
   if (!check_gev(function, dist, &result) || !check_probability(function, p, &result, Policy()))
   {
      return result;
   }
   RealType xi = dist.shape();
   RealType mu = dist.location();
   RealType sigma = dist.scale();
   // P = cdf(x):
   bool lower_end = complement ? p == 1 : p == 0;
   bool upper_end = complement ? p == 0 : p == 1;
   if (lower_end)
   {
      return xi > 0 ? RealType(mu - sigma / xi) : RealType(-policies::raise_overflow_error<RealType>(function, nullptr, Policy()));
   }
   if (upper_end)
   {
      return xi < 0 ? RealType(mu - sigma / xi) : policies::raise_overflow_error<RealType>(function, nullptr, Policy());
   }
   // t = -log P:
   RealType y = complement ? RealType(-BOOST_MATH_NAMESPACE::log1p(-p, Policy())) : RealType(-log(p));
   return mu + sigma * gev_standard_quantile(y, RealType(log(y)), xi, Policy());
}

} // namespace detail

BOOST_MATH_EXPORT template <class RealType, class Policy>
inline RealType cdf(const generalized_extreme_value_distribution<RealType, Policy>& dist, const RealType& x)
{
   return detail::gev_cdf_imp(dist, x, false, false, "boost::math::cdf(const generalized_extreme_value_distribution<%1%>&, %1%)");
}

BOOST_MATH_EXPORT template <class RealType, class Policy>
inline RealType cdf(const complemented2_type<generalized_extreme_value_distribution<RealType, Policy>, RealType>& c)
{
   return detail::gev_cdf_imp(c.dist, c.param, true, false, "boost::math::cdf(const generalized_extreme_value_distribution<%1%>&, %1%)");
}

BOOST_MATH_EXPORT template <class RealType, class Policy>
inline RealType logcdf(const generalized_extreme_value_distribution<RealType, Policy>& dist, const RealType& x)
{
   return detail::gev_cdf_imp(dist, x, false, true, "boost::math::logcdf(const generalized_extreme_value_distribution<%1%>&, %1%)");
}

BOOST_MATH_EXPORT template <class RealType, class Policy>
inline RealType logcdf(const complemented2_type<generalized_extreme_value_distribution<RealType, Policy>, RealType>& c)
{
   return detail::gev_cdf_imp(c.dist, c.param, true, true, "boost::math::logcdf(const generalized_extreme_value_distribution<%1%>&, %1%)");
}

BOOST_MATH_EXPORT template <class RealType, class Policy>
inline RealType quantile(const generalized_extreme_value_distribution<RealType, Policy>& dist, const RealType& p)
{
   return detail::gev_quantile_imp(dist, p, false, "boost::math::quantile(const generalized_extreme_value_distribution<%1%>&, %1%)");
}

BOOST_MATH_EXPORT template <class RealType, class Policy>
inline RealType quantile(const complemented2_type<generalized_extreme_value_distribution<RealType, Policy>, RealType>& c)
{
   return detail::gev_quantile_imp(c.dist, c.param, true, "boost::math::quantile(const generalized_extreme_value_distribution<%1%>&, %1%)");
}

BOOST_MATH_EXPORT template <class RealType, class Policy>
inline RealType mean(const generalized_extreme_value_distribution<RealType, Policy>& dist)
{
   constexpr auto function = "boost::math::mean(const generalized_extreme_value_distribution<%1%>&)";
   RealType result = 0;
   if (!detail::check_gev(function, dist, &result))
   {
      return result;
   }
   RealType xi = dist.shape();
   if (xi >= 1)
   {
      // The right tail decays like x^(-1/xi), so the mean is +infinity:
      return detail::gev_infinity<RealType>();
   }
   if (detail::gev_is_gumbel(xi))
   {
      return mean(extreme_value_distribution<RealType, Policy>(dist.location(), dist.scale()));
   }
   // mu + sigma (Gamma(1 - xi) - 1)/xi:
   return dist.location() + dist.scale() * BOOST_MATH_NAMESPACE::tgamma1pm1(-xi, Policy()) / xi;
}

BOOST_MATH_EXPORT template <class RealType, class Policy>
inline RealType variance(const generalized_extreme_value_distribution<RealType, Policy>& dist)
{
   constexpr auto function = "boost::math::variance(const generalized_extreme_value_distribution<%1%>&)";
   RealType result = 0;
   if (!detail::check_gev(function, dist, &result))
   {
      return result;
   }
   RealType xi = dist.shape();
   if (2 * xi >= 1)
   {
      // The second moment diverges, so the variance is +infinity:
      return detail::gev_infinity<RealType>();
   }
   RealType sigma = dist.scale();
   if (detail::gev_is_gumbel(xi))
   {
      return sigma * sigma * constants::pi_sqr<RealType>() / 6;
   }
   // sigma^2 (Gamma(1 - 2 xi) - Gamma(1 - xi)^2)/xi^2 = (sigma Gamma(1 - xi)/xi)^2 expm1(A_2):
   const int a2[3] = { 1, 0, 0 };
   RealType g1 = BOOST_MATH_NAMESPACE::tgamma(1 - xi, Policy());
   RealType s = sigma * g1 / xi;
   return s * s * BOOST_MATH_NAMESPACE::expm1(detail::gev_log_moment_combination(a2, xi, Policy()), Policy());
}

BOOST_MATH_EXPORT template <class RealType, class Policy>
inline RealType standard_deviation(const generalized_extreme_value_distribution<RealType, Policy>& dist)
{
   BOOST_MATH_STD_USING
   return sqrt(variance(dist));
}

BOOST_MATH_EXPORT template <class RealType, class Policy>
inline RealType skewness(const generalized_extreme_value_distribution<RealType, Policy>& dist)
{
   BOOST_MATH_STD_USING
   constexpr auto function = "boost::math::skewness(const generalized_extreme_value_distribution<%1%>&)";
   RealType result = 0;
   if (!detail::check_gev(function, dist, &result))
   {
      return result;
   }
   RealType xi = dist.shape();
   if (3 * xi >= 1)
   {
      return detail::gev_undefined_moment<RealType, Policy>(function, "The skewness is undefined for shape >= 1/3, but got %1%.", xi);
   }
   if (detail::gev_is_gumbel(xi))
   {
      return skewness(extreme_value_distribution<RealType, Policy>());
   }
   // sign(xi) (g3 - 3 g1 g2 + 2 g1^3)/(g2 - g1^2)^(3/2) with g_k = Gamma(1 - k xi); over g1^3 the numerator
   // is expm1(A_3) - 3 expm1(A_2) = (A_3 - 3 A_2) + expm1mx(A_3) - 3 expm1mx(A_2), which is O(xi^3):
   const int a2[3] = { 1, 0, 0 };
   const int a3[3] = { 0, 1, 0 };
   const int b[3] = { -3, 1, 0 };
   RealType A2 = detail::gev_log_moment_combination(a2, xi, Policy());
   RealType A3 = detail::gev_log_moment_combination(a3, xi, Policy());
   RealType numerator = detail::gev_log_moment_combination(b, xi, Policy()) + detail::gev_expm1mx(A3, Policy()) - 3 * detail::gev_expm1mx(A2, Policy());
   RealType v = BOOST_MATH_NAMESPACE::expm1(A2, Policy());
   RealType skew = numerator / (v * sqrt(v));
   return xi > 0 ? skew : RealType(-skew);
}

BOOST_MATH_EXPORT template <class RealType, class Policy>
inline RealType kurtosis_excess(const generalized_extreme_value_distribution<RealType, Policy>& dist)
{
   constexpr auto function = "boost::math::kurtosis_excess(const generalized_extreme_value_distribution<%1%>&)";
   RealType result = 0;
   if (!detail::check_gev(function, dist, &result))
   {
      return result;
   }
   RealType xi = dist.shape();
   if (4 * xi >= 1)
   {
      return detail::gev_undefined_moment<RealType, Policy>(function, "The kurtosis is undefined for shape >= 1/4, but got %1%.", xi);
   }
   if (detail::gev_is_gumbel(xi))
   {
      return kurtosis_excess(extreme_value_distribution<RealType, Policy>());
   }
   // (g4 - 4 g1 g3 + 6 g1^2 g2 - 3 g1^4)/(g2 - g1^2)^2 - 3; over g1^4 the numerator is
   // expm1(A_4) - 4 expm1(A_3) + 6 expm1(A_2), whose part linear in the A_k is O(xi^4):
   const int a2[3] = { 1, 0, 0 };
   const int a3[3] = { 0, 1, 0 };
   const int a4[3] = { 0, 0, 1 };
   const int c[3] = { 6, -4, 1 };
   RealType A2 = detail::gev_log_moment_combination(a2, xi, Policy());
   RealType A3 = detail::gev_log_moment_combination(a3, xi, Policy());
   RealType A4 = detail::gev_log_moment_combination(a4, xi, Policy());
   RealType numerator = detail::gev_log_moment_combination(c, xi, Policy())
      + detail::gev_expm1mx(A4, Policy()) - 4 * detail::gev_expm1mx(A3, Policy()) + 6 * detail::gev_expm1mx(A2, Policy());
   RealType v = BOOST_MATH_NAMESPACE::expm1(A2, Policy());
   return numerator / (v * v) - 3;
}

BOOST_MATH_EXPORT template <class RealType, class Policy>
inline RealType kurtosis(const generalized_extreme_value_distribution<RealType, Policy>& dist)
{
   return kurtosis_excess(dist) + 3;
}

BOOST_MATH_EXPORT template <class RealType, class Policy>
inline RealType mode(const generalized_extreme_value_distribution<RealType, Policy>& dist)
{
   BOOST_MATH_STD_USING
   constexpr auto function = "boost::math::mode(const generalized_extreme_value_distribution<%1%>&)";
   RealType result = 0;
   if (!detail::check_gev(function, dist, &result))
   {
      return result;
   }
   RealType xi = dist.shape();
   if (xi <= -1)
   {
      // The density increases to the upper endpoint:
      return dist.location() - dist.scale() / xi;
   }
   // mu + sigma ((1 + xi)^(-xi) - 1)/xi, where t = 1 + xi:
   return dist.location() + dist.scale() * detail::gev_standard_quantile(RealType(1 + xi), RealType(BOOST_MATH_NAMESPACE::log1p(xi, Policy())), xi, Policy());
}

BOOST_MATH_EXPORT template <class RealType, class Policy>
inline RealType median(const generalized_extreme_value_distribution<RealType, Policy>& dist)
{
   constexpr auto function = "boost::math::median(const generalized_extreme_value_distribution<%1%>&)";
   RealType result = 0;
   if (!detail::check_gev(function, dist, &result))
   {
      return result;
   }
   // t = log 2:
   return dist.location() + dist.scale() * detail::gev_standard_quantile(constants::ln_two<RealType>(), constants::ln_ln_two<RealType>(), dist.shape(), Policy());
}

BOOST_MATH_EXPORT template <class RealType, class Policy>
inline RealType entropy(const generalized_extreme_value_distribution<RealType, Policy>& dist)
{
   BOOST_MATH_STD_USING
   constexpr auto function = "boost::math::entropy(const generalized_extreme_value_distribution<%1%>&)";
   RealType result = 0;
   if (!detail::check_gev(function, dist, &result))
   {
      return result;
   }
   return log(dist.scale()) + constants::euler<RealType>() * (1 + dist.shape()) + 1;
}

BOOST_MATH_NAMESPACE_END

// This include must be at the end, *after* the accessors
// for this distribution have been defined, in order to
// keep compilers that support two-phase lookup happy.
#include <boost/math/distributions/detail/derived_accessors.hpp>

#endif // BOOST_MATH_DISTRIBUTIONS_GENERALIZED_EXTREME_VALUE_HPP
