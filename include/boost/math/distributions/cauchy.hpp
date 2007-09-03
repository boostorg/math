//  Copyright John Maddock 2006.
// Copyright Paul A. Bristow 2007.

//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_STATS_CAUCHY_HPP
#define BOOST_STATS_CAUCHY_HPP

#include <boost/math/distributions/fwd.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/math/distributions/complement.hpp>
#include <boost/math/distributions/detail/common_error_handling.hpp>
#include <cmath>

#ifdef BOOST_MSVC
# pragma warning(push)
# pragma warning(disable: 4702) // unreachable code (return after raise_domain_error throw).
#endif

#include <utility>

namespace boost{ namespace math
{

template <class RealType, class Policy>
class cauchy_distribution;

namespace detail
{

template <class RealType, class Policy>
inline bool check_cauchy_scale(const char* func, RealType scale, RealType* result, const Policy& pol)
{
   if(scale <= 0)
   {
      *result = policies::raise_domain_error<RealType>(
         func,
         "The scale parameter for the Cauchy distribution must be > 0 but got %1%.",
         scale, pol);
      return false;
   }
   return true;
}

template <class RealType, class Policy>
RealType cdf_imp(const cauchy_distribution<RealType, Policy>& dist, const RealType& x, bool complement)
{
   //
   // This calculates the cdf of the Cauchy distribution and/or its complement.
   //
   // The usual formula for the Cauchy cdf is:
   //
   // cdf = 0.5 + atan(x)/pi
   //
   // But that suffers from cancellation error as x -> -INF.
   //
   // Recall that for x < 0:
   //
   // atan(x) = -pi/2 - atan(1/x)
   //
   // Substituting into the above we get:
   //
   // CDF = -atan(1/x)  ; x < 0
   //
   // So the proceedure is to calculate the cdf for -fabs(x)
   // using the above formula, and then subtract from 1 when required
   // to get the result.
   //
   BOOST_MATH_STD_USING // for ADL of std functions

   RealType result;
   RealType loc = dist.location();
   RealType scale = dist.scale();
   if(0 == detail::check_cauchy_scale("boost::math::cdf(cauchy<%1%>&, %1%)", scale, &result, Policy()))
      return result;
   RealType mx = -fabs((x - loc) / scale);

   // special case first:
   if(mx > -tools::epsilon<RealType>() / 8)
      return 0.5;

   result = -atan(1 / mx) / constants::pi<RealType>();

   return (((x > loc) != complement) ? 1	- result : result);
} // cdf

template <class RealType, class Policy>
RealType quantile_imp(
      const cauchy_distribution<RealType, Policy>& dist,
      const RealType& p,
      bool complement)
{
   //
   // This routine implements the quantile for the Cauchy distibution,
   // the value p may be the probability or it's complement if complement=true.
   //
   // The proceedure first performs argument reduction on p to avoid error
   // when calculating the tangent, then calulates the distance from the
   // mid-point of the distribution.  This is either added or subtracted
   // from the location parameter depending on whether `complement` is true.
   //
   static const char* function = "boost::math::quantile(cauchy<%1%>&, %1%)";
   BOOST_MATH_STD_USING // for ADL of std functions
   // Special cases:
   if(p == 1)
      return (complement ? -1 : 1) * policies::raise_overflow_error<RealType>(function, 0, Policy());
   if(p == 0)
      return (complement ? 1 : -1) * policies::raise_overflow_error<RealType>(function, 0, Policy());

   RealType result;
   RealType loc = dist.location();
   RealType scale = dist.scale();
   if(0 == detail::check_cauchy_scale(function, scale, &result, Policy()))
      return result;
   if(0 == detail::check_probability(function, p, &result, Policy()))
      return result;

   // argument reduction of p:
   RealType P = p - floor(p);
   if(P > 0.5)
      P = P - 1;
   // special case:
   if(P == 0.5)
      return loc;
   result = -scale / tan(constants::pi<RealType>() * P);
   return complement ? loc - result : loc + result;
} // quantile

} // namespace detail

template <class RealType = double, class Policy = policies::policy<> >
class cauchy_distribution
{
public:
   typedef RealType value_type;
   typedef Policy policy_type;

   cauchy_distribution(RealType a = 0, RealType b = 1)
      : m_a(a), m_hg(b)
   {
      RealType r;
      detail::check_cauchy_scale("boost::math::cauchy<%1%>::cauchy", b, &r, Policy());
   } // cauchy_distribution

   RealType location()const
   {
      return m_a;
   }
   RealType scale()const
   {
      return m_hg;
   }

private:
   RealType m_a;    // The location, this is the median of the distribution
   RealType m_hg;   // The scale, this is the half width at half height.
};

typedef cauchy_distribution<double> cauchy;

template <class RealType, class Policy>
inline const std::pair<RealType, RealType> range(const cauchy_distribution<RealType, Policy>&)
{ // Range of permissible values for random variable x.
	using boost::math::tools::max_value;
	return std::pair<RealType, RealType>(-max_value<RealType>(), max_value<RealType>()); // - to + infinity.
}

template <class RealType, class Policy>
inline const std::pair<RealType, RealType> support(const cauchy_distribution<RealType, Policy>& )
{ // Range of supported values for random variable x.
	// This is range where cdf rises from 0 to 1, and outside it, the pdf is zero.
   return std::pair<RealType, RealType>(-tools::max_value<RealType>(), tools::max_value<RealType>()); // - to + infinity.
}

template <class RealType, class Policy>
inline RealType pdf(const cauchy_distribution<RealType, Policy>& dist, const RealType& x)
{
   RealType result;
   RealType loc = dist.location();
   RealType scale = dist.scale();
   if(0 == detail::check_cauchy_scale("boost::math::pdf(cauchy<%1%>&, %1%)", scale, &result, Policy()))
      return result;

   RealType xs = (x - loc) / scale;

   result = 1 / (constants::pi<RealType>() * scale * (1 + xs * xs));
   return result;
} // pdf

template <class RealType, class Policy>
inline RealType cdf(const cauchy_distribution<RealType, Policy>& dist, const RealType& x)
{
   return detail::cdf_imp(dist, x, false);
} // cdf

template <class RealType, class Policy>
inline RealType quantile(const cauchy_distribution<RealType, Policy>& dist, const RealType& p)
{
   return detail::quantile_imp(dist, p, false);
} // quantile

template <class RealType, class Policy>
inline RealType cdf(const complemented2_type<cauchy_distribution<RealType, Policy>, RealType>& c)
{
   return detail::cdf_imp(c.dist, c.param, true);
}

template <class RealType, class Policy>
inline RealType quantile(const complemented2_type<cauchy_distribution<RealType, Policy>, RealType>& c)
{
   return detail::quantile_imp(c.dist, c.param, true);
}

template <class RealType, class Policy>
inline RealType mean(const cauchy_distribution<RealType, Policy>&)
{  // There is no mean:
   typedef typename Policy::assert_undefined_type assert_type;
   BOOST_STATIC_ASSERT(assert_type::value == 0);

   return policies::raise_domain_error<RealType>(
      "boost::math::mean(cauchy<%1%>&)",
      "The Cauchy distribution does not have a mean: "
      "the only possible return value is %1%.",
      std::numeric_limits<RealType>::quiet_NaN(), Policy());
}

template <class RealType, class Policy>
inline RealType variance(const cauchy_distribution<RealType, Policy>& /*dist*/)
{
   // There is no variance:
   typedef typename Policy::assert_undefined_type assert_type;
   BOOST_STATIC_ASSERT(assert_type::value == 0);

   return policies::raise_domain_error<RealType>(
      "boost::math::variance(cauchy<%1%>&)",
      "The Cauchy distribution does not have a variance: "
      "the only possible return value is %1%.",
      std::numeric_limits<RealType>::quiet_NaN(), Policy());
}

template <class RealType, class Policy>
inline RealType mode(const cauchy_distribution<RealType, Policy>& dist)
{
   return dist.location();
}

template <class RealType, class Policy>
inline RealType median(const cauchy_distribution<RealType, Policy>& dist)
{
   return dist.location();
}
template <class RealType, class Policy>
inline RealType skewness(const cauchy_distribution<RealType, Policy>& /*dist*/)
{
   // There is no skewness:
   typedef typename Policy::assert_undefined_type assert_type;
   BOOST_STATIC_ASSERT(assert_type::value == 0);

   return policies::raise_domain_error<RealType>(
      "boost::math::skewness(cauchy<%1%>&)",
      "The Cauchy distribution does not have a skewness: "
      "the only possible return value is %1%.",
      std::numeric_limits<RealType>::quiet_NaN(), Policy()); // infinity?
}

template <class RealType, class Policy>
inline RealType kurtosis(const cauchy_distribution<RealType, Policy>& /*dist*/)
{
   // There is no kurtosis:
   typedef typename Policy::assert_undefined_type assert_type;
   BOOST_STATIC_ASSERT(assert_type::value == 0);

   return policies::raise_domain_error<RealType>(
      "boost::math::kurtosis(cauchy<%1%>&)",
      "The Cauchy distribution does not have a kurtosis: "
      "the only possible return value is %1%.",
      std::numeric_limits<RealType>::quiet_NaN(), Policy());
}

template <class RealType, class Policy>
inline RealType kurtosis_excess(const cauchy_distribution<RealType, Policy>& /*dist*/)
{
   // There is no kurtosis excess:
   typedef typename Policy::assert_undefined_type assert_type;
   BOOST_STATIC_ASSERT(assert_type::value == 0);

   return policies::raise_domain_error<RealType>(
      "boost::math::kurtosis_excess(cauchy<%1%>&)",
      "The Cauchy distribution does not have a kurtosis: "
      "the only possible return value is %1%.",
      std::numeric_limits<RealType>::quiet_NaN(), Policy());
}

} // namespace math
} // namespace boost

#ifdef BOOST_MSVC
# pragma warning(pop)
#endif

// This include must be at the end, *after* the accessors
// for this distribution have been defined, in order to
// keep compilers that support two-phase lookup happy.
#include <boost/math/distributions/detail/derived_accessors.hpp>

#endif // BOOST_STATS_CAUCHY_HPP
