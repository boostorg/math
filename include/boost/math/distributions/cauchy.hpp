//  Copyright John Maddock 2006.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_STATS_CAUCHY_HPP
#define BOOST_STATS_CAUCHY_HPP

#include <boost/math/constants/constants.hpp>
#include <boost/math/distributions/complement.hpp>
#include <boost/math/distributions/detail/common_error_handling.hpp>
#include <cmath>

#ifdef BOOST_MSVC
# pragma warning(push)
# pragma warning(disable: 4702) // unreachable code (return after domain_error throw).
#endif

namespace boost{ namespace math{

template <class RealType>
class cauchy_distribution;

namespace detail{

template <class RealType>
bool check_cauchy_scale(const char* func, RealType scale, RealType* result)
{
   if(scale <= 0)
   {
      *result = tools::domain_error<RealType>(
         func,
         "The scale parameter for the Cauchy distribution must be > 0 but got %1%.",
         scale);
      return false;
   }
   return true;
}

template <class RealType>
RealType cdf_imp(const cauchy_distribution<RealType>& dist, const RealType& x, bool complement)
{
   using namespace std; // for ADL of std functions

   RealType result;
   RealType loc = dist.location();
   RealType scale = dist.scale();
   if(0 == detail::check_cauchy_scale(BOOST_CURRENT_FUNCTION, scale, &result))
      return result;
   RealType mx = -fabs((x - loc) / scale);

   // special case first:
   if(mx > -tools::epsilon<RealType>() / 8)
      return 0.5;

   result = -atan(1 / mx) / constants::pi<RealType>();

   return (((x > loc) != complement) ? 1	- result : result);
} // cdf

template <class RealType>
RealType quantile_imp(
      const cauchy_distribution<RealType>& dist, 
      const RealType& p, 
      bool complement)
{
   using namespace std; // for ADL of std functions
   // Special cases:
   if(p == 1)
      return (complement ? -1 : 1) * tools::overflow_error<RealType>(BOOST_CURRENT_FUNCTION, 0);
   if(p == 0)
      return (complement ? 1 : -1) * tools::overflow_error<RealType>(BOOST_CURRENT_FUNCTION, 0);

   RealType result;
   RealType loc = dist.location();
   RealType scale = dist.scale();
   if(0 == detail::check_cauchy_scale(BOOST_CURRENT_FUNCTION, scale, &result))
      return result;
   if(0 == detail::check_probability(BOOST_CURRENT_FUNCTION, p, &result))
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

}

template <class RealType>
class cauchy_distribution
{
public:
   typedef RealType value_type;

   cauchy_distribution(RealType a = 0, RealType b = 1)
      : m_a(a), m_hg(b)
   {
      RealType r;
      detail::check_cauchy_scale(BOOST_CURRENT_FUNCTION, b, &r);
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

template <class RealType>
RealType pdf(const cauchy_distribution<RealType>& dist, const RealType& x)
{
   RealType result;
   RealType loc = dist.location();
   RealType scale = dist.scale();
   if(0 == detail::check_cauchy_scale(BOOST_CURRENT_FUNCTION, scale, &result))
      return result;

   RealType xs = (x - loc) / scale;

   result = 1 / (constants::pi<RealType>() * scale * (1 + xs * xs));
   return result;
} // pdf

template <class RealType>
inline RealType cdf(const cauchy_distribution<RealType>& dist, const RealType& x)
{
   return detail::cdf_imp(dist, x, false);
} // cdf

template <class RealType>
inline RealType quantile(const cauchy_distribution<RealType>& dist, const RealType& p)
{
   return detail::quantile_imp(dist, p, false);
} // quantile

template <class RealType>
inline RealType cdf(const complemented2_type<cauchy_distribution<RealType>, RealType>& c)
{
   return detail::cdf_imp(c.dist, c.param, true);
}

template <class RealType>
inline RealType quantile(const complemented2_type<cauchy_distribution<RealType>, RealType>& c)
{
   return detail::quantile_imp(c.dist, c.param, true);
}

template <class RealType>
inline RealType mean(const cauchy_distribution<RealType>& )
{
   return tools::domain_error<RealType>(
      BOOST_CURRENT_FUNCTION,
      "The Cauchy distribution does not have a mean: "
      "the only possible return value is %1%.",
      std::numeric_limits<RealType>::quiet_NaN());
}

template <class RealType>
inline RealType variance(const cauchy_distribution<RealType>& dist)
{
   return tools::domain_error<RealType>(
      BOOST_CURRENT_FUNCTION,
      "The Cauchy distribution does not have a variance: "
      "the only possible return value is %1%.",
      std::numeric_limits<RealType>::quiet_NaN());
}

template <class RealType>
inline RealType mode(const cauchy_distribution<RealType>& dist)
{
   return dist.location();
}

template <class RealType>
inline RealType skewness(const cauchy_distribution<RealType>& dist)
{
   return tools::domain_error<RealType>(
      BOOST_CURRENT_FUNCTION,
      "The Cauchy distribution does not have a skewness: "
      "the only possible return value is %1%.",
      std::numeric_limits<RealType>::quiet_NaN());
}

template <class RealType>
inline RealType kurtosis(const cauchy_distribution<RealType>& dist)
{
   return tools::domain_error<RealType>(
      BOOST_CURRENT_FUNCTION,
      "The Cauchy distribution does not have a kurtosis: "
      "the only possible return value is %1%.",
      std::numeric_limits<RealType>::quiet_NaN());
}

template <class RealType>
inline RealType kurtosis_excess(const cauchy_distribution<RealType>& dist)
{
   return tools::domain_error<RealType>(
      BOOST_CURRENT_FUNCTION,
      "The Cauchy distribution does not have a kurtosis: "
      "the only possible return value is %1%.",
      std::numeric_limits<RealType>::quiet_NaN());
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
