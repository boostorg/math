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
class cauchy_distribution
{
public:
   typedef RealType value_type;

   cauchy_distribution()
   {
   } // cauchy_distribution

private:
};

typedef cauchy_distribution<double> cauchy;

template <class RealType>
RealType pdf(const cauchy_distribution<RealType>& dist, const RealType& x)
{
   RealType result = 1 / (constants::pi<RealType>() * (1 + x * x));
   return result;
} // pdf

template <class RealType>
RealType cdf(const cauchy_distribution<RealType>& dist, const RealType& x)
{
   using namespace std; // for ADL of std functions

   RealType mx = -fabs(x);

   // special case first:
   if(mx > -tools::epsilon<RealType>() / 8)
      return 0.5;

   RealType result = -atan(1 / mx) / constants::pi<RealType>();

   return (x > 0 ? 1	- result : result);
} // cdf

template <class RealType>
RealType quantile(const cauchy_distribution<RealType>& dist, const RealType& p)
{
   using namespace std; // for ADL of std functions
   // Special cases:
   if(p == 1)
      return tools::overflow_error<RealType>(BOOST_CURRENT_FUNCTION, 0);
   if(p == 0)
      return -tools::overflow_error<RealType>(BOOST_CURRENT_FUNCTION, 0);

   RealType result;
   if(0 == detail::check_probability(BOOST_CURRENT_FUNCTION, p, &result))
      return result;

   // argument reduction of p:
   RealType P = p - floor(p);
   if(P > 0.5)
      P = P - 1;
   // special case:
   if(P == 0.5)
      return 0;
   result = -1 / tan(constants::pi<RealType>() * P);
   return result;
} // quantile

template <class RealType>
RealType cdf(const complemented2_type<cauchy_distribution<RealType>, RealType>& c)
{
   return cdf(c.dist, -c.param);
}

template <class RealType>
RealType quantile(const complemented2_type<cauchy_distribution<RealType>, RealType>& c)
{
   return -quantile(c.dist, c.param);
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
   return 0;
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
