//  Copyright John Maddock 2006.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_STATS_NORMAL_HPP
#define BOOST_STATS_NORMAL_HPP

// http://en.wikipedia.org/wiki/Normal_distribution
// http://www.itl.nist.gov/div898/handbook/eda/section3/eda3661.htm
// Also:
// Weisstein, Eric W. "Normal Distribution."
// From MathWorld--A Wolfram Web Resource.
// http://mathworld.wolfram.com/NormalDistribution.html

#include <boost/math/distributions/fwd.hpp>
#include <boost/math/special_functions/erf.hpp> // for erf/erfc.
#include <boost/math/distributions/complement.hpp>

#include <utility>

namespace boost{ namespace math{

template <class RealType = double, class Policy = policies::policy<> >
class normal_distribution
{
public:
   typedef RealType value_type;
   typedef Policy policy_type;

   normal_distribution(RealType mean = 0, RealType sd = 1)
      : m_mean(mean), m_sd(sd) {}

   RealType mean()const
   { // location
      return m_mean;
   }

   RealType standard_deviation()const
   { // scale
      return m_sd;
   }
private:
   //
   // Data members:
   //
   RealType m_mean;  // distribution mean
   RealType m_sd;    // distribution standard deviation
};

typedef normal_distribution<double> normal;

template <class RealType, class Policy>
inline const std::pair<RealType, RealType> range(const normal_distribution<RealType, Policy>& /*dist*/)
{ // Range of permissible values for random variable x.
	using boost::math::tools::max_value;
	return std::pair<RealType, RealType>(-max_value<RealType>(), max_value<RealType>()); // - to + infinity.
}

template <class RealType, class Policy>
inline const std::pair<RealType, RealType> support(const normal_distribution<RealType, Policy>& /*dist*/)
{ // Range of supported values for random variable x.
	// This is range where cdf rises from 0 to 1, and outside it, the pdf is zero.
	using boost::math::tools::max_value;
	return std::pair<RealType, RealType>(-max_value<RealType>(),  max_value<RealType>()); // - to + infinity.
}

template <class RealType, class Policy>
inline RealType pdf(const normal_distribution<RealType, Policy>& dist, const RealType& x)
{
   using namespace std;  // for ADL of std functions

   RealType sd = dist.standard_deviation();
   RealType mean = dist.mean();

   RealType exponent = x - mean;
   exponent *= -exponent;
   exponent /= 2 * sd * sd;

   RealType result = exp(exponent);
   result /= sd * sqrt(2 * constants::pi<RealType>());

   return result;
}

template <class RealType, class Policy>
inline RealType cdf(const normal_distribution<RealType, Policy>& dist, const RealType& x)
{
   using namespace std;  // for ADL of std functions

   RealType sd = dist.standard_deviation();
   RealType mean = dist.mean();

   RealType diff = (x - mean) / (sd * constants::root_two<RealType>());
   RealType result;

   result = boost::math::erfc(-diff, Policy()) / 2;
   return result;
}

template <class RealType, class Policy>
inline RealType quantile(const normal_distribution<RealType, Policy>& dist, const RealType& p)
{
   using namespace std;  // for ADL of std functions

   RealType sd = dist.standard_deviation();
   RealType mean = dist.mean();

   RealType r;

   r = boost::math::erfc_inv(2 * p, Policy());
   r = -r;
   r *= sd * constants::root_two<RealType>();
   r += mean;

   return r;
}

template <class RealType, class Policy>
inline RealType cdf(const complemented2_type<normal_distribution<RealType, Policy>, RealType>& c)
{
   using namespace std;  // for ADL of std functions

   RealType sd = c.dist.standard_deviation();
   RealType mean = c.dist.mean();
   RealType x = c.param;

   RealType diff = (x - mean) / (sd * constants::root_two<RealType>());
   RealType result;

   result = boost::math::erfc(diff, Policy()) / 2;

   return result;
}

template <class RealType, class Policy>
inline RealType quantile(const complemented2_type<normal_distribution<RealType, Policy>, RealType>& c)
{
   using namespace std;  // for ADL of std functions

   RealType sd = c.dist.standard_deviation();
   RealType mean = c.dist.mean();
   RealType q = c.param;

   RealType r;
   r = boost::math::erfc_inv(2 * q, Policy());
   r *= sd * constants::root_two<RealType>();
   r += mean;
   return r;
}

template <class RealType, class Policy>
inline RealType mean(const normal_distribution<RealType, Policy>& dist)
{
   return dist.mean();
}

template <class RealType, class Policy>
inline RealType standard_deviation(const normal_distribution<RealType, Policy>& dist)
{
   return dist.standard_deviation();
}

template <class RealType, class Policy>
inline RealType mode(const normal_distribution<RealType, Policy>& dist)
{
   return dist.mean();
}

template <class RealType, class Policy>
inline RealType median(const normal_distribution<RealType, Policy>& dist)
{
   return dist.mean();
}

template <class RealType, class Policy>
inline RealType skewness(const normal_distribution<RealType, Policy>& /*dist*/)
{
   return 0;
}

template <class RealType, class Policy>
inline RealType kurtosis(const normal_distribution<RealType, Policy>& /*dist*/)
{
   return 3;
}

template <class RealType, class Policy>
inline RealType kurtosis_excess(const normal_distribution<RealType, Policy>& /*dist*/)
{
   return 0;
}

} // namespace math
} // namespace boost

// This include must be at the end, *after* the accessors
// for this distribution have been defined, in order to
// keep compilers that support two-phase lookup happy.
#include <boost/math/distributions/detail/derived_accessors.hpp>

#endif // BOOST_STATS_NORMAL_HPP


