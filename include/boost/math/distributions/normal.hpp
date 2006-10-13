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

#include <boost/math/special_functions/erf.hpp> // for erf/erfc.
#include <boost/math/distributions/complement.hpp>

namespace boost{ namespace math{

template <class RealType>
class normal_distribution
{
public:
   typedef RealType value_type;

   normal_distribution(RealType mean = 0, RealType sd = 1)
      : m_mean(mean), m_sd(sd) {}

   RealType mean()const
   {
      return m_mean;
   }

   RealType standard_deviation()const
   {
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

template <class RealType>
RealType pdf(const normal_distribution<RealType>& dist, const RealType& x)
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

template <class RealType>
RealType cdf(const normal_distribution<RealType>& dist, const RealType& x)
{
   using namespace std;  // for ADL of std functions

   RealType sd = dist.standard_deviation();
   RealType mean = dist.mean();

   RealType diff = (x - mean) / (sd * constants::root_two<RealType>());
   RealType result;

   result = boost::math::erfc(-diff) / 2;

   return result;
}

template <class RealType>
RealType quantile(const normal_distribution<RealType>& dist, const RealType& p)
{
   using namespace std;  // for ADL of std functions

   RealType sd = dist.standard_deviation();
   RealType mean = dist.mean();

   RealType r;

   r = boost::math::erfc_inv(2 * p);
   r = -r;

   r *= sd * constants::root_two<RealType>();
   r += mean;

   return r;
}

template <class RealType>
RealType cdf(const complemented2_type<normal_distribution<RealType>, RealType>& c)
{
   using namespace std;  // for ADL of std functions

   RealType sd = c.dist.standard_deviation();
   RealType mean = c.dist.mean();
   RealType x = c.param;

   RealType diff = (x - mean) / (sd * constants::root_two<RealType>());
   RealType result;

   result = boost::math::erfc(diff) / 2;

   return result;
}

template <class RealType>
RealType quantile(const complemented2_type<normal_distribution<RealType>, RealType>& c)
{
   using namespace std;  // for ADL of std functions

   RealType sd = c.dist.standard_deviation();
   RealType mean = c.dist.mean();
   RealType q = c.param;

   RealType r;

   r = boost::math::erfc_inv(2 * q);

   r *= sd * constants::root_two<RealType>();
   r += mean;

   return r;
}

template <class RealType>
inline RealType mean(const normal_distribution<RealType>& dist)
{
   return dist.mean();
}

template <class RealType>
inline RealType standard_deviation(const normal_distribution<RealType>& dist)
{
   return dist.standard_deviation();
}

template <class RealType>
inline RealType mode(const normal_distribution<RealType>& dist)
{
   return dist.mean();
}

template <class RealType>
inline RealType skewness(const normal_distribution<RealType>& dist)
{
   return 0;
}

template <class RealType>
inline RealType kurtosis(const normal_distribution<RealType>& dist)
{
   return 3;
}

template <class RealType>
inline RealType kurtosis_excess(const normal_distribution<RealType>& dist)
{
   return 0;
}

} // namespace math
} // namespace boost

// This include must be at the end, *after* the accessors
// for this distribution have been defined, in order to
// keep compilers that support two-phase lookup happy.
#include <boost/math/distributions/detail/derived_accessors.hpp>

#endif // BOOST_STATS_STUDENTS_T_HPP

