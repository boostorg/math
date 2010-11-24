//  Copyright John Maddock 2010.
//  Copyright Paul A. Bristow 2010.

//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_STATS_INVERSE_NORMAL_HPP
#define BOOST_STATS_INVERSE_NORMAL_HPP

// http://en.wikipedia.org/wiki/Normal-inverse_Gaussian_distribution
// http://mathworld.wolfram.com/InverseGaussianDistribution.html

// The normal-inverse Gaussian distribution (also called the Wald distribution when mean = 1)
// is the continuous probability distribution
// that is defined as the normal variance-mean mixture where the mixing density is the 
// inverse Gaussian distribution. The tails of the distribution decrease more slowly
// than the normal distribution. It is therefore suitable to model phenomena
// where numerically large values are more probable than is the case for the normal distribution.

// Examples are returns from financial assets and turbulent wind speeds. 
// The normal-inverse Gaussian distributions form
// a subclass of the generalised hyperbolic distributions.

// See also

// http://en.wikipedia.org/wiki/Normal_distribution
// http://www.itl.nist.gov/div898/handbook/eda/section3/eda3661.htm
// Also:
// Weisstein, Eric W. "Normal Distribution."
// From MathWorld--A Wolfram Web Resource.
// http://mathworld.wolfram.com/NormalDistribution.html

// http://www.jstatsoft.org/v26/i04/paper General class of inverse Gaussian distributions.
// ig package - withdrawn but at http://cran.r-project.org/src/contrib/Archive/ig/

// http://www.stat.ucl.ac.be/ISdidactique/Rhelp/library/SuppDists/html/invGauss.html
// R package for dinvGauss, ...

#include <boost/math/distributions/fwd.hpp>
#include <boost/math/special_functions/erf.hpp> // for erf/erfc.
#include <boost/math/distributions/complement.hpp>
#include <boost/math/distributions/detail/common_error_handling.hpp>
#include <boost/math/distributions/normal.hpp>

#include <utility>

namespace boost{ namespace math{

template <class RealType = double, class Policy = policies::policy<> >
class inverse_normal_distribution
{
public:
   typedef RealType value_type;
   typedef Policy policy_type;

   inverse_normal_distribution(RealType mean = 1, RealType sd = 1)
      : m_mean(mean), m_sd(sd)
   { // Default is a 1,1 inverse_normal distribution.
     static const char* function = "boost::math::inverse_normal_distribution<%1%>::inverse_normal_distribution";

     RealType result;
     detail::check_scale(function, sd, &result, Policy());
     detail::check_location(function, mean, &result, Policy());
   }

   RealType mean()const
   { // alias for location.
      return m_mean; // aka mu
   }

   RealType standard_deviation()const
   { // alias for scale.
      return m_sd; // aka lambda.
   }

   // Synonyms, provided to allow generic use of find_location and find_scale.
   RealType location()const
   { // location, aka mu.
      return m_mean;
   }
   RealType scale()const
   { // scale, aka lambda.
      return m_sd;
   }

private:
   //
   // Data members:
   //
   RealType m_mean;  // distribution mean or location, aka mu.
   RealType m_sd;    // distribution standard deviation or scale, aka lambda.
}; // class normal_distribution

typedef inverse_normal_distribution<double> inverse_normal;

template <class RealType, class Policy>
inline const std::pair<RealType, RealType> range(const inverse_normal_distribution<RealType, Policy>& /*dist*/)
{ // Range of permissible values for random variable x, zero to max.
   using boost::math::tools::max_value;
   return std::pair<RealType, RealType>(static_cast<RealType>(0), max_value<RealType>()); // - to + max value.
}

template <class RealType, class Policy>
inline const std::pair<RealType, RealType> support(const inverse_normal_distribution<RealType, Policy>& /*dist*/)
{ // Range of supported values for random variable x, zero to max.
   // This is range where cdf rises from 0 to 1, and outside it, the pdf is zero.

   using boost::math::tools::max_value;
   return std::pair<RealType, RealType>(static_cast<RealType>(0),  max_value<RealType>()); // - to + max value.
}

template <class RealType, class Policy>
inline RealType pdf(const inverse_normal_distribution<RealType, Policy>& dist, const RealType& x)
{ // Probability Density Function
   BOOST_MATH_STD_USING  // for ADL of std functions

   RealType scale = dist.scale();
   RealType mean = dist.mean();
   RealType result;
   static const char* function = "boost::math::pdf(const inverse_normal_distribution<%1%>&, %1%)";
   if(false == detail::check_scale(function, scale, &result, Policy()))
   {
      return result;
   }
   if(false == detail::check_location(function, mean, &result, Policy()))
   {
      return result;
   }
   if(false == detail::check_x_gt0(function, x, &result, Policy()))
   {
      return numeric_limits<RealType>::quiet_NaN();
   }

   //result =
   //  sqrt(scale / (2 * constants::pi<RealType>() * x * x * x))
   // * exp(-scale * (x - mean) * (x - mean) / (2 * x * mean * mean));

   result =
     sqrt(scale / (constants::two_pi<RealType>() * x * x * x))
    * exp(-scale * (x - mean) * (x - mean) / (2 * x * mean * mean));
   return result;
} // pdf

template <class RealType, class Policy>
inline RealType cdf(const inverse_normal_distribution<RealType, Policy>& dist, const RealType& x)
{ // Cumulative Density Function.
   BOOST_MATH_STD_USING  // for ADL of std functions

   RealType scale = dist.scale();
   RealType mean = dist.mean();
   static const char* function = "boost::math::cdf(const inverse_normal_distribution<%1%>&, %1%)";
   RealType result;
   if(false == detail::check_scale(function, scale, &result, Policy()))
   {
      return result;
   }
   if(false == detail::check_location(function, mean, &result, Policy()))
   {
      return result;
   }
   if(false == detail::check_x_gt0(function, x, &result, Policy()))
   {
     return result;
   }

   //result = 0.5 * (erf(sqrt(scale / x) * (x / mean -1) / sqrt(2.L), Policy()) + 1)
   //  + exp(2 * scale / mean) / 2 
   //  * (1 - erf(sqrt(scale / x) * (x / mean +1) / sqrt(2.L), Policy()));

   result = 0.5 * (erf(sqrt(scale / x) * (x / mean - 1) / constants::root_two<RealType>(), Policy()) + 1)
     + exp(2 * scale / mean) / 2 
     * (1 - erf(sqrt(scale / x) * (x / mean + 1) / constants::root_two<RealType>(), Policy()));

   return result;
} // cdf

template <class RealType, class Policy>
inline RealType quantile(const inverse_normal_distribution<RealType, Policy>& dist, const RealType& x)
{
   BOOST_MATH_STD_USING  // for ADL of std functions

   RealType mean = dist.mean();
   RealType scale = dist.scale();
   static const char* function = "boost::math::quantile(const inverse_normal_distribution<%1%>&, %1%)";

   RealType result;
   if(false == detail::check_scale(function, scale, &result, Policy()))
      return result;
   if(false == detail::check_location(function, mean, &result, Policy()))
      return result;
   if(false == detail::check_probability(function, x, &result, Policy()))
      return result;

   cout << "x " << x << endl;
   RealType a = sqrt(scale / x); // a scale = lambda/x
   RealType b = x / mean; // b = x/mu

   // pnorm q, mean, sd, lower.tail = true;

  //double q=1.0-pnorm(+a*(b-1.0), 0, 1, true, false);
  //double p=    pnorm(-a*(b+1.0), 0, 1, true, false);

   //boost::math::normal_distribution<RealType> norm01;
   using boost::math::normal;
   normal norm01;

   double qx = a * (b - 1.0);
   RealType q = 1 - ((qx <= 0) ? 0 : cdf(norm01, qx));
   cout << "a = " << a << ", b = " <<  b << ", qx = " << qx << ", pnorm= " << pnorm01(qx) << ", cdf= "  << cdf(norm01, qx) << " q = " << q << endl;

   //cout << "1 - pnorm01(qx) " << 1.0 - pnorm01(qx) << endl;



   RealType px =  -a * (b + 1.0);
   RealType p =  pnorm01(px);
   RealType cdfpx = (px <= 0) ? 0 : cdf(norm01, px);
   cout << "-a*(b+1.0) == px = " << px <<", pnorm01(p) = " << p << ", cdfpx = " << cdfpx << endl;

   if (p == 0)
   {
     result = q;
   }
   else
   {
     RealType r2 = 2 * scale / mean;
     if (r2 >= numeric_limits<RealType>::max() )
     {
       result = numeric_limits<RealType>::quiet_NaN();
     }
     else
     {
       result = q - exp(r2) * p;
     }
   }
   return result;
} // quantile

template <class RealType, class Policy>
inline RealType cdf(const complemented2_type<inverse_normal_distribution<RealType, Policy>, RealType>& c)
{
   BOOST_MATH_STD_USING  // for ADL of std functions

   RealType sd = c.dist.standard_deviation();
   RealType mean = c.dist.mean();
   RealType x = c.param;
   static const char* function = "boost::math::cdf(const complement(inverse_normal_distribution<%1%>&), %1%)";

   if((boost::math::isinf)(x))
   {
     if(x < 0) return 1; // cdf complement -infinity is unity.
     return 0; // cdf complement +infinity is zero
   }
   // These produce MSVC 4127 warnings, so the above used instead.
   //if(std::numeric_limits<RealType>::has_infinity && x == std::numeric_limits<RealType>::infinity())
   //{ // cdf complement +infinity is zero.
   //  return 0;
   //}
   //if(std::numeric_limits<RealType>::has_infinity && x == -std::numeric_limits<RealType>::infinity())
   //{ // cdf complement -infinity is unity.
   //  return 1;
   //}
   RealType result;
   if(false == detail::check_scale(function, sd, &result, Policy()))
      return result;
   if(false == detail::check_location(function, mean, &result, Policy()))
      return result;
   if(false == detail::check_x(function, x, &result, Policy()))
      return result;

   RealType diff = (x - mean) / (sd * constants::root_two<RealType>());
   result = boost::math::erfc(diff, Policy()) / 2;
   return result;
} // cdf complement

template <class RealType, class Policy>
inline RealType quantile(const complemented2_type<inverse_normal_distribution<RealType, Policy>, RealType>& c)
{
   BOOST_MATH_STD_USING  // for ADL of std functions

   RealType sd = c.dist.standard_deviation();
   RealType mean = c.dist.mean();
   static const char* function = "boost::math::quantile(const complement(inverse_normal_distribution<%1%>&), %1%)";
   RealType result;
   if(false == detail::check_scale(function, sd, &result, Policy()))
      return result;
   if(false == detail::check_location(function, mean, &result, Policy()))
      return result;
   RealType q = c.param;
   if(false == detail::check_probability(function, q, &result, Policy()))
      return result;
   result = boost::math::erfc_inv(2 * q, Policy());
   result *= sd * constants::root_two<RealType>();
   result += mean;
   return result;
} // quantile

template <class RealType, class Policy>
inline RealType mean(const inverse_normal_distribution<RealType, Policy>& dist)
{
   return dist.mean();
}

template <class RealType, class Policy>
inline RealType standard_deviation(const inverse_normal_distribution<RealType, Policy>& dist)
{
  RealType scale = dist.scale();
  RealType mean = dist.mean();
  RealType result = sqrt(mean * mean * mean / scale)
  return result;
}

template <class RealType, class Policy>
inline RealType mode(const inverse_normal_distribution<RealType, Policy>& dist)
{
  RealType scale = dist.scale();
  RealType  mean = dist.mean();
  RealType result = mean * (sqrt(1 + (9 * mean * mean)/(4 * scale * scale)) 
      - 3 * mean / (2 * scale));
  return result;
}

template <class RealType, class Policy>
inline RealType skewness(const inverse_normal_distribution<RealType, Policy>& /*dist*/)
{
   RealType scale = dist.scale();
  RealType  mean = dist.mean();
  RealType result = 3 * sqrt(mean/scale);
  return result;
}

template <class RealType, class Policy>
inline RealType kurtosis(const inverse_normal_distribution<RealType, Policy>& /*dist*/)
{
  RealType scale = dist.scale();
  RealType  mean = dist.mean();
  RealType result = 12 * mean / scale ;
  return result;
}

template <class RealType, class Policy>
inline RealType kurtosis_excess(const inverse_normal_distribution<RealType, Policy>& /*dist*/)
{
  RealType scale = dist.scale();
  RealType  mean = dist.mean();
  RealType result = 15 * mean / scale;
  return result;
}

} // namespace math
} // namespace boost

// This include must be at the end, *after* the accessors
// for this distribution have been defined, in order to
// keep compilers that support two-phase lookup happy.
#include <boost/math/distributions/detail/derived_accessors.hpp>

#endif // BOOST_STATS_INVERSE_NORMAL_HPP


