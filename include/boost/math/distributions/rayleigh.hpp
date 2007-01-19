//  Copyright Paul A. Bristow 2007.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_STATS_rayleigh_HPP
#define BOOST_STATS_rayleigh_HPP

#include <boost/math/constants/constants.hpp>
#include <boost/math/special_functions/log1p.hpp>
#include <boost/math/special_functions/expm1.hpp>
#include <boost/math/distributions/complement.hpp>
#include <boost/math/distributions/detail/common_error_handling.hpp>
#include <cmath>

#ifdef BOOST_MSVC
# pragma warning(push)
# pragma warning(disable: 4702) // unreachable code (return after domain_error throw).
#endif

#include <utility>

namespace boost{ namespace math{

namespace detail
{ // Error checks:
  template <class RealType>
  inline bool verify_sigma(const char* function, RealType sigma, RealType* presult)
  {
     if(sigma <= 0)
     {
        *presult = tools::domain_error<RealType>(
           function,
           "The scale parameter \"sigma\" must be > 0, but was: %1%.", sigma);
        return false;
     }
     return true;
  } // bool verify_sigma

  template <class RealType>
  inline bool verify_rayleigh_x(const char* function, RealType x, RealType* presult)
  {
     if(x < 0)
     {
        *presult = tools::domain_error<RealType>(
           function,
           "The random variable must be >= 0, but was: %1%.", x);
        return false;
     }
     return true;
  } // bool verify_rayleigh_x
} // namespace detail

template <class RealType = double>
class rayleigh_distribution
{
public:
   typedef RealType value_type;

   rayleigh_distribution(RealType sigma = 1)
      : m_sigma(sigma)
   {
      RealType err;
      detail::verify_sigma(BOOST_CURRENT_FUNCTION, sigma, &err);
   } // rayleigh_distribution

   RealType sigma()const
   { // Accessor.
     return m_sigma;
   }

private:
   RealType m_sigma;
}; // class rayleigh_distribution

typedef rayleigh_distribution<double> rayleigh;

template <class RealType>
const std::pair<RealType, RealType> range(const rayleigh_distribution<RealType>& /*dist*/)
{ // Range of permissible values for random variable x.
	using boost::math::tools::max_value;
	return std::pair<RealType, RealType>(static_cast<RealType>(1), max_value<RealType>());
}

template <class RealType>
const std::pair<RealType, RealType> support(const rayleigh_distribution<RealType>& /*dist*/)
{ // Range of supported values for random variable x.
	// This is range where cdf rises from 0 to 1, and outside it, the pdf is zero.
	using boost::math::tools::max_value;
	return std::pair<RealType, RealType>((1),  max_value<RealType>());
}

template <class RealType>
RealType pdf(const rayleigh_distribution<RealType>& dist, const RealType& x)
{
   using namespace std; // for ADL of std function exp.

   RealType sigma = dist.sigma();
   RealType result;
   if(false == detail::verify_sigma(BOOST_CURRENT_FUNCTION, sigma, &result))
   {
      return result;
   }
   if(false == detail::verify_rayleigh_x(BOOST_CURRENT_FUNCTION, x, &result))
   {
      return result;
   }
   RealType sigmasqr = sigma * sigma;
   result = x * (exp(-(x * x) / ( 2 * sigmasqr))) / sigmasqr; 
   return result;
} // pdf

template <class RealType>
RealType cdf(const rayleigh_distribution<RealType>& dist, const RealType& x)
{
   using namespace std; // for ADL of std functions

   RealType result;
   RealType sigma = dist.sigma();
   if(false == detail::verify_sigma(BOOST_CURRENT_FUNCTION, sigma, &result))
   {
      return result;
   }
   if(false == detail::verify_rayleigh_x(BOOST_CURRENT_FUNCTION, x, &result))
   {
      return result;
   }
   result = -boost::math::expm1(-x * x / ( 2 * sigma * sigma));
   return result;
} // cdf

template <class RealType>
RealType quantile(const rayleigh_distribution<RealType>& dist, const RealType& p)
{
   using namespace std; // for ADL of std functions

   RealType result;
   RealType sigma = dist.sigma();
   if(false == detail::verify_sigma(BOOST_CURRENT_FUNCTION, sigma, &result))
      return result;
   if(false == detail::check_probability(BOOST_CURRENT_FUNCTION, p, &result))
      return result;

   if(p == 0)
   {
      return 0;
   }
   if(p == 1)
   {
     return tools::overflow_error<RealType>(BOOST_CURRENT_FUNCTION);
   }
   result = sqrt(-2 * sigma * sigma * boost::math::log1p(-p));
   return result;
} // quantile

template <class RealType>
RealType cdf(const complemented2_type<rayleigh_distribution<RealType>, RealType>& c)
{
   using namespace std; // for ADL of std functions

   RealType result;
   RealType sigma = c.dist.sigma();
   if(false == detail::verify_sigma(BOOST_CURRENT_FUNCTION, sigma, &result))
   {
      return result;
   }
   RealType x = c.param;
   if(false == detail::verify_rayleigh_x(BOOST_CURRENT_FUNCTION, x, &result))
   {
      return result;
   }
   result =  exp(-x * x / ( 2 * sigma * sigma));
   return result;
} // cdf complement

template <class RealType>
RealType quantile(const complemented2_type<rayleigh_distribution<RealType>, RealType>& c)
{
   using namespace std; // for ADL of std functions, log & sqrt.

   RealType result;
   RealType sigma = c.dist.sigma();
   if(false == detail::verify_sigma(BOOST_CURRENT_FUNCTION, sigma, &result))
   {
      return result;
   }
   RealType q = c.param;
   if(false == detail::check_probability(BOOST_CURRENT_FUNCTION, q, &result))
   {
      return result;
   }
   if(q == 1)
   {
      return 0;
   }
   if(q == 0)
   {
     return tools::overflow_error<RealType>(BOOST_CURRENT_FUNCTION);
   }
   result = sqrt(-2 * sigma * sigma * log(q));
   return result;
} // quantile complement

template <class RealType>
inline RealType mean(const rayleigh_distribution<RealType>& dist)
{
   RealType result;
   RealType sigma = dist.sigma();
   if(false == detail::verify_sigma(BOOST_CURRENT_FUNCTION, sigma, &result))
   {
      return result;
   }
   using boost::math::constants::root_half_pi;
   return sigma * root_half_pi<RealType>();
} // mean

template <class RealType>
inline RealType variance(const rayleigh_distribution<RealType>& dist)
{
   RealType result;
   RealType sigma = dist.sigma();
   if(false == detail::verify_sigma(BOOST_CURRENT_FUNCTION, sigma, &result))
   {
      return result;
   }
   using boost::math::constants::four_minus_pi;
   return four_minus_pi<RealType>() * sigma * sigma / 2;
} // variance

template <class RealType>
inline RealType mode(const rayleigh_distribution<RealType>& dist)
{
   return dist.sigma();
}

template <class RealType>
inline RealType median(const rayleigh_distribution<RealType>& dist)
{
   using boost::math::constants::root_ln_four;
   return root_ln_four<RealType>() * dist.sigma();
}

template <class RealType>
inline RealType skewness(const rayleigh_distribution<RealType>& /*dist*/)
{
  // using namespace boost::math::constants;
  return static_cast<RealType>(0.63111065781893713819189935154422777984404221106391L);
  // Computed using NTL at 150 bit, about 50 decimal digits.
  // return 2 * root_pi<RealType>() * pi_minus_three<RealType>() / pow23_four_minus_pi<RealType>();
}

template <class RealType>
inline RealType kurtosis(const rayleigh_distribution<RealType>& /*dist*/)
{
  // using namespace boost::math::constants;
  return static_cast<RealType>(3.2450893006876380628486604106197544154170667057995L);
  // Computed using NTL at 150 bit, about 50 decimal digits.
  // return 3 - (6 * pi<RealType>() * pi<RealType>() - 24 * pi<RealType>() + 16) /
  // (four_minus_pi<RealType>() * four_minus_pi<RealType>());
}

template <class RealType>
inline RealType kurtosis_excess(const rayleigh_distribution<RealType>& /*dist*/)
{
  //using namespace boost::math::constants;
  // Computed using NTL at 150 bit, about 50 decimal digits.
  return static_cast<RealType>(0.2450893006876380628486604106197544154170667057995L);
  // return -(6 * pi<RealType>() * pi<RealType>() - 24 * pi<RealType>() + 16) /
  //   (four_minus_pi<RealType>() * four_minus_pi<RealType>());
} // kurtosis

} // namespace math
} // namespace boost

#ifdef BOOST_MSVC
# pragma warning(pop)
#endif

// This include must be at the end, *after* the accessors
// for this distribution have been defined, in order to
// keep compilers that support two-phase lookup happy.
#include <boost/math/distributions/detail/derived_accessors.hpp>

#endif // BOOST_STATS_rayleigh_HPP
