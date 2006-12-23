//  Copyright John Maddock 2006.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_STATS_GAMMA_HPP
#define BOOST_STATS_GAMMA_HPP

// http://www.itl.nist.gov/div898/handbook/eda/section3/eda366b.htm
// http://mathworld.wolfram.com/GammaDistribution.html
// http://en.wikipedia.org/wiki/Gamma_distribution

#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/distributions/detail/common_error_handling.hpp>
#include <boost/math/distributions/complement.hpp>

#include <utility>

namespace boost{ namespace math
{
namespace detail
{

template <class RealType>
bool check_gamma_shape(
      const char* function,
      RealType shape,
      RealType* result)
{
   if((shape <= 0) || !(boost::math::isfinite)(shape))
   {
      *result = tools::domain_error<RealType>(
         function,
         "Shape parameter is %1%, but must be > 0 !", shape);
      return false;
   }
   return true;
}

template <class RealType>
bool check_gamma_x(
      const char* function,
      RealType const& x,
      RealType* result)
{
   if((x < 0) || !(boost::math::isfinite)(x))
   {
      *result = tools::domain_error<RealType>(
         function,
         "Random variate is %1% but must be >= 0 !", x);
      return false;
   }
   return true;
}

template <class RealType>
inline bool check_gamma(
      const char* function,
      RealType scale,
      RealType shape,
      RealType* result)
{
   return check_scale(function, scale, result) && check_gamma_shape(function, shape, result);
}

} // namespace detail

template <class RealType = double>
class gamma_distribution
{
public:
   typedef RealType value_type;

   gamma_distribution(RealType shape, RealType scale = 1)
      : m_shape(shape), m_scale(scale)
   {
      RealType result;
      detail::check_gamma(BOOST_CURRENT_FUNCTION, scale, shape, &result);
   }

   RealType shape()const
   {
      return m_shape;
   }

   RealType scale()const
   {
      return m_scale;
   }
private:
   //
   // Data members:
   //
   RealType m_shape;     // distribution shape
   RealType m_scale;     // distribution scale
};

// NO typedef because of clash with name of gamma function.

template <class RealType>
const std::pair<RealType, RealType> range(const gamma_distribution<RealType>& dist)
{ // Range of permissible values for random variable x.
	using boost::math::tools::max_value;
	return std::pair<RealType, RealType>(0, max_value<RealType>());
}

template <class RealType>
const std::pair<RealType, RealType> support(const gamma_distribution<RealType>& dist)
{ // Range of supported values for random variable x.
	// This is range where cdf rises from 0 to 1, and outside it, the pdf is zero.
	using boost::math::tools::max_value;
	return std::pair<RealType, RealType>(0,  max_value<RealType>());
}

template <class RealType>
RealType pdf(const gamma_distribution<RealType>& dist, const RealType& x)
{
   using namespace std;  // for ADL of std functions

   RealType shape = dist.shape();
   RealType scale = dist.scale();

   RealType result;
   if(false == detail::check_gamma(BOOST_CURRENT_FUNCTION, scale, shape, &result))
      return result;
   if(false == detail::check_gamma_x(BOOST_CURRENT_FUNCTION, x, &result))
      return result;

   if(x == 0)
   {
      return 0;
   }
   result = gamma_p_derivative(shape, x / scale) / scale;
   return result;
} // pdf

template <class RealType>
inline RealType cdf(const gamma_distribution<RealType>& dist, const RealType& x)
{
   using namespace std;  // for ADL of std functions

   RealType shape = dist.shape();
   RealType scale = dist.scale();

   RealType result;
   if(false == detail::check_gamma(BOOST_CURRENT_FUNCTION, scale, shape, &result))
      return result;
   if(false == detail::check_gamma_x(BOOST_CURRENT_FUNCTION, x, &result))
      return result;

   result = boost::math::gamma_p(shape, x / scale);
   return result;
} // cdf

template <class RealType>
RealType quantile(const gamma_distribution<RealType>& dist, const RealType& p)
{
   using namespace std;  // for ADL of std functions

   RealType shape = dist.shape();
   RealType scale = dist.scale();

   RealType result;
   if(false == detail::check_gamma(BOOST_CURRENT_FUNCTION, scale, shape, &result))
      return result;
   if(false == detail::check_probability(BOOST_CURRENT_FUNCTION, p, &result))
      return result;

   if(p == 1)
      return tools::overflow_error<RealType>(BOOST_CURRENT_FUNCTION, 0);

   result = gamma_p_inv(shape, p) * scale;

   return result;
}

template <class RealType>
RealType cdf(const complemented2_type<gamma_distribution<RealType>, RealType>& c)
{
   using namespace std;  // for ADL of std functions

   RealType shape = c.dist.shape();
   RealType scale = c.dist.scale();

   RealType result;
   if(false == detail::check_gamma(BOOST_CURRENT_FUNCTION, scale, shape, &result))
      return result;
   if(false == detail::check_gamma_x(BOOST_CURRENT_FUNCTION, c.param, &result))
      return result;

   result = gamma_q(shape, c.param / scale);

   return result;
}

template <class RealType>
RealType quantile(const complemented2_type<gamma_distribution<RealType>, RealType>& c)
{
   using namespace std;  // for ADL of std functions

   RealType shape = c.dist.shape();
   RealType scale = c.dist.scale();
   RealType q = c.param;

   RealType result;
   if(false == detail::check_gamma(BOOST_CURRENT_FUNCTION, scale, shape, &result))
      return result;
   if(false == detail::check_probability(BOOST_CURRENT_FUNCTION, q, &result))
      return result;

   if(q == 0)
      return tools::overflow_error<RealType>(BOOST_CURRENT_FUNCTION, 0);

   result = gamma_q_inv(shape, q) * scale;

   return result;
}

template <class RealType>
inline RealType mean(const gamma_distribution<RealType>& dist)
{
   using namespace std;  // for ADL of std functions

   RealType shape = dist.shape();
   RealType scale = dist.scale();

   RealType result;
   if(false == detail::check_gamma(BOOST_CURRENT_FUNCTION, scale, shape, &result))
      return result;

   result = shape * scale;
   return result;
}

template <class RealType>
RealType variance(const gamma_distribution<RealType>& dist)
{
   using namespace std;  // for ADL of std functions

   RealType shape = dist.shape();
   RealType scale = dist.scale();

   RealType result;
   if(false == detail::check_gamma(BOOST_CURRENT_FUNCTION, scale, shape, &result))
      return result;

   result = shape * scale * scale;
   return result;
}

template <class RealType>
inline RealType mode(const gamma_distribution<RealType>& dist)
{
   using namespace std;  // for ADL of std functions

   RealType shape = dist.shape();
   RealType scale = dist.scale();

   RealType result;
   if(false == detail::check_gamma(BOOST_CURRENT_FUNCTION, scale, shape, &result))
      return result;

   if(shape < 1)
      return tools::domain_error<RealType>(
         BOOST_CURRENT_FUNCTION,
         "The mode of the gamma distribution is only defined for values of the shape parameter >= 1, but got %1%.",
         shape);

   result = (shape - 1) * scale;
   return result;
}

//template <class RealType>
//inline RealType median(const gamma_distribution<RealType>& dist)
//{  // Rely on default definition in derived accessors.
//}

template <class RealType>
inline RealType skewness(const gamma_distribution<RealType>& dist)
{
   using namespace std;  // for ADL of std functions

   RealType shape = dist.shape();
   RealType scale = dist.scale();

   RealType result;
   if(false == detail::check_gamma(BOOST_CURRENT_FUNCTION, scale, shape, &result))
      return result;

   result = 2 / sqrt(shape);
   return result;
}

template <class RealType>
inline RealType kurtosis_excess(const gamma_distribution<RealType>& dist)
{
   using namespace std;  // for ADL of std functions

   RealType shape = dist.shape();
   RealType scale = dist.scale();

   RealType result;
   if(false == detail::check_gamma(BOOST_CURRENT_FUNCTION, scale, shape, &result))
      return result;

   result = 6 / shape;
   return result;
}

template <class RealType>
inline RealType kurtosis(const gamma_distribution<RealType>& dist)
{
   return kurtosis_excess(dist) + 3;
}

} // namespace math
} // namespace boost

// This include must be at the end, *after* the accessors
// for this distribution have been defined, in order to
// keep compilers that support two-phase lookup happy.
#include <boost/math/distributions/detail/derived_accessors.hpp>

#endif // BOOST_STATS_GAMMA_HPP

