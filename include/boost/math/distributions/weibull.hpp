//  Copyright John Maddock 2006.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_STATS_WEIBULL_HPP
#define BOOST_STATS_WEIBULL_HPP

// http://www.itl.nist.gov/div898/handbook/eda/section3/eda3668.htm
// http://mathworld.wolfram.com/WeibullDistribution.html

#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/special_functions/expm1.hpp>
#include <boost/math/distributions/detail/common_error_handling.hpp>
#include <boost/math/distributions/complement.hpp>

#include <utility>

namespace boost{ namespace math
{
namespace detail{

template <class RealType>
inline bool check_weibull_shape(
      const char* function,
      RealType shape,
      RealType* result)
{
   if((shape < 0) || !(boost::math::isfinite)(shape))
   {
      *result = tools::domain_error<RealType>(
         function,
         "Shape parameter is %1%, but must be > 0 !", shape);
      return false;
   }
   return true;
}

template <class RealType>
inline bool check_weibull_x(
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
inline bool check_weibull(
      const char* function,
      RealType scale,
      RealType shape,
      RealType* result)
{
   return check_scale(function, scale, result) && check_weibull_shape(function, shape, result);
}

} // namespace detail

template <class RealType = double>
class weibull_distribution
{
public:
   typedef RealType value_type;

   weibull_distribution(RealType shape, RealType scale = 1)
      : m_shape(shape), m_scale(scale)
   {
      RealType result;
      detail::check_weibull(BOOST_CURRENT_FUNCTION, scale, shape, &result);
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

typedef weibull_distribution<double> weibull;

template <class RealType>
inline const std::pair<RealType, RealType> range(const weibull_distribution<RealType>& /*dist*/)
{ // Range of permissible values for random variable x.
	using boost::math::tools::max_value;
	return std::pair<RealType, RealType>(0, max_value<RealType>());
}

template <class RealType>
inline const std::pair<RealType, RealType> support(const weibull_distribution<RealType>& /*dist*/)
{ // Range of supported values for random variable x.
	// This is range where cdf rises from 0 to 1, and outside it, the pdf is zero.
	using boost::math::tools::max_value;
	return std::pair<RealType, RealType>(0,  max_value<RealType>());
}

template <class RealType>
inline RealType pdf(const weibull_distribution<RealType>& dist, const RealType& x)
{
   using namespace std;  // for ADL of std functions

   RealType shape = dist.shape();
   RealType scale = dist.scale();

   RealType result;
   if(false == detail::check_weibull(BOOST_CURRENT_FUNCTION, scale, shape, &result))
      return result;
   if(false == detail::check_weibull_x(BOOST_CURRENT_FUNCTION, x, &result))
      return result;

   if(x == 0)
      return 0;

   result = exp(-pow(x / scale, shape));
   result *= pow(x / scale, shape) * shape / x;

   return result;
}

template <class RealType>
inline RealType cdf(const weibull_distribution<RealType>& dist, const RealType& x)
{
   using namespace std;  // for ADL of std functions

   RealType shape = dist.shape();
   RealType scale = dist.scale();

   RealType result;
   if(false == detail::check_weibull(BOOST_CURRENT_FUNCTION, scale, shape, &result))
      return result;
   if(false == detail::check_weibull_x(BOOST_CURRENT_FUNCTION, x, &result))
      return result;

   result = -boost::math::expm1(-pow(x / scale, shape));

   return result;
}

template <class RealType>
inline RealType quantile(const weibull_distribution<RealType>& dist, const RealType& p)
{
   using namespace std;  // for ADL of std functions

   RealType shape = dist.shape();
   RealType scale = dist.scale();

   RealType result;
   if(false == detail::check_weibull(BOOST_CURRENT_FUNCTION, scale, shape, &result))
      return result;
   if(false == detail::check_probability(BOOST_CURRENT_FUNCTION, p, &result))
      return result;

   if(p == 1)
      return tools::overflow_error<RealType>(BOOST_CURRENT_FUNCTION);

   result = scale * pow(-boost::math::log1p(-p), 1 / shape);

   return result;
}

template <class RealType>
inline RealType cdf(const complemented2_type<weibull_distribution<RealType>, RealType>& c)
{
   using namespace std;  // for ADL of std functions

   RealType shape = c.dist.shape();
   RealType scale = c.dist.scale();

   RealType result;
   if(false == detail::check_weibull(BOOST_CURRENT_FUNCTION, scale, shape, &result))
      return result;
   if(false == detail::check_weibull_x(BOOST_CURRENT_FUNCTION, c.param, &result))
      return result;

   result = exp(-pow(c.param / scale, shape));

   return result;
}

template <class RealType>
inline RealType quantile(const complemented2_type<weibull_distribution<RealType>, RealType>& c)
{
   using namespace std;  // for ADL of std functions

   RealType shape = c.dist.shape();
   RealType scale = c.dist.scale();
   RealType q = c.param;

   RealType result;
   if(false == detail::check_weibull(BOOST_CURRENT_FUNCTION, scale, shape, &result))
      return result;
   if(false == detail::check_probability(BOOST_CURRENT_FUNCTION, q, &result))
      return result;

   if(q == 0)
      return tools::overflow_error<RealType>(BOOST_CURRENT_FUNCTION);

   result = scale * pow(-log(q), 1 / shape);

   return result;
}

template <class RealType>
inline RealType mean(const weibull_distribution<RealType>& dist)
{
   using namespace std;  // for ADL of std functions

   RealType shape = dist.shape();
   RealType scale = dist.scale();

   RealType result;
   if(false == detail::check_weibull(BOOST_CURRENT_FUNCTION, scale, shape, &result))
      return result;

   result = scale * boost::math::tgamma(1 + 1 / shape);
   return result;
}

template <class RealType>
inline RealType variance(const weibull_distribution<RealType>& dist)
{
   RealType shape = dist.shape();
   RealType scale = dist.scale();

   RealType result;
   if(false == detail::check_weibull(BOOST_CURRENT_FUNCTION, scale, shape, &result))
   {
      return result;
   }
   result = boost::math::tgamma(1 + 1 / shape);
   result *= -result;
   result += boost::math::tgamma(1 + 2 / shape);
   result *= scale * scale;
   return result;
}

template <class RealType>
inline RealType mode(const weibull_distribution<RealType>& dist)
{
   using namespace std;  // for ADL of std function pow.

   RealType shape = dist.shape();
   RealType scale = dist.scale();

   RealType result;
   if(false == detail::check_weibull(BOOST_CURRENT_FUNCTION, scale, shape, &result))
   {
      return result;
   }
   result = scale * pow((shape - 1) / shape, 1 / shape);
   return result;
}

template <class RealType>
inline RealType median(const weibull_distribution<RealType>& dist)
{
   using namespace std;  // for ADL of std function pow.

   RealType shape = dist.shape(); // Wikipedia k
   RealType scale = dist.scale(); // Wikipedia lambda

   RealType result;
   if(false == detail::check_weibull(BOOST_CURRENT_FUNCTION, scale, shape, &result))
   {
      return result;
   }
   using boost::math::constants::ln_two;
   result = scale * pow(ln_two<RealType>(), 1 / shape);
   return result;
}

template <class RealType>
inline RealType skewness(const weibull_distribution<RealType>& dist)
{
   using namespace std;  // for ADL of std functions

   RealType shape = dist.shape();
   RealType scale = dist.scale();

   RealType result;
   if(false == detail::check_weibull(BOOST_CURRENT_FUNCTION, scale, shape, &result))
   {
      return result;
   }
   RealType g1, g2, g3, d;

   g1 = boost::math::tgamma(1 + 1 / shape);
   g2 = boost::math::tgamma(1 + 2 / shape);
   g3 = boost::math::tgamma(1 + 3 / shape);
   d = pow(g2 - g1 * g1, RealType(1.5));

   result = (2 * g1 * g1 * g1 - 3 * g1 * g2 + g3) / d;
   return result;
}

template <class RealType>
inline RealType kurtosis_excess(const weibull_distribution<RealType>& dist)
{
   using namespace std;  // for ADL of std functions

   RealType shape = dist.shape();
   RealType scale = dist.scale();

   RealType result;
   if(false == detail::check_weibull(BOOST_CURRENT_FUNCTION, scale, shape, &result))
      return result;

   RealType g1, g2, g3, g4, d, g1_2, g1_4;

   g1 = boost::math::tgamma(1 + 1 / shape);
   g2 = boost::math::tgamma(1 + 2 / shape);
   g3 = boost::math::tgamma(1 + 3 / shape);
   g4 = boost::math::tgamma(1 + 4 / shape);
   g1_2 = g1 * g1;
   g1_4 = g1_2 * g1_2;
   d = g2 - g1_2;
   d *= d;

   result = -6 * g1_4 + 12 * g1_2 * g2 - 3 * g2 * g2 - 4 * g1 * g3 + g4;
   result /= d;
   return result;
}

template <class RealType>
inline RealType kurtosis(const weibull_distribution<RealType>& dist)
{
   return kurtosis_excess(dist) + 3;
}

} // namespace math
} // namespace boost

// This include must be at the end, *after* the accessors
// for this distribution have been defined, in order to
// keep compilers that support two-phase lookup happy.
#include <boost/math/distributions/detail/derived_accessors.hpp>

#endif // BOOST_STATS_WEIBULL_HPP

