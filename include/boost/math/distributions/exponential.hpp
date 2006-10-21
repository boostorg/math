//  Copyright John Maddock 2006.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_STATS_EXPONENTIAL_HPP
#define BOOST_STATS_EXPONENTIAL_HPP

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

namespace boost{ namespace math{

namespace detail{
//
// Error check:
//
template <class RealType>
inline bool verify_lambda(const char* function, RealType l, RealType* presult)
{
   if(l <= 0)
   {
      *presult = tools::domain_error<RealType>(
         function,
         "The scale parameter \"lambda\" must be > 0, but was: %1%.", l);
      return false;
   }
   return true;
}

template <class RealType>
inline bool verify_exp_x(const char* function, RealType x, RealType* presult)
{
   if(x < 0)
   {
      *presult = tools::domain_error<RealType>(
         function,
         "The random variable must be >= 0, but was: %1%.", x);
      return false;
   }
   return true;
}

} // namespace detail

template <class RealType = double>
class exponential_distribution
{
public:
   typedef RealType value_type;

   exponential_distribution(RealType lambda = 1)
      : m_lambda(lambda)
   {
      RealType err;
      detail::verify_lambda(BOOST_CURRENT_FUNCTION, lambda, &err);
   } // exponential_distribution

   RealType lambda()const { return m_lambda; }

private:
   RealType m_lambda;
};

typedef exponential_distribution<double> exponential;

template <class RealType>
RealType pdf(const exponential_distribution<RealType>& dist, const RealType& x)
{
   using namespace std; // for ADL of std functions

   RealType lambda = dist.lambda();
   RealType result;
   if(0 == detail::verify_lambda(BOOST_CURRENT_FUNCTION, lambda, &result))
      return result;
   if(0 == detail::verify_exp_x(BOOST_CURRENT_FUNCTION, x, &result))
      return result;
   result = lambda * exp(-lambda * x);
   return result;
} // pdf

template <class RealType>
RealType cdf(const exponential_distribution<RealType>& dist, const RealType& x)
{
   using namespace std; // for ADL of std functions

   RealType result;
   RealType lambda = dist.lambda();
   if(0 == detail::verify_lambda(BOOST_CURRENT_FUNCTION, lambda, &result))
      return result;
   if(0 == detail::verify_exp_x(BOOST_CURRENT_FUNCTION, x, &result))
      return result;
   result = -boost::math::expm1(-x * lambda);

   return result;
} // cdf

template <class RealType>
RealType quantile(const exponential_distribution<RealType>& dist, const RealType& p)
{
   using namespace std; // for ADL of std functions

   RealType result;
   RealType lambda = dist.lambda();
   if(0 == detail::verify_lambda(BOOST_CURRENT_FUNCTION, lambda, &result))
      return result;
   if(0 == detail::check_probability(BOOST_CURRENT_FUNCTION, p, &result))
      return result;

   if(p == 0)
      return 0;
   if(p == 1)
      return tools::overflow_error<RealType>(BOOST_CURRENT_FUNCTION, 0);

   result = -boost::math::log1p(-p) / lambda;
   return result;
} // quantile

template <class RealType>
RealType cdf(const complemented2_type<exponential_distribution<RealType>, RealType>& c)
{
   using namespace std; // for ADL of std functions

   RealType result;
   RealType lambda = c.dist.lambda();
   if(0 == detail::verify_lambda(BOOST_CURRENT_FUNCTION, lambda, &result))
      return result;
   if(0 == detail::verify_exp_x(BOOST_CURRENT_FUNCTION, c.param, &result))
      return result;
   result = exp(-c.param * lambda);

   return result;
}

template <class RealType>
RealType quantile(const complemented2_type<exponential_distribution<RealType>, RealType>& c)
{
   using namespace std; // for ADL of std functions

   RealType result;
   RealType lambda = c.dist.lambda();
   if(0 == detail::verify_lambda(BOOST_CURRENT_FUNCTION, lambda, &result))
      return result;

   RealType q = c.param;
   if(0 == detail::check_probability(BOOST_CURRENT_FUNCTION, q, &result))
      return result;

   if(q == 1)
      return 0;
   if(q == 0)
      return tools::overflow_error<RealType>(BOOST_CURRENT_FUNCTION, 0);

   result = -log(q) / lambda;
   return result;
}

template <class RealType>
inline RealType mean(const exponential_distribution<RealType>& dist)
{
   RealType result;
   RealType lambda = dist.lambda();
   if(0 == detail::verify_lambda(BOOST_CURRENT_FUNCTION, lambda, &result))
      return result;
   return 1 / lambda;
}

template <class RealType>
inline RealType standard_deviation(const exponential_distribution<RealType>& dist)
{
   RealType result;
   RealType lambda = dist.lambda();
   if(0 == detail::verify_lambda(BOOST_CURRENT_FUNCTION, lambda, &result))
      return result;
   return 1 / lambda;
}

template <class RealType>
inline RealType mode(const exponential_distribution<RealType>& /*dist*/)
{
   return 0;
}

template <class RealType>
inline RealType skewness(const exponential_distribution<RealType>& /*dist*/)
{
   return 2;
}

template <class RealType>
inline RealType kurtosis(const exponential_distribution<RealType>& /*dist*/)
{
   return 9;
}

template <class RealType>
inline RealType kurtosis_excess(const exponential_distribution<RealType>& /*dist*/)
{
   return 6;
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

#endif // BOOST_STATS_EXPONENTIAL_HPP
