// Copyright John Maddock 2006.

// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_DISTRIBUTIONS_CHI_SQUARED_HPP
#define BOOST_MATH_DISTRIBUTIONS_CHI_SQUARED_HPP

#include <boost/math/special_functions/gamma.hpp> // for incomplete beta.
#include <boost/math/distributions/complement.hpp> // complements
#include <boost/math/distributions/detail/common_error_handling.hpp> // error checks
#include <boost/math/special_functions/fpclassify.hpp>

namespace boost{ namespace math{ 
   
namespace detail{


} // namespace detail

template <class RealType>
class chi_squared_distribution
{
public:
   typedef RealType value_type;

   chi_squared_distribution(RealType i) : m_df(i)
   {
      RealType result;
      detail::check_df(
         BOOST_CURRENT_FUNCTION, m_df, &result);
   } // chi_squared_distribution

   RealType degrees_of_freedom()const
   {
      return m_df;
   }

   // Parameter estimation:
   static RealType estimate_degrees_of_freedom(
      RealType difference_from_mean,
      RealType alpha,
      RealType beta,
      RealType sd,
      RealType hint = 100);

private:
   //
   // Data members:
   //
   RealType m_df;  // degrees of freedom are a real number.
};

typedef chi_squared_distribution<double> chi_squared;

template <class RealType>
RealType pdf(const chi_squared_distribution<RealType>& dist, const RealType& chi_square)
{
   using namespace std;  // for ADL of std functions
   RealType degrees_of_freedom = dist.degrees_of_freedom();
   // Error check:
   RealType error_result;
   if(false == detail::check_df(
         BOOST_CURRENT_FUNCTION, degrees_of_freedom, &error_result))
      return error_result;

   if((chi_square <= 0) || !(boost::math::isfinite)(chi_square))
   {
      return tools::domain_error<RealType>(
         BOOST_CURRENT_FUNCTION, "Chi Square parameter was %1%, but must be > 0 !", chi_square);
   }

   return gamma_P_derivative(degrees_of_freedom / 2, chi_square / 2) / 2;
} // pdf

template <class RealType>
RealType cdf(const chi_squared_distribution<RealType>& dist, const RealType& chi_square)
{
   RealType degrees_of_freedom = dist.degrees_of_freedom();
   // Error check:
   RealType error_result;
   if(false == detail::check_df(
         BOOST_CURRENT_FUNCTION, degrees_of_freedom, &error_result))
      return error_result;

   if((chi_square <= 0) || !(boost::math::isfinite)(chi_square))
   {
      return tools::domain_error<RealType>(
         BOOST_CURRENT_FUNCTION, "Chi Square parameter was %1%, but must be > 0 !", chi_square);
   }

   return boost::math::gamma_P(degrees_of_freedom / 2, chi_square / 2);
} // cdf

template <class RealType>
RealType quantile(const chi_squared_distribution<RealType>& dist, const RealType& p)
{
   RealType degrees_of_freedom = dist.degrees_of_freedom();
   // Error check:
   RealType error_result;
   if(false == detail::check_df(
         BOOST_CURRENT_FUNCTION, degrees_of_freedom, &error_result)
         && detail::check_probability(
            BOOST_CURRENT_FUNCTION, p, &error_result))
      return error_result;

   return 2 * boost::math::gamma_P_inv(degrees_of_freedom / 2, p);
} // quantile

template <class RealType>
RealType cdf(const complemented2_type<chi_squared_distribution<RealType>, RealType>& c)
{
   RealType const& degrees_of_freedom = c.dist.degrees_of_freedom();
   RealType const& chi_square = c.param;
   // Error check:
   RealType error_result;
   if(false == detail::check_df(
         BOOST_CURRENT_FUNCTION, degrees_of_freedom, &error_result))
      return error_result;

   if((chi_square <= 0) || !(boost::math::isfinite)(chi_square))
   {
      return tools::domain_error<RealType>(
         BOOST_CURRENT_FUNCTION, "Chi Square parameter was %1%, but must be > 0 !", chi_square);
   }

   return boost::math::gamma_Q(degrees_of_freedom / 2, chi_square / 2);
}

template <class RealType>
RealType quantile(const complemented2_type<chi_squared_distribution<RealType>, RealType>& c)
{
   RealType const& degrees_of_freedom = c.dist.degrees_of_freedom();
   RealType const& q = c.param;
   // Error check:
   RealType error_result;
   if(false == detail::check_df(
         BOOST_CURRENT_FUNCTION, degrees_of_freedom, &error_result)
         && detail::check_probability(
            BOOST_CURRENT_FUNCTION, q, &error_result))
      return error_result;

   return 2 * boost::math::gamma_Q_inv(degrees_of_freedom / 2, q);
}

template <class RealType>
inline RealType mean(const chi_squared_distribution<RealType>& dist)
{ // Mean of Chi-Squared distribution = v.
  return dist.degrees_of_freedom();
} // mean

template <class RealType>
inline RealType variance(const chi_squared_distribution<RealType>& dist)
{ // Variance of Chi-Squared distribution = 2v.
  return 2 * dist.degrees_of_freedom();
} // variance


} // namespace math
} // namespace boost

// This include must be at the end, *after* the accessors
// for this distribution have been defined, in order to
// keep compilers that support two-phase lookup happy.
#include <boost/math/distributions/detail/derived_accessors.hpp>

#endif // BOOST_MATH_DISTRIBUTIONS_CHI_SQUARED_HPP