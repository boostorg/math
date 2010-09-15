// Copyright John Maddock 2010.
// Copyright Paul A. Bristow 2010.

// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_DISTRIBUTIONS_INVERSE_CHI_SQUARED_HPP
#define BOOST_MATH_DISTRIBUTIONS_INVERSE_CHI_SQUARED_HPP

#include <boost/math/distributions/fwd.hpp>
#include <boost/math/special_functions/gamma.hpp> // for incomplete beta.
#include <boost/math/distributions/complement.hpp> // complements
#include <boost/math/distributions/detail/common_error_handling.hpp> // error checks
#include <boost/math/special_functions/fpclassify.hpp>

// See http://en.wikipedia.org/wiki/Scaled-inverse-chi-square_distribution
// for definitions of this scaled version.

// http://reference.wolfram.com/mathematica/ref/InverseChiSquareDistribution.html
// Weisstein, Eric W. "Inverse Chi-Squared Distribution." From MathWorld--A Wolfram Web Resource.
// http://mathworld.wolfram.com/InverseChi-SquaredDistribution.html

#include <utility>

namespace boost{ namespace math{

namespace detail
{
  template <class RealType, class Policy>
  inline bool check_inverse_chi_squared( // Check distribution parameters.
        const char* function,
        RealType degrees_of_freedom, // degrees_of_freedom (aka nu).
        RealType scale,  // scale (aka sigma^2)
        RealType* result,
        const Policy& pol)
  {
     return check_scale(function, scale, result, pol)
       && check_df(function, degrees_of_freedom,
       result, pol);
  } // bool check_inverse_chi_squared
} // namespace detail

template <class RealType = double, class Policy = policies::policy<> >
class inverse_chi_squared_distribution
{
public:
   typedef RealType value_type;
   typedef Policy policy_type;

   inverse_chi_squared_distribution(RealType df, RealType scale) : m_df(df), m_scale (scale)
   {
      RealType result;
      detail::check_df(
         "boost::math::inverse_chi_squared_distribution<%1%>::inverse_chi_squared_distribution", m_df, &result, Policy());
   } // inverse_chi_squared_distribution constructor 

   inverse_chi_squared_distribution(RealType df = 1) : m_df(df)
   {
      RealType result;
      m_scale = 1 / m_df ; // Default scale = 1 / degrees of freedom (Wikipedia definition 1).
      detail::check_df(
         "boost::math::inverse_chi_squared_distribution<%1%>::inverse_chi_squared_distribution", m_df, &result, Policy());
   } // inverse_chi_squared_distribution

   RealType degrees_of_freedom()const
   {
      return m_df;
   }
   RealType scale()const
   {
      return m_scale; // aka nu
   }

   // Parameter estimation:
   static RealType find_degrees_of_freedom(
      RealType difference_from_variance,
      RealType alpha,
      RealType beta,
      RealType variance,
      RealType hint = 100);

private:
   // Data members:
   RealType m_df;  // degrees of freedom are treated as a real number.
   RealType m_scale;     // distribution scale

}; // class chi_squared_distribution

typedef inverse_chi_squared_distribution<double> inverse_chi_squared;

template <class RealType, class Policy>
inline const std::pair<RealType, RealType> range(const inverse_chi_squared_distribution<RealType, Policy>& /*dist*/)
{ // Range of permissible values for random variable x.
   using boost::math::tools::max_value;
   return std::pair<RealType, RealType>(static_cast<RealType>(0), max_value<RealType>()); // 0 to + infinity.
}

template <class RealType, class Policy>
inline const std::pair<RealType, RealType> support(const inverse_chi_squared_distribution<RealType, Policy>& /*dist*/)
{ // Range of supported values for random variable x.
   // This is range where cdf rises from 0 to 1, and outside it, the pdf is zero.
   return std::pair<RealType, RealType>(static_cast<RealType>(0), tools::max_value<RealType>()); // 0 to + infinity.
}

template <class RealType, class Policy>
RealType pdf(const inverse_chi_squared_distribution<RealType, Policy>& dist, const RealType& x)
{
   BOOST_MATH_STD_USING  // for ADL of std functions.
   RealType df = dist.degrees_of_freedom();
   RealType scale = dist.scale();
   RealType error_result;

   static const char* function = "boost::math::pdf(const inverse_chi_squared_distribution<%1%>&, %1%)";

   if(false == detail::check_inverse_chi_squared
     (function, df, scale, &error_result, Policy())
     )
   { // Bad distribution.
      return error_result;
   }

   if((x < 0) || !(boost::math::isfinite)(x))
   { // Bad x.
      return policies::raise_domain_error<RealType>(
         function, "inverse Chi Square parameter was %1%, but must be >= 0 !", x, Policy());
   }

   if(x == 0)
   { // Treat as special case.
     return 0;
   }
   // Wikipedia scaled inverse chi sq (df, scale) related to inv gamma (df/2, df * scale /2) 
   // so use inverse gamma pdf with shape = df/2, scale df * scale /2 
   // RealType shape = df /2; // inv_gamma shape
   // RealType scale = df * scale/2; // inv_gamma scale
   //RealType result = gamma_p_derivative(shape, scale / x, Policy()) * scale / (x * x);
   RealType result = gamma_p_derivative(df/2, df * scale/2 / x, Policy()) * df * scale/2 / (x * x);
   return result;
} // pdf

template <class RealType, class Policy>
inline RealType cdf(const inverse_chi_squared_distribution<RealType, Policy>& dist, const RealType& x)
{
   static const char* function = "boost::math::cdf(const inverse_chi_squared_distribution<%1%>&, %1%)";
   RealType df = dist.degrees_of_freedom();
   RealType scale = dist.scale();
   RealType error_result;

   if(false ==
       detail::check_inverse_chi_squared(function, df, scale, &error_result, Policy())
     )
   { // Bad distribution.
      return error_result;
   }
   if((x < 0) || !(boost::math::isfinite)(x))
   { // Bad x.
      return policies::raise_domain_error<RealType>(
         function, "inverse Chi Square parameter was %1%, but must be >= 0 !", x, Policy());
   }
   if (x == 0)
   { // Treat zero as a special case.
     return 0;
   }
   // RealType shape = df /2; // inv_gamma shape
   // RealType scale = df * scale/2; // inv_gamma scale
   // result = boost::math::gamma_q(shape, scale / x, Policy()); // inverse_gamma code
   return boost::math::gamma_q(df / 2, (df * (scale / 2)) / x, Policy());
} // cdf

template <class RealType, class Policy>
inline RealType quantile(const inverse_chi_squared_distribution<RealType, Policy>& dist, const RealType& p)
{
   using boost::math::gamma_q_inv;
   RealType df = dist.degrees_of_freedom();
   RealType scale = dist.scale();

   static const char* function = "boost::math::quantile(const inverse_chi_squared_distribution<%1%>&, %1%)";
   // Error check:
   RealType error_result;
   if(false == detail::check_df(
         function, df, &error_result, Policy())
         && detail::check_probability(
            function, p, &error_result, Policy()))
      return error_result;
   // RealType shape = df /2; // inv_gamma shape
   // RealType scale = df * scale/2; // inv_gamma scale
   // result = scale / gamma_q_inv(shape, p, Policy());
      RealType result = df * scale/2 / gamma_q_inv(df /2, p, Policy());
      return result;
} // quantile

template <class RealType, class Policy>
inline RealType cdf(const complemented2_type<inverse_chi_squared_distribution<RealType, Policy>, RealType>& c)
{
   using boost::math::gamma_q_inv;
   RealType const& df = c.dist.degrees_of_freedom();
   RealType const& scale = c.dist.scale();
   RealType const& x = c.param;
   static const char* function = "boost::math::cdf(const inverse_chi_squared_distribution<%1%>&, %1%)";
   // Error check:
   RealType error_result;
   if(false == detail::check_df(
         function, degrees_of_freedom, &error_result, Policy()))
      return error_result;

   if((x < 0) || !(boost::math::isfinite)(x))
   {
      return policies::raise_domain_error<RealType>(
         function, "inverse Chi Square parameter was %1%, but must be > 0 !", x, Policy());
   }
   // RealType shape = df /2; // inv_gamma shape
   // RealType scale = df * scale/2; // inv_gamma scale
   //   result = gamma_p(shape, scale/c.param, Policy()); inv_gamma

   return gamma_q(df / 2, (df * scale/2)/ x, Policy());
}

template <class RealType, class Policy>
inline RealType quantile(const complemented2_type<inverse_chi_squared_distribution<RealType, Policy>, RealType>& c)
{
   using boost::math::gamma_q_inv;

   RealType const& df = c.dist.degrees_of_freedom();
   RealType const& scale = c.dist.scale();
   RealType const& q = c.param;
   static const char* function = "boost::math::quantile(const inverse_chi_squared_distribution<%1%>&, %1%)";
   // Error check:
   RealType error_result;
   if(false == detail::check_df(
         function, df, &error_result, Policy())
         && detail::check_probability(
            function, q, &error_result, Policy()))
      return error_result;
   // RealType shape = df /2; // inv_gamma shape
   // RealType scale = df * scale/2; // inv_gamma scale
   // result = scale / gamma_p_inv(shape, q, Policy());  // inv_gamma
   return df * scale/2 * gamma_p_inv(df /2, q, Policy());
}
/**/

template <class RealType, class Policy>
inline RealType mean(const inverse_chi_squared_distribution<RealType, Policy>& dist)
{ // Mean of inverse Chi-Squared distribution.
   RealType df = dist.degrees_of_freedom();
   RealType scale = dist.scale();

   static const char* function = "boost::math::mode(const inverse_chi_squared_distribution<%1%>&)";
   if(df <= 2)
      return policies::raise_domain_error<RealType>(
         function,
         "inverse Chi-Squared distribution only has a mode for degrees of freedom > 2, but got degrees of freedom = %1%.",
         df, Policy());
  return (df * scale) / (df - 2);
} // mean

template <class RealType, class Policy>
inline RealType variance(const inverse_chi_squared_distribution<RealType, Policy>& dist)
{ // Variance of inverse Chi-Squared distribution..
   RealType df = dist.degrees_of_freedom();
   RealType scale = dist.scale();
   static const char* function = "boost::math::variance(const inverse_chi_squared_distribution<%1%>&)";
   if(df <= 4)
   {
      return policies::raise_domain_error<RealType>(
         function,
         "inverse Chi-Squared distribution only has a variance for degrees of freedom > 4, but got degrees of freedom = %1%.",
         df, Policy());  return 2 * dist.degrees_of_freedom();
   }
   return 2 * df * df * scale * scale / ((df - 2)*(df - 2) * (df - 4));

} // variance

template <class RealType, class Policy>
inline RealType mode(const inverse_chi_squared_distribution<RealType, Policy>& dist)
{
   RealType df = dist.degrees_of_freedom();
   RealType scale = dist.scale();
   static const char* function = "boost::math::mode(const inverse_chi_squared_distribution<%1%>&)";
   if(df < 0)
      return policies::raise_domain_error<RealType>(
         function,
         "inverse Chi-Squared distribution only has a mode for degrees of freedom >= 0, but got degrees of freedom = %1%.",
         df, Policy());
   return (df * scale) / (df + 2);
}

//template <class RealType, class Policy>
//inline RealType median(const inverse_chi_squared_distribution<RealType, Policy>& dist)
//{ // Median is given by Quantile[dist, 1/2]
//   RealType df = dist.degrees_of_freedom();
//   if(df <= 1)
//      return tools::domain_error<RealType>(
//         BOOST_CURRENT_FUNCTION,
//         "The inverse_Chi-Squared distribution only has a median for degrees of freedom >= 0, but got degrees of freedom = %1%.",
//         df);
//   return df;
//}
// Now implemented via quantile(half) in derived accessors.

template <class RealType, class Policy>
inline RealType skewness(const inverse_chi_squared_distribution<RealType, Policy>& dist)
{
   BOOST_MATH_STD_USING // For ADL
   RealType df = dist.degrees_of_freedom();
   static const char* function = "boost::math::mode(const inverse_chi_squared_distribution<%1%>&)";
   if(df <= 6)
      return policies::raise_domain_error<RealType>(
         function,
         "inverse Chi-Squared distribution only has a skewness for degrees of freedom > 6, but got degrees of freedom = %1%.",
         df, Policy());

   return 4 * sqrt (2 * (df - 4)) / (df - 6);  // Not a function of scale.
}

template <class RealType, class Policy>
inline RealType kurtosis(const inverse_chi_squared_distribution<RealType, Policy>& dist)
{
   RealType df = dist.degrees_of_freedom();
   static const char* function = "boost::math::mode(const inverse_chi_squared_distribution<%1%>&)";
   if(df <= 8)
      return policies::raise_domain_error<RealType>(
         function,
         "inverse Chi-Squared distribution only has a kurtosis for degrees of freedom > 8, but got degrees of freedom = %1%.",
         df, Policy());

   return kurtosis_excess(dist) + 3;
}

template <class RealType, class Policy>
inline RealType kurtosis_excess(const inverse_chi_squared_distribution<RealType, Policy>& dist)
{
   RealType df = dist.degrees_of_freedom();
   static const char* function = "boost::math::mode(const inverse_chi_squared_distribution<%1%>&)";
   if(df <= 8)
      return policies::raise_domain_error<RealType>(
         function,
         "inverse Chi-Squared distribution only has a kurtosis excess for degrees of freedom > 8, but got degrees of freedom = %1%.",
         df, Policy());

   return 12 * (5 * df - 22) / ((df - 6 )*(df - 8));  // Not a function of scale.
}

/*
//
// Parameter estimation comes last:
//
namespace detail
{

template <class RealType, class Policy>
struct df_estimator
{
   df_estimator(RealType a, RealType b, RealType variance, RealType delta)
      : alpha(a), beta(b), ratio(delta/variance) {}

   RealType operator()(const RealType& df)
   {
      if(df <= tools::min_value<RealType>())
         return 1;
      inverse_chi_squared_distribution<RealType, Policy> cs(df);

      RealType result;
      if(ratio > 0)
      {
         RealType r = 1 + ratio;
         result = cdf(cs, quantile(complement(cs, alpha)) / r) - beta;
      }
      else
      {
         RealType r = 1 + ratio;
         result = cdf(complement(cs, quantile(cs, alpha) / r)) - beta;
      }
      return result;
   }
private:
   RealType alpha, beta, ratio;
};

} // namespace detail

template <class RealType, class Policy>
RealType inverse_chi_squared_distribution<RealType, Policy>::find_degrees_of_freedom(
   RealType difference_from_variance,
   RealType alpha,
   RealType beta,
   RealType variance,
   RealType hint)
{
   static const char* function = "boost::math::inverse_chi_squared_distribution<%1%>::find_degrees_of_freedom(%1%,%1%,%1%,%1%,%1%)";
   // Check for domain errors:
   RealType error_result;
   if(false == detail::check_probability(
         function, alpha, &error_result, Policy())
         && detail::check_probability(function, beta, &error_result, Policy()))
      return error_result;

   if(hint <= 0)
      hint = 1;

   detail::df_estimator<RealType, Policy> f(alpha, beta, variance, difference_from_variance);
   tools::eps_tolerance<RealType> tol(policies::digits<RealType, Policy>());
   boost::uintmax_t max_iter = policies::get_max_root_iterations<Policy>();
   std::pair<RealType, RealType> r = tools::bracket_and_solve_root(f, hint, RealType(2), false, tol, max_iter, Policy());
   RealType result = r.first + (r.second - r.first) / 2;
   if(max_iter >= policies::get_max_root_iterations<Policy>())
   {
      policies::raise_evaluation_error<RealType>(function, "Unable to locate solution in a reasonable time:"
         " either there is no answer to how many degrees of freedom are required"
         " or the answer is infinite.  Current best guess is %1%", result, Policy());
   }
   return result;
}

*/

} // namespace math
} // namespace boost

// This include must be at the end, *after* the accessors
// for this distribution have been defined, in order to
// keep compilers that support two-phase lookup happy.
#include <boost/math/distributions/detail/derived_accessors.hpp>

#endif // BOOST_MATH_DISTRIBUTIONS_INVERSE_CHI_SQUARED_HPP
