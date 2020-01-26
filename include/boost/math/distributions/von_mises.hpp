//  Copyright John Maddock 2006, 2007.
//  Copyright Paul A. Bristow 2006, 2007.

//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_STATS_VON_MISES_HPP
#define BOOST_STATS_VON_MISES_HPP

// https://en.wikipedia.org/wiki/Von_Mises_distribution
// From MathWorld--A Wolfram Web Resource.
// http://mathworld.wolfram.com/VonMisesDistribution.html

#include <boost/math/distributions/fwd.hpp>
#include <boost/math/special_functions/bessel.hpp> 
#include <boost/math/special_functions/erf.hpp> // for erf/erfc.
#include <boost/math/distributions/complement.hpp>
#include <boost/math/distributions/detail/common_error_handling.hpp>

#include <utility>

namespace boost{ namespace math{

template <class RealType = double, class Policy = policies::policy<> >
class von_mises_distribution
{
public:
   typedef RealType value_type;
   typedef Policy policy_type;

   von_mises_distribution(RealType l_mean = 0, RealType concentration = 1)
      : m_mean(l_mean), m_concentration(concentration)
   { // Default is a 'standard' von Mises distribution vM01.
     static const char* function = "boost::math::von_mises_distribution<%1%>::von_mises_distribution";

     RealType result;
     detail::check_positive_x(function, concentration, &result, Policy());
     detail::check_angle(function, l_mean, &result, Policy());
   }

   RealType mean()const
   { // alias for location.
      return m_mean;
   }

   RealType concentration() const
   { // alias for scale.
      return m_concentration;
   }

   // Synonyms, provided to allow generic use of find_location and find_scale.
   RealType location()const
   { // location.
      return m_mean;
   }
   RealType scale()const
   { // scale.
      return m_concentration;
   }

private:
   //
   // Data members:
   //
   RealType m_mean;  // distribution mean or location.
   RealType m_concentration;    // distribution standard deviation or scale.
}; // class von_mises_distribution

typedef von_mises_distribution<double> von_mises;

#ifdef BOOST_MSVC
#pragma warning(push)
#pragma warning(disable:4127)
#endif

template <class RealType, class Policy>
inline const std::pair<RealType, RealType> range(const von_mises_distribution<RealType, Policy>& /*dist*/)
{ // Range of permissible values for random variable x.
  if (std::numeric_limits<RealType>::has_infinity)
  { 
     return std::pair<RealType, RealType>(-std::numeric_limits<RealType>::infinity(),
                                          +std::numeric_limits<RealType>::infinity()); // - to + infinity.
  }
  else
  { // Can only use max_value.
    using boost::math::tools::max_value;
    return std::pair<RealType, RealType>(-max_value<RealType>(),
                                         +max_value<RealType>()); // - to + max value.
  }
}

template <class RealType, class Policy>
inline const std::pair<RealType, RealType> support(const von_mises_distribution<RealType, Policy>& /*dist*/)
{ // This is range values for random variable x where cdf rises from 0 to 1, and outside it, the pdf is zero.
   return std::pair<RealType, RealType>(-constants::pi<RealType>(),
                                        +constants::pi<RealType>()); // -pi to +pi
}

#ifdef BOOST_MSVC
#pragma warning(pop)
#endif

// float version of log_bessel_i0
template <class RealType, class Policy>
inline RealType pdf_impl(const von_mises_distribution<RealType, Policy>& dist, RealType x,
                         const boost::integral_constant<int, 24>&)
{
  RealType mean = dist.mean();
  RealType conc = dist.concentration();

  BOOST_MATH_STD_USING
  if(conc < 7.75)
  {
    RealType bessel_i0 = cyl_bessel_i(0, conc, Policy());
    return exp(conc * cos(x - mean)) / bessel_i0 / boost::math::constants::two_pi<RealType>();
  }
  else if(conc < 50)
  {
    // Max error in interpolated form: 5.195e-08
    // Max Error found at float precision = Poly: 8.502534e-08
    static const float P[] = {
        3.98942651588301770e-01f,
        4.98327234176892844e-02f,
        2.91866904423115499e-02f,
        1.35614940793742178e-02f,
        1.31409251787866793e-01f
    };
    RealType exponent = conc * (cos(x - mean) - 1.f);
    exponent -= std::log(boost::math::tools::evaluate_polynomial(P, RealType(1 / conc)) / sqrt(conc));
    return exp(exponent) / boost::math::constants::two_pi<RealType>();
  }
  else
  {
    // Max error in interpolated form: 1.782e-09
    // Max Error found at float precision = Poly: 6.473568e-08
    static const float P[] = {
        3.98942280401432677e-01f,
        4.98677850501790847e-02f,
        2.80506290907257351e-02f,
        2.92194053028393074e-02f,
        4.47422143699726895e-02f,
    };
    RealType exponent = conc * (cos(x - mean) - 1.f);
    exponent -= std::log(boost::math::tools::evaluate_polynomial(P, RealType(1 / conc)) / sqrt(conc));
    return exp(exponent) / boost::math::constants::two_pi<RealType>();
  }
}

// double version of log_bessel_i0
template <typename RealType>
inline RealType log_bessel_i0_impl(RealType x, const boost::integral_constant<int, 53>&)
{
   BOOST_MATH_STD_USING
   if(x < 7.75)
   {
      static const double P[] = {
         1.00000000000000000e+00,
         2.49999999999999909e-01,
         2.77777777777782257e-02,
         1.73611111111023792e-03,
         6.94444444453352521e-05,
         1.92901234513219920e-06,
         3.93675991102510739e-08,
         6.15118672704439289e-10,
         7.59407002058973446e-12,
         7.59389793369836367e-14,
         6.27767773636292611e-16,
         4.34709704153272287e-18,
         2.63417742690109154e-20,
         1.13943037744822825e-22,
         9.07926920085624812e-25
      };
      RealType a = x * x / 4;
      return std::log(a * boost::math::tools::evaluate_polynomial(P, a) + 1);
   }
   else if(x < 500)
   {
      static const double P[] = {
         3.98942280401425088e-01,
         4.98677850604961985e-02,
         2.80506233928312623e-02,
         2.92211225166047873e-02,
         4.44207299493659561e-02,
         1.30970574605856719e-01,
         -3.35052280231727022e+00,
         2.33025711583514727e+02,
         -1.13366350697172355e+04,
         4.24057674317867331e+05,
         -1.23157028595698731e+07,
         2.80231938155267516e+08,
         -5.01883999713777929e+09,
         7.08029243015109113e+10,
         -7.84261082124811106e+11,
         6.76825737854096565e+12,
         -4.49034849696138065e+13,
         2.24155239966958995e+14,
         -8.13426467865659318e+14,
         2.02391097391687777e+15,
         -3.08675715295370878e+15,
         2.17587543863819074e+15
      };
      return x + std::log(boost::math::tools::evaluate_polynomial(P, RealType(1 / x)) / sqrt(x));
   }
   else
   {
      static const double P[] = {
         3.98942280401432905e-01,
         4.98677850491434560e-02,
         2.80506308916506102e-02,
         2.92179096853915176e-02,
         4.53371208762579442e-02
      };
      RealType ex = x / 2;
      RealType result = ex + std::log(boost::math::tools::evaluate_polynomial(P, RealType(1 / x)) / sqrt(x));
      result += ex;
      return result;
   }
}

template <class RealType>
inline RealType log_bessel_i0(RealType x)
{
  typedef boost::integral_constant<int,
      ((std::numeric_limits<RealType>::digits == 0) || (std::numeric_limits<RealType>::radix != 2)) ?
        0 :
        std::numeric_limits<RealType>::digits <= 24 ?
        24 :
        std::numeric_limits<RealType>::digits <= 53 ?
        53 :
        std::numeric_limits<RealType>::digits <= 64 ?
        64 :
        std::numeric_limits<RealType>::digits <= 113 ?
        113 : -1
      > tag_type;

  return log_bessel_i0_impl(x, tag_type());
}

template <class RealType, class Policy>
inline RealType pdf(const von_mises_distribution<RealType, Policy>& dist, const RealType& x)
{
   BOOST_MATH_STD_USING  // for ADL of std functions

   RealType conc = dist.concentration();
   RealType mean = dist.mean();

   static const char* function = "boost::math::pdf(const von_mises_distribution<%1%>&, %1%)";

   RealType result = 0;
   if(false == detail::check_positive_x(function, conc, &result, Policy()))
   {
      return result;
   }
   if(false == detail::check_angle(function, mean, &result, Policy()))
   {
      return result;
   }
   if(false == detail::check_x(function, x, &result, Policy()))
   {
     return result;
   }
   // Below produces MSVC 4127 warnings, so the above used instead.
   //if(std::numeric_limits<RealType>::has_infinity && abs(x) == std::numeric_limits<RealType>::infinity())
   //{ // pdf + and - infinity is zero.
   //  return 0;
   //}
   typedef boost::integral_constant<int,
       ((std::numeric_limits<RealType>::digits == 0) || (std::numeric_limits<RealType>::radix != 2)) ?
         0 :
         std::numeric_limits<RealType>::digits <= 24 ?
         24 :
         std::numeric_limits<RealType>::digits <= 53 ?
         53 :
         std::numeric_limits<RealType>::digits <= 64 ?
         64 :
         std::numeric_limits<RealType>::digits <= 113 ?
         113 : -1
       > tag_type;

   return pdf_impl(dist, x, tag_type());
} // pdf

template <class RealType, class Policy>
inline RealType cdf(const von_mises_distribution<RealType, Policy>& dist, const RealType& x)
{
   BOOST_MATH_STD_USING  // for ADL of std functions

   RealType conc = dist.concentration();
   RealType mean = dist.mean();
   static const char* function = "boost::math::cdf(const von_mises_distribution<%1%>&, %1%)";
   RealType result = 0;
   if(false == detail::check_positive_x(function, conc, &result, Policy()))
   {
      return result;
   }
   if(false == detail::check_angle(function, mean, &result, Policy()))
   {
      return result;
   }
   if(false == detail::check_angle(function, x, &result, Policy()))
   {
     return result;
   }
   
   RealType diff = (x - mean) / (conc * constants::root_two<RealType>());
   result = boost::math::erfc(-diff, Policy()) / 2;
   return result;
} // cdf

// template <class RealType, class Policy>
// inline RealType quantile(const von_mises_distribution<RealType, Policy>& dist, const RealType& p)
// {
   // BOOST_MATH_STD_USING  // for ADL of std functions

   // RealType conc = dist.concentration();
   // RealType mean = dist.mean();
   // static const char* function = "boost::math::quantile(const von_mises_distribution<%1%>&, %1%)";

   // RealType result = 0;
   // if(false == detail::check_positive_x(function, conc, &result, Policy()))
      // return result;
   // if(false == detail::check_angle(function, mean, &result, Policy()))
      // return result;
   // if(false == detail::check_probability(function, p, &result, Policy()))
      // return result;

   // result = boost::math::erfc_inv(2 * p, Policy());
   // result = -result;
   // result *= sd * constants::root_two<RealType>();
   // result += mean;
   // return result;
// } // quantile

template <class RealType, class Policy>
inline RealType cdf(const complemented2_type<von_mises_distribution<RealType, Policy>, RealType>& c)
{
   BOOST_MATH_STD_USING  // for ADL of std functions

   RealType conc = c.dist.concentration();
   RealType mean = c.dist.mean();
   RealType x = c.param;
   static const char* function = "boost::math::cdf(const complement(von_mises_distribution<%1%>&), %1%)";

   RealType result = 0;
   if(false == detail::check_positive_x(function, conc, &result, Policy()))
      return result;
   if(false == detail::check_angle(function, mean, &result, Policy()))
      return result;
   if(false == detail::check_angle(function, x, &result, Policy()))
      return result;

   RealType diff = (x - mean) / (conc * constants::root_two<RealType>());
   result = boost::math::erfc(diff, Policy()) / 2;
   return result;
} // cdf complement

// template <class RealType, class Policy>
// inline RealType quantile(const complemented2_type<von_mises_distribution<RealType, Policy>, RealType>& c)
// {
   // BOOST_MATH_STD_USING  // for ADL of std functions

   // RealType conc = c.dist.concentration();
   // RealType mean = c.dist.mean();
   // static const char* function = "boost::math::quantile(const complement(von_mises_distribution<%1%>&), %1%)";
   // RealType result = 0;
   // if(false == detail::check_positive_x(function, conc, &result, Policy()))
      // return result;
   // if(false == detail::check_angle(function, mean, &result, Policy()))
      // return result;
   // RealType q = c.param;
   // if(false == detail::check_probability(function, q, &result, Policy()))
      // return result;
   // result = boost::math::erfc_inv(2 * q, Policy());
   // result *= conc * constants::root_two<RealType>();
   // result += mean;
   // return result;
// } // quantile

template <class RealType, class Policy>
inline RealType mean(const von_mises_distribution<RealType, Policy>& dist)
{
   return dist.mean();
}

template <class RealType, class Policy>
inline RealType standard_deviation(const von_mises_distribution<RealType, Policy>& dist)
{
   return dist.concentration();
}

template <class RealType, class Policy>
inline RealType mode(const von_mises_distribution<RealType, Policy>& dist)
{
   return dist.mean();
}

template <class RealType, class Policy>
inline RealType median(const von_mises_distribution<RealType, Policy>& dist)
{
   return dist.mean();
}

template <class RealType, class Policy>
inline RealType skewness(const von_mises_distribution<RealType, Policy>& /*dist*/)
{
   return 0;
}

// template <class RealType, class Policy>
// inline RealType kurtosis(const von_mises_distribution<RealType, Policy>& /*dist*/)
// {
//   return 3;
// }

// template <class RealType, class Policy>
// inline RealType kurtosis_excess(const von_mises_distribution<RealType, Policy>& /*dist*/)
// {
   // return 0;
// } 

template <class RealType, class Policy>
inline RealType entropy(const von_mises_distribution<RealType, Policy> & dist)
{
   using std::log;
   RealType arg = constants::two_pi<RealType>() * cyl_bessel_i(0, dist.concentration(), Policy());
   return log(arg) - dist.concentration() * cyl_bessel_i(1, dist.concentration(), Policy()) / cyl_bessel_i(0, dist.concentration(), Policy());
}

} // namespace math
} // namespace boost

// This include must be at the end, *after* the accessors
// for this distribution have been defined, in order to
// keep compilers that support two-phase lookup happy.
#include <boost/math/distributions/detail/derived_accessors.hpp>

#endif // BOOST_STATS_VON_MISES_HPP


