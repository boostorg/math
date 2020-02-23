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
#include <boost/math/distributions/complement.hpp>
#include <boost/math/distributions/detail/common_error_handling.hpp>
#include <boost/math/special_functions/bessel.hpp> // for besseli0 and besseli1
#include <boost/math/special_functions/erf.hpp> // for erf
#include <boost/math/tools/roots.hpp>

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

  RealType mean() const
  { // alias for location.
    return m_mean;
  }

  RealType concentration() const
  { // alias for scale.
    return m_concentration;
  }

  // Synonyms, provided to allow generic use of find_location and find_scale.
  RealType location() const
  { // location.
    return m_mean;
  }
  RealType scale() const
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
inline const std::pair<RealType, RealType> support(const von_mises_distribution<RealType, Policy>& dist)
{ // This is range values for random variable x where cdf rises from 0 to 1, and outside it, the pdf is zero.
  return std::pair<RealType, RealType>(dist.mean() - constants::pi<RealType>(),
                                       dist.mean() + constants::pi<RealType>()); //  [µ-π, µ+π)
}

#ifdef BOOST_MSVC
#pragma warning(pop)
#endif

namespace detail {
// float version of pdf_impl
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
    // Polynomial coefficients from boost/math/special_functions/detail/bessel_i0.hpp
    // > Max error in interpolated form: 5.195e-08
    // > Max Error found at float precision = Poly: 8.502534e-08
    static const float P[] = {
        3.98942651588301770e-01f,
        4.98327234176892844e-02f,
        2.91866904423115499e-02f,
        1.35614940793742178e-02f,
        1.31409251787866793e-01f
    };
    RealType result = exp(conc * (cos(x - mean) - 1.f));
    result /= boost::math::tools::evaluate_polynomial(P, RealType(1.f / conc)) / sqrt(conc)
              * boost::math::constants::two_pi<RealType>();
    return result;
  }
  else
  {
    // Polynomial coefficients from boost/math/special_functions/detail/bessel_i0.hpp
    // > Max error in interpolated form: 1.782e-09
    // > Max Error found at float precision = Poly: 6.473568e-08
    static const float P[] = {
        3.98942280401432677e-01f,
        4.98677850501790847e-02f,
        2.80506290907257351e-02f,
        2.92194053028393074e-02f,
        4.47422143699726895e-02f,
    };
    RealType result = exp(conc * (cos(x - mean) - 1.f));
    result /= boost::math::tools::evaluate_polynomial(P, RealType(1.f / conc)) / sqrt(conc)
              * boost::math::constants::two_pi<RealType>();
    return result;
  }
}

// double version of pdf_impl
template <typename RealType, class Policy>
inline RealType pdf_impl(const von_mises_distribution<RealType, Policy>& dist, RealType x,
                         const boost::integral_constant<int, 53>&)
{
  RealType mean = dist.mean();
  RealType conc = dist.concentration();

  BOOST_MATH_STD_USING
  if(conc < 7.75)
  {
    RealType bessel_i0 = cyl_bessel_i(0, conc, Policy());
    return exp(conc * cos(x - mean)) / bessel_i0 / boost::math::constants::two_pi<RealType>();
  }
  else if(conc < 500)
  {
    // Polynomial coefficients from boost/math/special_functions/detail/bessel_i0.hpp
    // > Max error in interpolated form : 1.685e-16
    // > Max Error found at double precision = Poly : 2.575063e-16 Cheb : 2.247615e+00
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
    RealType result = exp(conc * (cos(x - mean) - 1.0));
    result /= boost::math::tools::evaluate_polynomial(P, RealType(1.0 / conc)) / sqrt(conc)
              * boost::math::constants::two_pi<RealType>();
    return result;
  }
  else
  {
    // Polynomial coefficients from boost/math/special_functions/detail/bessel_i0.hpp
    // > Max error in interpolated form : 2.437e-18
    // > Max Error found at double precision = Poly : 1.216719e-16
    static const double P[] = {
        3.98942280401432905e-01,
        4.98677850491434560e-02,
        2.80506308916506102e-02,
        2.92179096853915176e-02,
        4.53371208762579442e-02
    };
    RealType result = exp(conc * (cos(x - mean) - 1.0));
    result /= boost::math::tools::evaluate_polynomial(P, RealType(1.0 / conc)) / sqrt(conc)
              * boost::math::constants::two_pi<RealType>();
    return result;
  }
}

// long double version of pdf_impl
template <typename RealType, class Policy>
inline RealType pdf_impl(const von_mises_distribution<RealType, Policy>& dist, RealType x,
                         const boost::integral_constant<int, 64>&)
{
  RealType mean = dist.mean();
  RealType conc = dist.concentration();

  BOOST_MATH_STD_USING
  if (conc < 15)
  {
    RealType bessel_i0 = cyl_bessel_i(0, conc, Policy());
    return exp(conc * cos(x - mean)) / bessel_i0 / boost::math::constants::two_pi<RealType>();
  }
  else if (x < 50)
  {
    // Max error in interpolated form: 1.035e-21
    // Max Error found at float80 precision = Poly: 1.885872e-21
    static const RealType Y = 4.011702537536621093750e-01f;
    static const RealType P[] = {
        BOOST_MATH_BIG_CONSTANT(RealType, 64, -2.227973351806078464328e-03),
        BOOST_MATH_BIG_CONSTANT(RealType, 64, 4.986778486088017419036e-02),
        BOOST_MATH_BIG_CONSTANT(RealType, 64, 2.805066823812285310011e-02),
        BOOST_MATH_BIG_CONSTANT(RealType, 64, 2.921443721160964964623e-02),
        BOOST_MATH_BIG_CONSTANT(RealType, 64, 4.517504941996594744052e-02),
        BOOST_MATH_BIG_CONSTANT(RealType, 64, 6.316922639868793684401e-02),
        BOOST_MATH_BIG_CONSTANT(RealType, 64, 1.535891099168810015433e+00),
        BOOST_MATH_BIG_CONSTANT(RealType, 64, -4.706078229522448308087e+01),
        BOOST_MATH_BIG_CONSTANT(RealType, 64, 1.351015763079160914632e+03),
        BOOST_MATH_BIG_CONSTANT(RealType, 64, -2.948809013999277355098e+04),
        BOOST_MATH_BIG_CONSTANT(RealType, 64, 4.967598958582595361757e+05),
        BOOST_MATH_BIG_CONSTANT(RealType, 64, -6.346924657995383019558e+06),
        BOOST_MATH_BIG_CONSTANT(RealType, 64, 5.998794574259956613472e+07),
        BOOST_MATH_BIG_CONSTANT(RealType, 64, -4.016371355801690142095e+08),
        BOOST_MATH_BIG_CONSTANT(RealType, 64, 1.768791455631826490838e+09),
        BOOST_MATH_BIG_CONSTANT(RealType, 64, -4.441995678177349895640e+09),
        BOOST_MATH_BIG_CONSTANT(RealType, 64, 4.482292669974971387738e+09)
    };
    RealType result = exp(conc * (cos(x - mean) - 1.0));
    result /= (boost::math::tools::evaluate_polynomial(P, RealType(1.0 / conc)) + Y) / sqrt(conc);
              * boost::math::constants::two_pi<RealType>();
    return result;
  }
  else
  {
    // Bessel I0 over[50, INF]
    // Max error in interpolated form : 5.587e-20
    // Max Error found at float80 precision = Poly : 8.776852e-20
    static const RealType P[] = {
        BOOST_MATH_BIG_CONSTANT(RealType, 64, +3.98942280401432677955074061e-01),
        BOOST_MATH_BIG_CONSTANT(RealType, 64, +4.98677850501789875615574058e-02),
        BOOST_MATH_BIG_CONSTANT(RealType, 64, +2.80506290908675604202206833e-02),
        BOOST_MATH_BIG_CONSTANT(RealType, 64, +2.92194052159035901631494784e-02),
        BOOST_MATH_BIG_CONSTANT(RealType, 64, +4.47422430732256364094681137e-02),
        BOOST_MATH_BIG_CONSTANT(RealType, 64, +9.05971614435738691235525172e-02),
        BOOST_MATH_BIG_CONSTANT(RealType, 64, +2.29180522595459823234266708e-01),
        BOOST_MATH_BIG_CONSTANT(RealType, 64, +6.15122547776140254569073131e-01),
        BOOST_MATH_BIG_CONSTANT(RealType, 64, +7.48491812136365376477357324e+00),
        BOOST_MATH_BIG_CONSTANT(RealType, 64, -2.45569740166506688169730713e+02),
        BOOST_MATH_BIG_CONSTANT(RealType, 64, +9.66857566379480730407063170e+03),
        BOOST_MATH_BIG_CONSTANT(RealType, 64, -2.71924083955641197750323901e+05),
        BOOST_MATH_BIG_CONSTANT(RealType, 64, +5.74276685704579268845870586e+06),
        BOOST_MATH_BIG_CONSTANT(RealType, 64, -8.89753803265734681907148778e+07),
        BOOST_MATH_BIG_CONSTANT(RealType, 64, +9.82590905134996782086242180e+08),
        BOOST_MATH_BIG_CONSTANT(RealType, 64, -7.30623197145529889358596301e+09),
        BOOST_MATH_BIG_CONSTANT(RealType, 64, +3.27310000726207055200805893e+10),
        BOOST_MATH_BIG_CONSTANT(RealType, 64, -6.64365417189215599168817064e+10)
    };
    RealType result = exp(conc * (cos(x - mean) - 1.0));
    result /= boost::math::tools::evaluate_polynomial(P, RealType(1.0 / conc)) / sqrt(conc)
              * boost::math::constants::two_pi<RealType>();
    return result;
  }
}
} // namespace detail

template <class RealType, class Policy>
inline RealType pdf(const von_mises_distribution<RealType, Policy>& dist, const RealType& x)
{
  BOOST_MATH_STD_USING  // for ADL of std functions

  RealType conc = dist.concentration();
  RealType mean = dist.mean();
  
  static const char* function = "boost::math::pdf(const von_mises_distribution<%1%>&, %1%)";
  
  RealType result = 0;
  if (false == detail::check_positive_x(function, conc, &result, Policy()))
    return result;
  if (false == detail::check_angle(function, mean, &result, Policy()))
    return result;
  if (false == detail::check_angle(function, x - mean, &result, Policy()))
    return result;
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

  return detail::pdf_impl(dist, x, tag_type());
} // pdf

namespace detail {

template <class RealType, class Policy>
inline RealType cdf_impl(const von_mises_distribution<RealType, Policy>& dist, const RealType& x) 
{
  BOOST_MATH_STD_USING  // for ADL of std functions
  
  RealType conc = dist.concentration();
  RealType mean = dist.mean();
  
  RealType u = x - mean;
  if (u <= -boost::math::constants::pi<RealType>())
    return 0;
  if (u >= +boost::math::constants::pi<RealType>())
    return 1;  

  // We use the Fortran algorithm designed by Geoffrey W. Hill in
  // "Algorithm 518: Incomplete Bessel Function I0. The Von Mises Distribution", 1977, ACM
  // DOI: 10.1145/355744.355753
  RealType result = 0;
  if (conc > RealType{10.5}) {
    auto c = 24.0 * conc;
    auto v = c - 56.0;
    auto r = sqrt((54 / (347 / v + 26 - c) - 6 + c) / 12);    
    auto z = sin(u / 2) * r;
    auto s = z * z * 2;
    v = v - s + 3;
    auto y = (c - s - s - 16) / 3;
    y = ((s + 1.75) * s + 83.5) / v - y;
    result = boost::math::erf(z - s / (y * y) * z) / 2 + 0.5;
  }
  else {
    RealType v = 0;
    if(conc > 0) {
      int iterations = static_cast<int>(ceil(conc * 0.8 - 8 / (conc + 1) + 12));
      RealType r = 0;
      RealType z = 2 / conc;
      for (int j = iterations - 1; j > 0; --j) {        
        RealType sj = sin(j * u);
        r = 1 / (j * z + r);
        v = (sj / j + v) * r;
      }
    }
    result = (x - mean + boost::math::constants::pi<RealType>()) / 2;
    result = (result + v) / boost::math::constants::pi<RealType>();
  }
  
  return result <= 0 ? 0 : (1 <= result ? 1 : result);
}
} // namespace detail

template <class RealType, class Policy>
inline RealType cdf(const von_mises_distribution<RealType, Policy>& dist, const RealType& x)
{
  RealType conc = dist.concentration();
  RealType mean = dist.mean();
  static const char* function = "boost::math::cdf(const von_mises_distribution<%1%>&, %1%)";
  RealType result = 0;
  if (false == detail::check_positive_x(function, conc, &result, Policy()))
    return result;
  if (false == detail::check_angle(function, mean, &result, Policy()))
    return result;
  if (false == detail::check_angle(function, x - mean, &result, Policy()))
    return result;
   
  return detail::cdf_impl(dist, x);
} // cdf

template <class RealType, class Policy>
inline RealType quantile(const von_mises_distribution<RealType, Policy>& dist, const RealType& p)
{
  BOOST_MATH_STD_USING  // for ADL of std functions

  RealType conc = dist.concentration();
  RealType mean = dist.mean();
  static const char* function = "boost::math::quantile(const von_mises_distribution<%1%>&, %1%)";

  RealType result = 0;
  if (false == detail::check_positive_x(function, conc, &result, Policy()))
    return result;
  if (false == detail::check_angle(function, mean, &result, Policy()))
    return result;
  if (false == detail::check_probability(function, p, &result, Policy()))
    return result;
    
  if (p <= 0)
    return -boost::math::constants::pi<RealType>();
  if (p >= 1)
    return +boost::math::constants::pi<RealType>();
  
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
     
  auto step_func = [&](RealType x) {
      return std::make_pair(detail::cdf_impl(dist, x) - p,            // f(x)
                            detail::pdf_impl(dist, x, tag_type()));   // f'(x)
  };
  RealType lower = mean - boost::math::constants::pi<RealType>();
  RealType upper = mean + boost::math::constants::pi<RealType>();
  RealType zero = boost::math::tools::newton_raphson_iterate(
      step_func, mean, lower, upper, 15 /* digits */);
  
  return zero;
} // quantile

template <class RealType, class Policy>
inline RealType cdf(const complemented2_type<von_mises_distribution<RealType, Policy>, RealType>& c)
{
  RealType conc = c.dist.concentration();
  RealType mean = c.dist.mean();
  RealType x = c.param;
  static const char* function = "boost::math::cdf(const complement(von_mises_distribution<%1%>&), %1%)";

  RealType result = 0;
  if (false == detail::check_positive_x(function, conc, &result, Policy()))
    return result;
  if (false == detail::check_angle(function, mean, &result, Policy()))
    return result;
  if (false == detail::check_angle(function, x - mean, &result, Policy()))
    return result;

  return detail::cdf_impl(c.dist, 2 * mean - x);
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


