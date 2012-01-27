// (C) Benjamin Sobotta 2012

//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_STATS_SKEW_NORMAL_HPP
#define BOOST_STATS_SKEW_NORMAL_HPP

// http://en.wikipedia.org/wiki/Skew_normal_distribution
// http://azzalini.stat.unipd.it/SN/
// Also:
// Azzalini, A. (1985). "A class of distributions which includes the normal ones".
// Scand. J. Statist. 12: 171-178.

#include <boost/math/distributions/fwd.hpp> // TODO add skew_normal distribution to fwd.hpp!
#include <boost/math/special_functions/owens_t.hpp> // Owen's T function
#include <boost/math/distributions/complement.hpp>
#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/detail/common_error_handling.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/math/tools/tuple.hpp>
#include <boost/math/tools/roots.hpp> // Newton-Raphson
#include <boost/math/distributions/detail/generic_mode.hpp> // pdf max finder.

#include <utility>

namespace boost{ namespace math{

  namespace detail
  {
    template <class RealType, class Policy>
    inline bool check_skew_normal_shape(
      const char* function,
      RealType location,
      RealType* result,
      const Policy& pol)
    {
      if(!(boost::math::isfinite)(location))
      {
        *result =
          policies::raise_domain_error<RealType>(function,
          "Shape parameter is %1%, but must be finite!",
          location, pol);
        return false;
      }
      return true;
    }

  } // namespace detail

  template <class RealType = double, class Policy = policies::policy<> >
  class skew_normal_distribution
  {
  public:
    typedef RealType value_type;
    typedef Policy policy_type;

    skew_normal_distribution(RealType location = 0, RealType scale = 1, RealType shape = 0)
      : location_(location), scale_(scale), shape_(shape)
    { // Default is a 'standard' normal distribution N01. (shape=0 results in the normal distribution with no skew)
      static const char* function = "boost::math::skew_normal_distribution<%1%>::skew_normal_distribution";

      RealType result;
      detail::check_scale(function, scale, &result, Policy());
      detail::check_location(function, location, &result, Policy());
      detail::check_skew_normal_shape(function, shape, &result, Policy());
    }

    RealType location()const
    { 
      return location_;
    }

    RealType scale()const
    { 
      return scale_;
    }

    RealType shape()const
    { 
      return shape_;
    }


  private:
    //
    // Data members:
    //
    RealType location_;  // distribution location.
    RealType scale_;    // distribution scale.
    RealType shape_;    // distribution shape.
  }; // class skew_normal_distribution

  typedef skew_normal_distribution<double> skew_normal;

  template <class RealType, class Policy>
  inline const std::pair<RealType, RealType> range(const skew_normal_distribution<RealType, Policy>& /*dist*/)
  { // Range of permissible values for random variable x.
    using boost::math::tools::max_value;
    return std::pair<RealType, RealType>(-max_value<RealType>(), max_value<RealType>()); // - to + max value.
  }

  template <class RealType, class Policy>
  inline const std::pair<RealType, RealType> support(const skew_normal_distribution<RealType, Policy>& /*dist*/)
  { // Range of supported values for random variable x.
    // This is range where cdf rises from 0 to 1, and outside it, the pdf is zero.

    using boost::math::tools::max_value;
    return std::pair<RealType, RealType>(-max_value<RealType>(),  max_value<RealType>()); // - to + max value.
  }

  template <class RealType, class Policy>
  inline RealType pdf(const skew_normal_distribution<RealType, Policy>& dist, const RealType& x)
  {
    const RealType scale = dist.scale();
    const RealType location = dist.location();
    const RealType shape = dist.shape();

    static const char* function = "boost::math::pdf(const skew_normal_distribution<%1%>&, %1%)";
    if((boost::math::isinf)(x))
    {
      return 0; // pdf + and - infinity is zero.
    }
    // Below produces MSVC 4127 warnings, so the above used instead.
    //if(std::numeric_limits<RealType>::has_infinity && abs(x) == std::numeric_limits<RealType>::infinity())
    //{ // pdf + and - infinity is zero.
    //  return 0;
    //}

    RealType result = 0;
    if(false == detail::check_scale(function, scale, &result, Policy()))
    {
      return result;
    }
    if(false == detail::check_location(function, location, &result, Policy()))
    {
      return result;
    }
    if(false == detail::check_skew_normal_shape(function, shape, &result, Policy()))
    {
      return result;
    }
    if(false == detail::check_x(function, x, &result, Policy()))
    {
      return result;
    }

    const RealType transformed_x = (x-location)/scale;

    normal_distribution<RealType, Policy> std_normal;

    result = pdf(std_normal, transformed_x) * cdf(std_normal, shape*transformed_x) * 2 / scale;

    return result;
  } // pdf

  template <class RealType, class Policy>
  inline RealType cdf(const skew_normal_distribution<RealType, Policy>& dist, const RealType& x)
  {
    const RealType scale = dist.scale();
    const RealType location = dist.location();
    const RealType shape = dist.shape();

    static const char* function = "boost::math::cdf(const skew_normal_distribution<%1%>&, %1%)";
    RealType result = 0;
    if(false == detail::check_scale(function, scale, &result, Policy()))
    {
      return result;
    }
    if(false == detail::check_location(function, location, &result, Policy()))
    {
      return result;
    }
    if(false == detail::check_skew_normal_shape(function, shape, &result, Policy()))
    {
      return result;
    }
    if((boost::math::isinf)(x))
    {
      if(x < 0) return 0; // -infinity
      return 1; // + infinity
    }
    // These produce MSVC 4127 warnings, so the above used instead.
    //if(std::numeric_limits<RealType>::has_infinity && x == std::numeric_limits<RealType>::infinity())
    //{ // cdf +infinity is unity.
    //  return 1;
    //}
    //if(std::numeric_limits<RealType>::has_infinity && x == -std::numeric_limits<RealType>::infinity())
    //{ // cdf -infinity is zero.
    //  return 0;
    //}
    if(false == detail::check_x(function, x, &result, Policy()))
    {
      return result;
    }

    const RealType transformed_x = (x-location)/scale;

    normal_distribution<RealType, Policy> std_normal;

    result = cdf(std_normal, transformed_x) - owens_t(transformed_x, shape)*2;

    return result;
  } // cdf

  template <class RealType, class Policy>
  inline RealType cdf(const complemented2_type<skew_normal_distribution<RealType, Policy>, RealType>& c)
  {
    const RealType scale = c.dist.scale();
    const RealType location = c.dist.location();
    const RealType shape = c.dist.shape();
    const RealType x = c.param;

    static const char* function = "boost::math::cdf(const complement(skew_normal_distribution<%1%>&), %1%)";

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
    RealType result = 0;
    if(false == detail::check_scale(function, scale, &result, Policy()))
      return result;
    if(false == detail::check_location(function, location, &result, Policy()))
      return result;
    if(false == detail::check_skew_normal_shape(function, shape, &result, Policy()))
      return result;
    if(false == detail::check_x(function, x, &result, Policy()))
      return result;

    const RealType transformed_x = (x-location)/scale;

    normal_distribution<RealType, Policy> std_normal;

    result = cdf(complement(std_normal, transformed_x)) + owens_t(transformed_x, shape)*2;
    return result;
  } // cdf complement

  template <class RealType, class Policy>
  inline RealType location(const skew_normal_distribution<RealType, Policy>& dist)
  {
    return dist.location();
  }

  template <class RealType, class Policy>
  inline RealType scale(const skew_normal_distribution<RealType, Policy>& dist)
  {
    return dist.scale();
  }

  template <class RealType, class Policy>
  inline RealType shape(const skew_normal_distribution<RealType, Policy>& dist)
  {
    return dist.shape();
  }

  template <class RealType, class Policy>
  inline RealType mean(const skew_normal_distribution<RealType, Policy>& dist)
  {
    BOOST_MATH_STD_USING  // for ADL of std functions

    using namespace boost::math::constants;
    // TODO replace by root_two_div_pi<RealType>() when in boost-trunk
    //static const RealType sqrt_two_over_pi = root_two_div_pi<RealType>(); // = ( 2 / pi ) ^ 0.5
    static const RealType sqrt_two_over_pi = 7.978845608028654e-001; // = ( 2 / pi ) ^ 0.5
    const RealType delta = dist.shape() / sqrt(1+dist.shape()*dist.shape());

    return dist.location() + dist.scale()*delta*sqrt_two_over_pi;
    //return dist.location() + dist.scale() * delta * root_two_div_pi<RealType>();
  }

  template <class RealType, class Policy>
  inline RealType variance(const skew_normal_distribution<RealType, Policy>& dist)
  {
    BOOST_MATH_STD_USING  // for ADL of std functions

    using namespace boost::math::constants;
    //static const RealType two_over_pi = two_div_pi<RealType>(); // = 2 / pi
    // TODO replace by two_div_pi<RealType>() when in boost-trunk
    static const RealType two_over_pi = 6.366197723675814e-001; // = 2 / pi
    const RealType delta = dist.shape() / sqrt(1+dist.shape()*dist.shape());
    RealType variance = dist.scale()*dist.scale()*(1-two_over_pi*delta*delta);
    //RealType variance = dist.scale()*dist.scale()*(1-two_div_pi<RealType>() *delta*delta);
    return variance;
  }

  /*
     TODO No closed expression, so use max of pdf.
  */

  template <class RealType, class Policy>
  inline RealType mode(const skew_normal_distribution<RealType, Policy>& dist)
  { // mode.
      static const char* function = "mode(skew_normal_distribution<%1%> const&)";
      RealType scale = dist.scale();
      RealType shape = dist.shape();
      RealType r;
      if(!detail::check_scale(
        function,
        scale, &r, Policy())
        ||
      !detail::check_finite(
        function,
        shape,
        &r,
        Policy()))
      return (RealType)r;

      BOOST_MATH_STD_USING

         //RealType m = v < 3 ? 0 : detail::mean(v, l, Policy());
         //RealType var = v < 4 ? 1 : detail::variance(v, l, Policy());

      RealType guess = (scale < 3) ? 0 : mean(dist);
      RealType step = (scale < 4) ? 0 : variance(dist); //

      cout << "guess= " << guess << ", step = " << step << endl;

      // TODO problem here in finding lower and upper bounds.
      // Guess looks reasonable but What should step be?

      return detail::generic_find_mode(
        dist, 
        guess, 
        function);
  } // mode

  template <class RealType, class Policy>
  inline RealType skewness(const skew_normal_distribution<RealType, Policy>& dist)
  {
    BOOST_MATH_STD_USING  // for ADL of std functions
    // TODO replace by two_div_pi<RealType>() when in boost-trunk
    //    static const RealType two_over_pi = two_div_pi<RealType>(); // = 2 / pi
    static const RealType two_over_pi = 6.366197723675814e-001; // = 2 / pi
    static const RealType sqrt_two_over_pi = 7.978845608028654e-001; // = ( 2 / pi ) ^ 0.5
    static const RealType factor = 4.292036732051034e-001; // = (4-pi)/2
    const RealType delta = dist.shape() / sqrt(1+dist.shape()*dist.shape());

    return factor * pow(sqrt_two_over_pi * delta, 3) /
      pow(1 - two_over_pi * delta * delta, static_cast<RealType>(1.5));
  }

  template <class RealType, class Policy>
  inline RealType kurtosis(const skew_normal_distribution<RealType, Policy>& dist)
  {
    return kurtosis_excess(dist)+3;
  }

  template <class RealType, class Policy>
  inline RealType kurtosis_excess(const skew_normal_distribution<RealType, Policy>& dist)
  {
    BOOST_MATH_STD_USING  // for ADL of std functions
     // TODO replace by two_div_pi<RealType>() etc when in boost-trunk
    using namespace boost::math::constants;

    static const RealType two_over_pi = 6.366197723675814e-001; // = 2 / pi
    static const RealType sqrt_two_over_pi = 7.978845608028654e-001; // = ( 2 / pi ) ^ 0.5 
    // static const RealType factor = 2.831853071795862e-001; // = 2*(pi-3)  ???  not just 2 * pi_minus_three = 0.141593
    static const RealType factor = pi_minus_three<RealType>()*static_cast<RealType>(2);

    const RealType delta = dist.shape() / sqrt(1+dist.shape()*dist.shape());

    const RealType x = 1-two_over_pi *delta*delta;

    //return 2 * pi_minus_three<RealType>() * pow(sqrt_two_over_pi * delta, 4) / (x*x) ;
     return factor * pow(sqrt_two_over_pi * delta, 4) / (x*x) ;
  }

  namespace detail
  {

    template <class RealType, class Policy>
    struct skew_normal_quantile_functor
    { 
      skew_normal_quantile_functor(const boost::math::skew_normal_distribution<RealType, Policy> dist, RealType const& p)
        : distribution(dist), prob(p)
      {
      }

      boost::math::tuple<RealType, RealType> operator()(RealType const& x)
      {
        RealType c = cdf(distribution, x);
        RealType fx = c - prob;  // Difference cdf - value - to minimize.
        RealType dx = pdf(distribution, x); // pdf is 1st derivative.
        // return both function evaluation difference f(x) and 1st derivative f'(x).
        return boost::math::make_tuple(fx, dx);
      }
    private:
      const boost::math::skew_normal_distribution<RealType, Policy> distribution;
      RealType prob; 
    };

  } // namespace detail

  template <class RealType, class Policy>
  inline RealType quantile(const skew_normal_distribution<RealType, Policy>& dist, const RealType& p)
  {
    const RealType scale = dist.scale();
    const RealType location = dist.location();
    const RealType shape = dist.shape();

    static const char* function = "boost::math::quantile(const skew_normal_distribution<%1%>&, %1%)";

    RealType result = 0;
    if(false == detail::check_scale(function, scale, &result, Policy()))
      return result;
    if(false == detail::check_location(function, location, &result, Policy()))
      return result;
    if(false == detail::check_skew_normal_shape(function, shape, &result, Policy()))
      return result;
    if(false == detail::check_probability(function, p, &result, Policy()))
      return result;

    // compute initial guess via Cornish-Fisher expansion
    RealType x = -boost::math::erfc_inv(2 * p, Policy()) * constants::root_two<RealType>();

    // avoid unnecessary computations if there is no skew
    if(shape != 0)
    {
      const RealType skew = skewness(dist);
      const RealType exk = kurtosis_excess(dist);

      x = x + (x*x-1)*skew/6 + x*(x*x-3)*exk/24 - x*(2*x*x-5)*skew*skew/36;
    } // if(shape != 0)

    result = standard_deviation(dist)*x+mean(dist);

    // handle special case of non-skew normal distribution
    if(shape == 0)
      return result;

    // refine the result by numerically searching the root of (p-cdf)

    const RealType search_min = range(dist).first;
    const RealType search_max = range(dist).second;

    const int get_digits = policies::digits<RealType, Policy>();// get digits from policy, 
    boost::uintmax_t m = policies::get_max_root_iterations<Policy>(); // and max iterations.

    result = tools::newton_raphson_iterate(detail::skew_normal_quantile_functor<RealType, Policy>(dist, p), result,
      search_min, search_max, get_digits, m);

    return result;
  } // quantile

  template <class RealType, class Policy>
  inline RealType quantile(const complemented2_type<skew_normal_distribution<RealType, Policy>, RealType>& c)
  {	
    const RealType scale = c.dist.scale();
    const RealType location = c.dist.location();
    const RealType shape = c.dist.shape();

    static const char* function = "boost::math::quantile(const complement(skew_normal_distribution<%1%>&), %1%)";
    RealType result = 0;
    if(false == detail::check_scale(function, scale, &result, Policy()))
      return result;
    if(false == detail::check_location(function, location, &result, Policy()))
      return result;
    if(false == detail::check_skew_normal_shape(function, shape, &result, Policy()))
      return result;
    RealType q = c.param;
    if(false == detail::check_probability(function, q, &result, Policy()))
      return result;

    skew_normal_distribution<RealType, Policy> D(-location, scale, -shape);

    result = -quantile(D, q);

    return result;
  } // quantile


} // namespace math
} // namespace boost

// This include must be at the end, *after* the accessors
// for this distribution have been defined, in order to
// keep compilers that support two-phase lookup happy.
#include <boost/math/distributions/detail/derived_accessors.hpp>

#endif // BOOST_STATS_SKEW_NORMAL_HPP


