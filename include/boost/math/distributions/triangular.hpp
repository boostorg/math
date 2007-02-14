//  Copyright John Maddock 2006.
//  Copyright Paul A. Bristow 2006.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_STATS_TRIANGULAR_HPP
#define BOOST_STATS_TRIANGULAR_HPP

// http://mathworld.wolfram.com/TriangularDistribution.html
// http://en.wikipedia.org/wiki/Triangular_distribution

#include <boost/math/special_functions/expm1.hpp>
#include <boost/math/distributions/detail/common_error_handling.hpp>
#include <boost/math/distributions/complement.hpp>
#include <boost/math/constants/constants.hpp>

#include <utility>

namespace boost{ namespace math
{
  namespace detail
  {
    template <class RealType>
    bool check_triangular_lower(
      const char* function,
      RealType lower,
      RealType* result)
    {
      if((boost::math::isfinite)(lower))
      { // Any finite value is OK.
        return true;
      }
      else
      { // Not finite: infinity or NaN.
        *result = tools::domain_error<RealType>(
          function,
          "Lower parameter is %1%, but must be finite!", lower);
        return false;
      }
    } // bool check_triangular_lower(

    template <class RealType>
    bool check_triangular_mode(
      const char* function,
      RealType mode,
      RealType* result)
    {
      if((boost::math::isfinite)(mode))
      { // any finite value is OK.
        return true;
      }
      else
      { // Not finite: infinity or NaN.
        *result = tools::domain_error<RealType>(
          function,
          "Mode parameter is %1%, but must be finite!", mode);
        return false;
      }
    } // bool check_triangular_mode(

    template <class RealType>
    bool check_triangular_upper(
      const char* function,
      RealType upper,
      RealType* result)
    {
      if((boost::math::isfinite)(upper))
      { // any finite value is OK.
        return true;
      }
      else
      { // Not finite: infinity or NaN.
        *result = tools::domain_error<RealType>(
          function,
          "Upper parameter is %1%, but must be finite!", upper);
        return false;
      }
    } // bool check_triangular_upper(

    template <class RealType>
    bool check_triangular_x(
      const char* function,
      RealType const& x,
      RealType* result)
    {
      if((boost::math::isfinite)(x))
      { // Any finite value is OK
        return true;
      }
      else
      { // Not finite: infinity or NaN.
        *result = tools::domain_error<RealType>(
          function,
          "x parameter is %1%, but must be finite!", x);
        return false;
      }
    } // bool check_triangular_x

    template <class RealType>
    inline bool check_triangular(
      const char* function,
      RealType lower,
      RealType mode,
      RealType upper,
      RealType* result)
    {
      if(check_triangular_lower(function, lower, result)
        && check_triangular_mode(function, mode, result)
        && check_triangular_upper(function, upper, result)
        && (lower < upper) // lower == upper NOT useful.
        )
      {
        if (mode < lower)
        {
          *result = tools::domain_error<RealType>(
            function,
            "mode parameter is %1%, but must be >= than lower!", lower);
          return false;
        }
        if (mode > upper )
        {
          *result = tools::domain_error<RealType>(
            function,
            "mode parameter is %1%, but must be <= than upper!", upper);
          return false;
        }
        return true;
      }
      else
      { // upper and lower have each been checked before, so must be lower >= upper.
        *result = tools::domain_error<RealType>(
          function,
          "lower parameter is %1%, but must be less than upper!", lower);
        return false;
      }
    } // bool check_triangular
  } // namespace detail

  template <class RealType = double>
  class triangular_distribution
  {
  public:
    typedef RealType value_type;

    triangular_distribution(RealType lower = -1, RealType mode = 0, RealType upper = 1)
      : m_lower(lower), m_mode(mode), m_upper(upper) // Constructor.
    { // Evans says standard is lower 0, mode 1/2, upper 1,
      // has median sqrt(c/2) for c <=1/2 and 1 - sqrt(1-c)/2 for c >= 1/2
      // But this is more useful in most applications to approximate normal distribution,
      // where the central value is the most likely and deviations either side equally likely.
      RealType result;
      detail::check_triangular(BOOST_CURRENT_FUNCTION,lower, mode, upper, &result);
    }
    // Accessor functions.
    RealType lower()const
    {
      return m_lower;
    }
    RealType mode()const
    {
      return m_mode;
    }
    RealType upper()const
    {
      return m_upper;
    }
  private:
    // Data members:
    RealType m_lower;  // distribution lower aka a
    RealType m_mode;  // distribution mode aka c
    RealType m_upper;  // distribution upper aka b
  }; // class triangular_distribution

  typedef triangular_distribution<double> triangular;

  template <class RealType>
  const std::pair<RealType, RealType> range(const triangular_distribution<RealType>& /* dist */)
  { // Range of permissible values for random variable x.
    using boost::math::tools::max_value;
    return std::pair<RealType, RealType>(-max_value<RealType>(), max_value<RealType>());
  }

  template <class RealType>
  const std::pair<RealType, RealType> support(const triangular_distribution<RealType>& dist)
  { // Range of supported values for random variable x.
    // This is range where cdf rises from 0 to 1, and outside it, the pdf is zero.
    return std::pair<RealType, RealType>(dist.lower(), dist.upper());
  }

  template <class RealType>
  RealType pdf(const triangular_distribution<RealType>& dist, const RealType& x)
  {
    RealType lower = dist.lower();
    RealType mode = dist.mode();
    RealType upper = dist.upper();
    RealType result; // of checks.
    if(false == detail::check_triangular(BOOST_CURRENT_FUNCTION, lower, mode, upper, &result))
    {
      return result;
    }
    if(false == detail::check_triangular_x(BOOST_CURRENT_FUNCTION, x, &result))
    {
      return result;
    }
    if((x < lower) || (x > upper))
    {
      return 0;
    }
    if (x == lower)
    { // (mode - lower) == 0 which would lead to divide by zero!
      return (mode == lower) ? 2 / (upper - lower) : 0;
    }
    else if (x == upper)
    {
      return (mode == upper) ? 2 / (upper - lower) : 0;
    }
    else if (x <= mode)
    {
      return 2 * (x - lower) / ((upper - lower) * (mode - lower));
    }
    else
    {  // (x > mode)
      return 2 * (upper - x) / ((upper - lower) * (upper - mode));
    }
  } // RealType pdf(const triangular_distribution<RealType>& dist, const RealType& x)

  template <class RealType>
  inline RealType cdf(const triangular_distribution<RealType>& dist, const RealType& x)
  {
    RealType lower = dist.lower();
    RealType mode = dist.mode();
    RealType upper = dist.upper();
    RealType result; // of checks.
    if(false == detail::check_triangular(BOOST_CURRENT_FUNCTION, lower, mode, upper, &result))
    {
      return result;
    }
    if(false == detail::check_triangular_x(BOOST_CURRENT_FUNCTION, x, &result))
    {
      return result;
    }
    if((x <= lower))
    {
      return 0;
    }
    if (x >= upper)
    {
      return 1;
    }
    // else lower < x < upper
    if (x <= mode)
    {
      return ((x - lower) * (x - lower)) / ((upper - lower) * (mode - lower));
    }
    else
    {
      return 1 - (upper - x) *  (upper - x) / ((upper - lower) * (upper - mode));
    }
  } // RealType cdf(const triangular_distribution<RealType>& dist, const RealType& x)

  template <class RealType>
  RealType quantile(const triangular_distribution<RealType>& dist, const RealType& p)
  {
    using namespace std;  // for ADL of std functions (sqrt).
    RealType lower = dist.lower();
    RealType mode = dist.mode();
    RealType upper = dist.upper();
    RealType result; // of checks
    if(false == detail::check_triangular(BOOST_CURRENT_FUNCTION,lower, mode, upper, &result))
    {
      return result;
    }
    if(false == detail::check_probability(BOOST_CURRENT_FUNCTION, p, &result))
    {
      return result;
    }
    if(p == 0)
    {
      return lower;
    }
    if(p == 1)
    {
      return upper;
    }
    RealType p0 = (mode - lower) / (upper - lower);
    RealType q = 1 - p;
    if (p < p0)
    {
      result = sqrt((upper - lower) * (mode - lower) * p) + lower;
    }
    else if (p == p0)
    {
      result = mode;
    }
    else // p > p0
    {
      result = upper - sqrt((upper - lower) * (upper - mode) * q);
    }
    return result;

  } // RealType quantile(const triangular_distribution<RealType>& dist, const RealType& q)

  template <class RealType>
  RealType cdf(const complemented2_type<triangular_distribution<RealType>, RealType>& c)
  {
    RealType lower = c.dist.lower();
    RealType mode = c.dist.mode();
    RealType upper = c.dist.upper();
    RealType x = c.param;
    RealType result; // of checks.
    if(false == detail::check_triangular(BOOST_CURRENT_FUNCTION, lower, mode, upper, &result))
    {
      return result;
    }
    if(false == detail::check_triangular_x(BOOST_CURRENT_FUNCTION, x, &result))
    {
      return result;
    }
    if (x <= lower)
    {
      return 1;
    }
    if (x >= upper)
    {
      return 0;
    }
    if (x <= mode)
    {
      return 1 - ((x - lower) * (x - lower)) / ((upper - lower) * (mode - lower));
    }
    else
    {
      return (upper - x) *  (upper - x) / ((upper - lower) * (upper - mode));
    }
  } // RealType cdf(const complemented2_type<triangular_distribution<RealType>, RealType>& c)

  template <class RealType>
  RealType quantile(const complemented2_type<triangular_distribution<RealType>, RealType>& c)
  {
    using namespace std;  // Aid ADL for sqrt.
    RealType l = c.dist.lower();
    RealType m = c.dist.mode();
    RealType u = c.dist.upper();
    RealType q = c.param; // probability 0 to 1.
    RealType result; // of checks.
    if(false == detail::check_triangular(BOOST_CURRENT_FUNCTION, l, m, u, &result))
    {
      return result;
    }
    if(false == detail::check_probability(BOOST_CURRENT_FUNCTION, q, &result))
    {
      return result;
    }
    if(q == 0)
    {
      return u;
    }
    if(q == 1)
    {
      return l;
    }
    RealType lower = c.dist.lower();
    RealType mode = c.dist.mode();
    RealType upper = c.dist.upper();

    RealType p = 1 - q;
    RealType p0 = (mode - lower) / (upper - lower);
    if(p < p0)
    {
      RealType s = (upper - lower) * (mode - lower);
      s *= p;
      result = sqrt((upper - lower) * (mode - lower) * p) + lower;
    }
    else if (p == p0)
    {
      result = mode;
    }
    else // p > p0
    {
      result = upper - sqrt((upper - lower) * (upper - mode) * q);
    }
    return result;
  } // RealType quantile(const complemented2_type<triangular_distribution<RealType>, RealType>& c)

  template <class RealType>
  inline RealType mean(const triangular_distribution<RealType>& dist)
  {
    RealType lower = dist.lower();
    RealType mode = dist.mode();
    RealType upper = dist.upper();
    RealType result;  // of checks.
    if(false == detail::check_triangular(BOOST_CURRENT_FUNCTION, lower, mode, upper, &result))
    {
      return result;
    }
    return (lower + upper + mode) / 3;
  } // RealType mean(const triangular_distribution<RealType>& dist)


  template <class RealType>
  RealType variance(const triangular_distribution<RealType>& dist)
  {
    RealType lower = dist.lower();
    RealType mode = dist.mode();
    RealType upper = dist.upper();
    RealType result; // of checks.
    if(false == detail::check_triangular(BOOST_CURRENT_FUNCTION, lower, mode, upper, &result))
    {
      return result;
    }
    return (lower * lower + upper * upper + mode * mode - lower * upper - lower * mode - upper * mode) / 18;
  } // RealType variance(const triangular_distribution<RealType>& dist)

  template <class RealType>
  inline RealType mode(const triangular_distribution<RealType>& dist)
  {
    RealType mode = dist.mode();
    RealType result; // of checks.
    if(false == detail::check_triangular_mode(BOOST_CURRENT_FUNCTION, mode, &result))
    { // This should never happen!
      return result;
    }
    return mode;
  } // RealType mode

  template <class RealType>
  inline RealType median(const triangular_distribution<RealType>& dist)
  {
    using namespace std; // ADL of std functions.
    RealType mode = dist.mode();
    RealType result; // of checks.
    if(false == detail::check_triangular_mode(BOOST_CURRENT_FUNCTION, mode, &result))
    { // This should never happen!
      return result;
    }
    RealType lower = dist.lower();
    RealType upper = dist.upper();
    if (mode < (upper - lower) / 2)
    {
      return lower + sqrt((upper - lower) * (mode - lower)) / constants::root_two<RealType>();
    }
    else
    {
      return upper - sqrt((upper - lower) * (upper - mode)) / constants::root_two<RealType>();
    }
  } // RealType mode

  template <class RealType>
  inline RealType skewness(const triangular_distribution<RealType>& dist)
  {
    using namespace std;  // for ADL of std functions
    using namespace boost::math::constants; // for root_two

    RealType lower = dist.lower();
    RealType mode = dist.mode();
    RealType upper = dist.upper();
    RealType result; // of checks.
    if(false == detail::check_triangular(BOOST_CURRENT_FUNCTION,lower, mode, upper, &result))
    {
      return result;
    }
    return root_two<RealType>() * (lower + upper - 2 * mode) * (2 * lower - upper - mode) * (lower - 2 * upper + mode) /
      (5 * pow((lower * lower + upper + upper + mode * mode - lower * upper - lower * mode - upper * mode), RealType(3)/RealType(2)));
  } // RealType skewness(const triangular_distribution<RealType>& dist)

  template <class RealType>
  inline RealType kurtosis(const triangular_distribution<RealType>& dist)
  { // These checks may be belt and braces as should have been checked on construction?
    RealType lower = dist.lower();
    RealType upper = dist.upper();
    RealType mode = dist.mode();
    RealType result;  // of checks.
    if(false == detail::check_triangular(BOOST_CURRENT_FUNCTION,lower, mode, upper, &result))
    {
      return result;
    }
    return static_cast<RealType>(12)/5; //  12/5 = 2.4;
  } // RealType kurtosis_excess(const triangular_distribution<RealType>& dist)

  template <class RealType>
  inline RealType kurtosis_excess(const triangular_distribution<RealType>& dist)
  { // These checks may be belt and braces as should have been checked on construction?
    RealType lower = dist.lower();
    RealType upper = dist.upper();
    RealType mode = dist.mode();
    RealType result;  // of checks.
    if(false == detail::check_triangular(BOOST_CURRENT_FUNCTION,lower, mode, upper, &result))
    {
      return result;
    }
    return static_cast<RealType>(-3)/5; // - 3/5 = -0.6
    // Assuming mathworld really means kurtosis excess?  Wikipedia corrected to match this.
  }

} // namespace math
} // namespace boost

// This include must be at the end, *after* the accessors
// for this distribution have been defined, in order to
// keep compilers that support two-phase lookup happy.
#include <boost/math/distributions/detail/derived_accessors.hpp>

#endif // BOOST_STATS_TRIANGULAR_HPP


