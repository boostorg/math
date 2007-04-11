//  Copyright John Maddock 2007.
//  Copyright Paul A. Bristow 2007
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_STATS_PARETO_HPP
#define BOOST_STATS_PARETO_HPP

// http://en.wikipedia.org/wiki/Pareto_distribution
// http://www.itl.nist.gov/div898/handbook/eda/section3/eda3661.htm
// Also:
// Weisstein, Eric W. "Pareto Distribution."
// From MathWorld--A Wolfram Web Resource.
// http://mathworld.wolfram.com/ParetoDistribution.html

#include <boost/math/distributions/complement.hpp>
#include <boost/math/distributions/detail/common_error_handling.hpp>
#include <boost/math/special_functions/powm1.hpp>


#include <utility> // for BOOST_CURRENT_VALUE?

namespace boost
{
  namespace math
  {
    namespace detail
    { // Parameter checking.
      template <class RealType>
      bool check_pareto_location(
        const char* function,
        RealType location,
        RealType* result)
      {
        if((boost::math::isfinite)(location))
        { // any > 0 finite value is OK.
          if (location > 0)
          {
            return true;
          }
          else
          {
            *result = tools::domain_error<RealType>(
              function,
              "Location parameter is %1%, but must be > 0!", location);
            return false;
          }
        }
        else
        { // Not finite.
          *result = tools::domain_error<RealType>(
            function,
            "Location parameter is %1%, but must be finite!", location);
          return false;
        }
      } // bool check_pareto_location

      template <class RealType>
      bool check_pareto_shape(
        const char* function,
        RealType shape,
        RealType* result)
      {
        if((boost::math::isfinite)(shape))
        { // Any finite value is OK.
          if (shape > 0)
          {
            return true;
          }
          else
          {
            *result = tools::domain_error<RealType>(
              function,
              "Shape parameter is %1%, but must be > 0!", shape);
            return false;
          }
        }
        else
        { // Not finite.
          *result = tools::domain_error<RealType>(
            function,
            "Shape parameter is %1%, but must be finite!", shape);
          return false;
        }
      } // bool check_pareto_shape(

      template <class RealType>
      bool check_pareto_x(
        const char* function,
        RealType const& x,
        RealType* result)
      {
        if((boost::math::isfinite)(x))
        { // 
          if (x > 0)
          {
            return true;
          }
          else
          {
            *result = tools::domain_error<RealType>(
              function,
              "x parameter is %1%, but must be > 0 !", x);
            return false;
          }
        }
        else
        { // Not finite..
          *result = tools::domain_error<RealType>(
            function,
            "x parameter is %1%, but must be finite!", x);
          return false;
        }
      } // bool check_pareto_x

      template <class RealType>
      inline bool check_pareto( // distribution parameters.
        const char* function,
        RealType location,
        RealType shape,
        RealType* result)
      {
        if(check_pareto_location(function, location, result)
        && check_pareto_shape(function, shape, result) )
        {
          return true;
        }
        else
        { 
          return false;
        }
      } // bool check_pareto(

    } // namespace detail

    template <class RealType = double>
    class pareto_distribution
    {
    public:
      typedef RealType value_type;

      pareto_distribution(RealType location = 1, RealType shape = 1)
        : m_location(location), m_shape(shape)
      { // Constructor.
        RealType result;
        detail::check_pareto(BOOST_CURRENT_FUNCTION, location, shape, &result);
      }

      RealType location()const
      { // AKA Xm and b
        return m_location;
      }

      RealType shape()const
      { // AKA k and a
        return m_shape;
      }
    private:
      // Data members:
      RealType m_location;  // distribution location (xm)
      RealType m_shape;  // distribution shape (k)
    };

    typedef pareto_distribution<double> pareto; // Convenience to allow pareto(2., 3.);

    template <class RealType>
    const std::pair<RealType, RealType> range(const pareto_distribution<RealType>& /*dist*/)
    { // Range of permissible values for random variable x.
      using boost::math::tools::max_value;
      return std::pair<RealType, RealType>(0, max_value<RealType>()); // location zero to + infinity.
    } // range

    template <class RealType>
    const std::pair<RealType, RealType> support(const pareto_distribution<RealType>& dist)
    { // Range of supported values for random variable x.
      // This is range where cdf rises from 0 to 1, and outside it, the pdf is zero.
      using boost::math::tools::max_value;
      return std::pair<RealType, RealType>(dist.location(), max_value<RealType>() ); // location to + infinity.
    } // support

    template <class RealType>
    RealType pdf(const pareto_distribution<RealType>& dist, const RealType& x)
    {
      using namespace std;  // for ADL of std function pow.
      RealType location = dist.location();
      RealType result;
      detail::check_pareto_x(BOOST_CURRENT_FUNCTION, x, &result);
      if (x < location)
      { // regardless of shape, pdf is zero.
        return 0; 
      }
      RealType shape = dist.shape();

      result = shape * pow(location, shape) / pow(x, shape+1);
      return result;
    } // pdf

    template <class RealType>
    RealType cdf(const pareto_distribution<RealType>& dist, const RealType& x)
    {
      using namespace std;  // for ADL of std function pow.
      RealType location = dist.location();
      RealType result;
      detail::check_pareto_x(BOOST_CURRENT_FUNCTION, x, &result);
      if (x <= location)
      { // regardless of shape, cdf is zero.
        return 0; 
      }

      RealType shape = dist.shape();
      // result = RealType(1) - pow((location / x), shape);
      result = -powm1(location/x, shape); // should be more accurate.
      return result;
    } // cdf

    template <class RealType>
    RealType quantile(const pareto_distribution<RealType>& dist, const RealType& p)
    {
      using namespace std;  // for ADL of std function pow.
      RealType result;
      if(false == detail::check_probability(BOOST_CURRENT_FUNCTION, p, &result))
      {
        return result;
      }
      RealType location = dist.location();
      if (p == 0)
      {
        return location; // x must be location (or less).
      }
      if (p == 1)
      {
        return tools::max_value<RealType>(); // x = + infinity.
      }
      RealType shape = dist.shape();
      result = location /
        (pow((1 - p), 1 / shape));
      // K. Krishnamoorthy,  ISBN 1-58488-635-8 eq 23.1.3
      return result;
    } // quantile

    template <class RealType>
    RealType cdf(const complemented2_type<pareto_distribution<RealType>, RealType>& c)
    {
       using namespace std;  // for ADL of std function pow.
       RealType result;
       RealType x = c.param;
       if(false == detail::check_pareto_x(BOOST_CURRENT_FUNCTION, x, &result))
       {
         return result;
       }
       RealType location = c.dist.location();
       if (x <= location)
       { // regardless of shape, cdf is zero, and complement is unity.
         return 1; 
       }
       RealType shape = c.dist.shape();
       result = pow((location/x), shape);
   
       return result;
    } // cdf complement
    
    template <class RealType>
    RealType quantile(const complemented2_type<pareto_distribution<RealType>, RealType>& c)
    {
      using namespace std;  // for ADL of std function pow.
      RealType result;
      RealType q = c.param;
      if(false == detail::check_probability(BOOST_CURRENT_FUNCTION, q, &result))
      {
        return result;
      }
      RealType location = c.dist.location();
      if (q == 1)
      {
        return location; // x must be location (or less).
      }
      if (q == 0)
      {
        return tools::max_value<RealType>(); // x = + infinity.
      }
      RealType shape = c.dist.shape();
      result = location / (pow(q, 1 / shape));
      // K. Krishnamoorthy,  ISBN 1-58488-635-8 eq 23.1.3
      return result;
    } // quantile complement

    template <class RealType>
    inline RealType mean(const pareto_distribution<RealType>& dist)
    {
      if (dist.shape() > RealType(1))
      {
        return dist.shape() * dist.location() / (dist.shape() - 1);
      }
      else
      {
        using boost::math::tools::max_value;
        return max_value<RealType>(); // +infinity.
      }
    } // mean

    template <class RealType>
    inline RealType mode(const pareto_distribution<RealType>& dist)
    {
      return dist.location();
    } // mode

    template <class RealType>
    inline RealType median(const pareto_distribution<RealType>& dist)
    {
      using namespace std;
      return dist.location() * pow(2, (1/dist.shape()));
    } // median

    template <class RealType>
    inline RealType variance(const pareto_distribution<RealType>& dist)
    {
      RealType result;
      RealType location = dist.location();
      RealType shape = dist.shape();
      if (shape > 2)
      {
        result = (location * location * shape) /
         ((shape - 1) *  (shape - 1) * (shape - 2));
      }
      else
      {
        result = tools::domain_error<RealType>(
          BOOST_CURRENT_FUNCTION,
          "variance is undefined for shape <= 2, but got %1%.", dist.shape());
      }
      return result;
    } // variance

    template <class RealType>
    inline RealType skewness(const pareto_distribution<RealType>& dist)
    {  
      using namespace std;
      RealType result;
      RealType shape = dist.shape();
      if (shape > 3)
      {
        result = sqrt((shape - 2) / shape) *
          2 * (shape + 1) /
          (shape - 3);
      }
      else
      {
        result = tools::domain_error<RealType>(
          BOOST_CURRENT_FUNCTION,
          "skewness is undefined for shape <= 3, but got %1%.", dist.shape());
      }
      return result;
    } // skewness

    template <class RealType>
    inline RealType kurtosis(const pareto_distribution<RealType>& dist)
    {
      RealType result;
      RealType shape = dist.shape();
      if (shape > 4)
      {
        result = 3 * ((shape - 2) * (3 * shape * shape + shape + 2)) /
          (shape * (shape - 3) * (shape - 4));
      }
      else
      {
        result = tools::domain_error<RealType>(
          BOOST_CURRENT_FUNCTION,
          "kurtosis_excess is undefined for shape <= 4, but got %1%.", shape);
      }
      return result;
    } // kurtosis

    template <class RealType>
    inline RealType kurtosis_excess(const pareto_distribution<RealType>& dist)
    {
      RealType result;
      RealType shape = dist.shape();
      if (shape > 4)
      {
        result = 6 * ((shape * shape * shape) + (shape * shape) - 6 * shape - 2) /
          (shape * (shape - 3) * (shape - 4));
      }
      else
      {
        result = tools::domain_error<RealType>(
          BOOST_CURRENT_FUNCTION,
          "kurtosis_excess is undefined for shape <= 4, but got %1%.", dist.shape());
      }
      return result;
    } // kurtosis_excess

    } // namespace math
  } // namespace boost

  // This include must be at the end, *after* the accessors
  // for this distribution have been defined, in order to
  // keep compilers that support two-phase lookup happy.
#include <boost/math/distributions/detail/derived_accessors.hpp>

#endif // BOOST_STATS_PARETO_HPP

