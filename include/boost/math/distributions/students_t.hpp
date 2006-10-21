//  Copyright John Maddock 2006.
//  Copyright Paul A. Bristow 2006.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_STATS_STUDENTS_T_HPP
#define BOOST_STATS_STUDENTS_T_HPP

// http://en.wikipedia.org/wiki/Student%27s_t_distribution
// http://www.itl.nist.gov/div898/handbook/eda/section3/eda3664.htm

#include <boost/math/special_functions/beta.hpp> // for ibeta(a, b, x).
#include <boost/math/distributions/complement.hpp>
#include <boost/math/distributions/detail/common_error_handling.hpp>

#ifdef BOOST_MSVC
# pragma warning(push)
# pragma warning(disable: 4702) // unreachable code (return after domain_error throw).
#endif

namespace boost{ namespace math{

template <class RealType = double>
class students_t_distribution
{
public:
   typedef RealType value_type;

   students_t_distribution(RealType i) : m_df(i)
   {
      RealType result;
      detail::check_df(
         BOOST_CURRENT_FUNCTION, m_df, &result);
   } // students_t_distribution

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

typedef students_t_distribution<double> students_t;

template <class RealType>
RealType pdf(const students_t_distribution<RealType>& dist, const RealType& t)
{
   BOOST_FPU_EXCEPTION_GUARD
   using namespace std;  // for ADL of std functions

   RealType degrees_of_freedom = dist.degrees_of_freedom();
   // Error check:
   RealType error_result;
   if(false == detail::check_df(
         BOOST_CURRENT_FUNCTION, degrees_of_freedom, &error_result))
      return error_result;
	 // Might conceivably permit df = +infinity and use normal distribution.
   // TODO fails for t == 0 and df >=1e16 for ALL fp types.
   // - probably need to use normal distribution - when available.
   RealType basem1 = t * t / degrees_of_freedom;
   RealType result;
   if(basem1 < 0.125)
   {
      result = exp(-boost::math::log1p(basem1) * (1+degrees_of_freedom) / 2);
   }
   else
   {
      result = pow(1 / (1 + basem1), (degrees_of_freedom + 1) / 2);
   }
   result /= sqrt(degrees_of_freedom) * boost::math::beta(degrees_of_freedom / 2, RealType(0.5f));

   return result;
} // pdf

template <class RealType>
RealType cdf(const students_t_distribution<RealType>& dist, const RealType& t)
{
   RealType degrees_of_freedom = dist.degrees_of_freedom();
   // Error check:
   RealType error_result;
   if(false == detail::check_df(
         BOOST_CURRENT_FUNCTION, degrees_of_freedom, &error_result))
      return error_result;

   if (t == 0)
   {
     return 0.5;
   }
   //
   // Calculate probability of Student's t using the incomplete beta function.
   // probability = ibeta(degrees_of_freedom / 2, 1/2, degrees_of_freedom / (degrees_of_freedom + t*t))
   //
   // However when t is small compared to the degrees of freedom, that formula
   // suffers from rounding error, use the identity formula to work around
   // the problem:
   //
   // I[x](a,b) = 1 - I[1-x](b,a)
   //
   // and:
   //
   //     x = df / (df + t^2)
   //
   // so:
   //
   // 1 - x = t^2 / (df + t^2)
   //
   RealType t2 = t * t;
   RealType probability;
   if(degrees_of_freedom > 2 * t2)
   {
      RealType z = t2 / (degrees_of_freedom + t * t);
      probability = ibetac(static_cast<RealType>(0.5), degrees_of_freedom / 2, z) / 2;
   }
   else
   {
      RealType z = degrees_of_freedom / (degrees_of_freedom + t * t);
      probability = ibeta(degrees_of_freedom / 2, static_cast<RealType>(0.5), z) / 2;
   }
   // Check 0 <= probability probability <= 1.
   // Numerical errors might cause probability to be slightly outside the range < 0 or > 1.
   // This might cause trouble downstream, so warn, possibly throw exception, but constrain to the limits.
   if (probability < static_cast<RealType>(0.))
   {
      tools::logic_error<RealType>(BOOST_CURRENT_FUNCTION, "probability %1% is < 0, so has been constrained to zero !", probability);
      return static_cast<RealType>(0.); // Constrain to zero if logic_error does not throw.
   }
   if(probability > static_cast<RealType>(1.))
   {
      tools::logic_error<RealType>(BOOST_CURRENT_FUNCTION, "probability %1% is > 1, so has been constrained to unity!", probability);
      return static_cast<RealType>(1.); // Constrain to unity if logic_error does not throw.
   }
   return (t > 0 ? 1	- probability : probability);
} // cdf

template <class RealType>
RealType quantile(const students_t_distribution<RealType>& dist, const RealType& p)
{
   using namespace std; // for ADL of std functions
   //
   // Obtain parameters:
   //
   RealType degrees_of_freedom = dist.degrees_of_freedom();
   RealType probability = p;
   //
   // Check for domain errors:
   //
   RealType error_result;
   if(false == detail::check_df(
         BOOST_CURRENT_FUNCTION, degrees_of_freedom, &error_result)
         && detail::check_probability(BOOST_CURRENT_FUNCTION, probability, &error_result))
      return error_result;

   // Special cases, regardless of degrees_of_freedom.
   if (probability == 0)
      return numeric_limits<RealType>::has_infinity ? -numeric_limits<RealType>::infinity() : -tools::max_value<RealType>();
   if (probability == 1)
     return numeric_limits<RealType>::has_infinity ? numeric_limits<RealType>::infinity() : tools::max_value<RealType>();
   if (probability == static_cast<RealType>(0.5))
     return 0;
   //
   // Calculate quantile of Student's t using the incomplete beta function inverse:
   //
   probability = (probability > 0.5) ? 1 - probability : probability;
   RealType t, x, y;
   x = ibeta_inv(degrees_of_freedom / 2, RealType(0.5), 2 * probability, &y);
   if(degrees_of_freedom * y > tools::max_value<RealType>() * x)
      t = numeric_limits<RealType>::has_infinity ? numeric_limits<RealType>::infinity() : tools::max_value<RealType>();
   else
      t = sqrt(degrees_of_freedom * y / x);
   //
   // Figure out sign based on the size of p:
   //
   if(p < 0.5)
      t = -t;

   return t;
} // quantile

template <class RealType>
RealType cdf(const complemented2_type<students_t_distribution<RealType>, RealType>& c)
{
   return cdf(c.dist, -c.param);
}

template <class RealType>
RealType quantile(const complemented2_type<students_t_distribution<RealType>, RealType>& c)
{
   return -quantile(c.dist, c.param);
}

//
// Parameter estimation follows:
//
namespace detail{
//
// Functors for finding degrees of freedom:
//
template <class RealType>
struct sample_size_func
{
   sample_size_func(RealType a, RealType b, RealType s, RealType d)
      : alpha(a), beta(b), ratio(s*s/(d*d)) {}

   RealType operator()(const RealType& df)
   {
      if(df <= tools::min_value<RealType>())
         return 1;
      students_t_distribution<RealType> t(df);
      RealType qa = quantile(complement(t, alpha));
      RealType qb = quantile(complement(t, beta));
      qa += qb;
      qa *= qa;
      qa *= ratio;
      qa -= (df + 1);
      return qa;
   }
   RealType alpha, beta, ratio;
};

}  // namespace detail

template <class RealType>
RealType students_t_distribution<RealType>::estimate_degrees_of_freedom(
      RealType difference_from_mean,
      RealType alpha,
      RealType beta,
      RealType sd,
      RealType hint)
{
   //
   // Check for domain errors:
   //
   RealType error_result;
   if(false == detail::check_probability(
         BOOST_CURRENT_FUNCTION, alpha, &error_result)
         && detail::check_probability(BOOST_CURRENT_FUNCTION, beta, &error_result))
      return error_result;

   if(hint <= 0)
      hint = 1;

   detail::sample_size_func<RealType> f(alpha, beta, sd, difference_from_mean);
   tools::eps_tolerance<RealType> tol(tools::digits<RealType>());
   boost::uintmax_t max_iter = 10000;
   std::pair<RealType, RealType> r = tools::bracket_and_solve_root(f, hint, RealType(2), false, tol, max_iter);
   RealType result = r.first + (r.second - r.first) / 2;
   if(max_iter == 10000)
   {
      tools::logic_error<RealType>(BOOST_CURRENT_FUNCTION, "Unable to locate solution in a reasonable time:"
         " either there is no answer to how many degrees of freedom are required"
         " or the answer is infinite.  Current best guess is %1%", result);
   }
   return result;
}

template <class RealType>
inline RealType mean(const students_t_distribution<RealType>& )
{
   return 0;
}

template <class RealType>
inline RealType variance(const students_t_distribution<RealType>& dist)
{
   // Error check:
   RealType error_result;
   if(false == detail::check_df(
         BOOST_CURRENT_FUNCTION, dist.degrees_of_freedom(), &error_result))
      return error_result;

   RealType v = dist.degrees_of_freedom();
   return v / (v - 2);
}

template <class RealType>
inline RealType mode(const students_t_distribution<RealType>& dist)
{
   return 0;
}

template <class RealType>
inline RealType skewness(const students_t_distribution<RealType>& dist)
{
   if(dist.degrees_of_freedom() <= 3)
   {
      tools::domain_error<RealType>(
         BOOST_CURRENT_FUNCTION,
         "Skewness is undefined for degrees of freedom <= 3, but got %1%.",
         dist.degrees_of_freedom());
   }
   return 0;
}

template <class RealType>
inline RealType kurtosis(const students_t_distribution<RealType>& dist)
{
   RealType df = dist.degrees_of_freedom();
   if(df <= 3)
   {
      tools::domain_error<RealType>(
         BOOST_CURRENT_FUNCTION,
         "Skewness is undefined for degrees of freedom <= 3, but got %1%.",
         df);
   }
   return 3 * (df - 2) / (df - 4);
}

template <class RealType>
inline RealType kurtosis_excess(const students_t_distribution<RealType>& dist)
{
   // see http://mathworld.wolfram.com/Kurtosis.html
   RealType df = dist.degrees_of_freedom();
   if(df <= 3)
   {
      tools::domain_error<RealType>(
         BOOST_CURRENT_FUNCTION,
         "Skewness is undefined for degrees of freedom <= 3, but got %1%.",
         df);
   }
   return 6 / (df - 4);
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

#endif // BOOST_STATS_STUDENTS_T_HPP
