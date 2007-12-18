//  (C) Copyright John Maddock 2005-2006.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_LOG1P_INCLUDED
#define BOOST_MATH_LOG1P_INCLUDED

#ifdef _MSC_VER
#pragma once
#endif

#include <cmath>
#include <math.h> // platform's ::log1p
#include <boost/limits.hpp>
#include <boost/math/tools/config.hpp>
#include <boost/math/tools/series.hpp>
#include <boost/math/policies/error_handling.hpp>
#include <boost/math/special_functions/math_fwd.hpp>

#ifndef BOOST_NO_LIMITS_COMPILE_TIME_CONSTANTS
#  include <boost/static_assert.hpp>
#else
#  include <boost/assert.hpp>
#endif

namespace boost{ namespace math{

namespace detail
{
  // Functor log1p_series returns the next term in the Taylor series
  //   pow(-1, k-1)*pow(x, k) / k
  // each time that operator() is invoked.
  //
  template <class T>
  struct log1p_series
  {
     typedef T result_type;

     log1p_series(T x)
        : k(0), m_mult(-x), m_prod(-1){}

     T operator()()
     {
        m_prod *= m_mult;
        return m_prod / ++k;
     }

     int count()const
     {
        return k;
     }

  private:
     int k;
     const T m_mult;
     T m_prod;
     log1p_series(const log1p_series&);
     log1p_series& operator=(const log1p_series&);
  };

} // namespace detail

// Algorithm log1p is part of C99, but is not yet provided by many compilers.
//
// This version uses a Taylor series expansion for 0.5 > x > epsilon, which may
// require up to std::numeric_limits<T>::digits+1 terms to be calculated. 
// It would be much more efficient to use the equivalence:
//   log(1+x) == (log(1+x) * x) / ((1-x) - 1)
// Unfortunately many optimizing compilers make such a mess of this, that 
// it performs no better than log(1+x): which is to say not very well at all.
//
template <class T, class Policy>
typename tools::promote_args<T>::type log1p(T x, const Policy& pol)
{ // The function returns the natural logarithm of 1 + x.
  // A domain error occurs if x < -1.  TODO should there be a check?
   typedef typename tools::promote_args<T>::type result_type;
   BOOST_MATH_STD_USING
   using std::abs;

   static const char* function = "boost::math::log1p<%1%>(%1%)";

   if(x < -1)
      return policies::raise_domain_error<T>(
         function, "log1p(x) requires x > -1, but got x = %1%.", x, pol);
   if(x == -1)
      return -policies::raise_overflow_error<T>(
         function, 0, pol);

   result_type a = abs(result_type(x));
   if(a > result_type(0.5L))
      return log(1 + result_type(x));
   // Note that without numeric_limits specialisation support, 
   // epsilon just returns zero, and our "optimisation" will always fail:
   if(a < tools::epsilon<result_type>())
      return x;
   detail::log1p_series<result_type> s(x);
   boost::uintmax_t max_iter = policies::get_max_series_iterations<Policy>();
#if !BOOST_WORKAROUND(__BORLANDC__, BOOST_TESTED_AT(0x582))
   result_type result = tools::sum_series(s, policies::digits<result_type, Policy>(), max_iter);
#else
   result_type zero = 0;
   result_type result = tools::sum_series(s, policies::digits<result_type, Policy>(), max_iter, zero);
#endif
   policies::check_series_iterations(function, max_iter, pol);
   return result;
}

#if BOOST_WORKAROUND(__BORLANDC__, BOOST_TESTED_AT(0x564))
// These overloads work around a type deduction bug:
inline float log1p(float z)
{
   return log1p<float>(z);
}
inline double log1p(double z)
{
   return log1p<double>(z);
}
#ifndef BOOST_MATH_NO_LONG_DOUBLE_MATH_FUNCTIONS
inline long double log1p(long double z)
{
   return log1p<long double>(z);
}
#endif
#endif

#ifdef log1p
#  ifndef BOOST_HAS_LOG1P
#     define BOOST_HAS_LOG1P
#  endif
#  undef log1p
#endif

#ifdef BOOST_HAS_LOG1P
#  if (defined(__STDC_VERSION__) && (__STDC_VERSION__ >= 199901)) \
   || ((defined(linux) || defined(__linux) || defined(__linux__)) && !defined(__SUNPRO_CC)) \
   || (defined(__hpux) && !defined(__hppa))
template <class Policy>
inline float log1p(float x, const Policy& pol)
{ 
   if(x < -1)
      return policies::raise_domain_error<float>(
         "log1p<%1%>(%1%)", "log1p(x) requires x > -1, but got x = %1%.", x, pol);
   if(x == -1)
      return -policies::raise_overflow_error<float>(
         "log1p<%1%>(%1%)", 0, pol);
   return ::log1pf(x); 
}
template <class Policy>
inline long double log1p(long double x, const Policy& pol)
{ 
   if(x < -1)
      return policies::raise_domain_error<long double>(
         "log1p<%1%>(%1%)", "log1p(x) requires x > -1, but got x = %1%.", x, pol);
   if(x == -1)
      return -policies::raise_overflow_error<long double>(
         "log1p<%1%>(%1%)", 0, pol);
   return ::log1pl(x); 
}
#else
template <class Policy>
inline float log1p(float x, const Policy& pol)
{ 
   if(x < -1)
      return policies::raise_domain_error<float>(
         "log1p<%1%>(%1%)", "log1p(x) requires x > -1, but got x = %1%.", x, pol);
   if(x == -1)
      return -policies::raise_overflow_error<float>(
         "log1p<%1%>(%1%)", 0, pol);
   return ::log1p(x); 
}
#endif
template <class Policy>
inline double log1p(double x, const Policy& pol)
{ 
   if(x < -1)
      return policies::raise_domain_error<double>(
         "log1p<%1%>(%1%)", "log1p(x) requires x > -1, but got x = %1%.", x, pol);
   if(x == -1)
      return -policies::raise_overflow_error<double>(
         "log1p<%1%>(%1%)", 0, pol);
   return ::log1p(x); 
}
#elif defined(_MSC_VER) && (BOOST_MSVC >= 1400)
//
// You should only enable this branch if you are absolutely sure
// that your compilers optimizer won't mess this code up!!
// Currently tested with VC8 and Intel 9.1.
//
template <class Policy>
inline double log1p(double x, const Policy& pol)
{
   if(x < -1)
      return policies::raise_domain_error<double>(
         "log1p<%1%>(%1%)", "log1p(x) requires x > -1, but got x = %1%.", x, pol);
   if(x == -1)
      return -policies::raise_overflow_error<double>(
         "log1p<%1%>(%1%)", 0, pol);
   double u = 1+x;
   if(u == 1.0) 
      return x; 
   else
      return log(u)*(x/(u-1.0));
}
template <class Policy>
inline float log1p(float x, const Policy& pol)
{
   return static_cast<float>(boost::math::log1p(static_cast<double>(x), pol));
}
template <class Policy>
inline long double log1p(long double x, const Policy& pol)
{
   if(x < -1)
      return policies::raise_domain_error<long double>(
         "log1p<%1%>(%1%)", "log1p(x) requires x > -1, but got x = %1%.", x, pol);
   if(x == -1)
      return -policies::raise_overflow_error<long double>(
         "log1p<%1%>(%1%)", 0, pol);
   long double u = 1+x;
   if(u == 1.0) 
      return x; 
   else
      return log(u)*(x/(u-1.0));
}
#endif

template <class T>
inline typename tools::promote_args<T>::type log1p(T x)
{
   return boost::math::log1p(x, policies::policy<>());
}
//
// Compute log(1+x)-x:
//
template <class T, class Policy>
inline typename tools::promote_args<T>::type 
   log1pmx(T x, const Policy& pol)
{
   typedef typename tools::promote_args<T>::type result_type;
   BOOST_MATH_STD_USING
   static const char* function = "boost::math::log1pmx<%1%>(%1%)";

   if(x < -1)
      return policies::raise_domain_error<T>(
         function, "log1pmx(x) requires x > -1, but got x = %1%.", x, pol);
   if(x == -1)
      return -policies::raise_overflow_error<T>(
         function, 0, pol);

   result_type a = abs(result_type(x));
   if(a > result_type(0.95L))
      return log(1 + result_type(x)) - result_type(x);
   // Note that without numeric_limits specialisation support, 
   // epsilon just returns zero, and our "optimisation" will always fail:
   if(a < tools::epsilon<result_type>())
      return -x * x / 2;
   boost::math::detail::log1p_series<T> s(x);
   s();
   boost::uintmax_t max_iter = policies::get_max_series_iterations<Policy>();
#if BOOST_WORKAROUND(__BORLANDC__, BOOST_TESTED_AT(0x582))
   T zero = 0;
   T result = boost::math::tools::sum_series(s, policies::digits<T, Policy>(), max_iter, zero);
#else
   T result = boost::math::tools::sum_series(s, policies::digits<T, Policy>(), max_iter);
#endif
   policies::check_series_iterations(function, max_iter, pol);
   return result;
}

template <class T>
inline T log1pmx(T x)
{
   return log1pmx(x, policies::policy<>());
}

} // namespace math
} // namespace boost

#endif // BOOST_MATH_LOG1P_INCLUDED



