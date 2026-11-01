
//  (C) Copyright John Maddock 2006.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_SPECIAL_SPHERICAL_HARMONIC_HPP
#define BOOST_MATH_SPECIAL_SPHERICAL_HARMONIC_HPP

#ifdef _MSC_VER
#pragma once
#endif

#include <boost/math/special_functions/math_fwd.hpp>
#include <boost/math/special_functions/legendre.hpp>
#include <boost/math/tools/workaround.hpp>
#ifndef BOOST_MATH_BUILD_MODULE
#include <complex>
#endif

BOOST_MATH_NAMESPACE_BEGIN

namespace detail{

//
// Calculates the prefix term that's common to the real
// and imaginary parts.  Does *not* fix up the sign of the result
// though.
//
template <class T, class Policy>
inline T spherical_harmonic_prefix(unsigned n, unsigned m, T theta, const Policy& pol)
{
   BOOST_MATH_STD_USING

   if(m > n)
      return 0;

   T sin_theta = sin(theta);
   T x = cos(theta);

   if(n + m > max_factorial<T>::value)
   {
      // (n-m)!/(n+m)! may underflow and P_n^m overflow, so recur in l on the
      // normalised sqrt((2l+1)/(4pi) (l-m)!/(l+m)!) P_l^m, starting from l = m.
      // sin^m(theta) may underflow too, so p0 and p1 are scaled by 2^-e2:
      int e2 = 0;
      const T s = frexp(fabs(sin_theta), &e2);
      e2 *= static_cast<int>(m);
      T p0 = sqrt((2 * m + 1) * BOOST_MATH_NAMESPACE::tgamma_delta_ratio(T(m) + 0.5f, T(0.5f), pol) / (4 * constants::pi<T>() * constants::root_pi<T>()));
      if(m & 1)
         p0 = -p0;
      for(unsigned k = 0; k < m; k += 64)
      {
         // s^m in pieces of at most s^64 >= 2^-64, moving the exponent of each into e2:
         int e = 0;
         p0 = frexp(p0 * pow(s, T((std::min)(m - k, 64u))), &e);
         e2 += e;
      }
      T a_prev = sqrt(T(2 * m + 3));
      T p1 = a_prev * x * p0;
      for(unsigned l = m + 2; l <= n; ++l)
      {
         const T a = sqrt((4 * T(l) * l - 1) / ((T(l) - m) * (T(l) + m)));
         std::swap(p0, p1);
         p1 = a * (x * p0 - p1 / a_prev);
         a_prev = a;
         if((e2 < 0) && (fabs(p1) > 1))
         {
            int e = 0;
            (void)frexp(p1, &e);
            e = (std::min)(e, -e2);
            p0 = ldexp(p0, -e);
            p1 = ldexp(p1, -e);
            e2 += e;
         }
      }
      return ldexp(n == m ? p0 : p1, e2);
   }

   T leg = detail::legendre_p_imp(n, m, x, static_cast<T>(pow(fabs(sin_theta), T(m))), pol);
   
   T prefix = BOOST_MATH_NAMESPACE::tgamma_delta_ratio(static_cast<T>(n - m + 1), static_cast<T>(2 * m), pol);
   prefix *= (2 * n + 1) / (4 * constants::pi<T>());
   prefix = sqrt(prefix);
   return prefix * leg;
}
//
// Real Part:
//
template <class T, class Policy>
T spherical_harmonic_r(unsigned n, int m, T theta, T phi, const Policy& pol)
{
   BOOST_MATH_STD_USING  // ADL of std functions

   bool sign = false;
   if(m < 0)
   {
      // Reflect and adjust sign if m < 0:
      sign = m&1;
      m = abs(m);
   }
   if(m&1)
   {
      // Check phase if theta is outside [0, PI]:
      T mod = BOOST_MATH_NAMESPACE::tools::fmod_workaround(theta, T(2 * constants::pi<T>()));
      if(mod < 0)
         mod += 2 * constants::pi<T>();
      if(mod > constants::pi<T>())
         sign = !sign;
   }
   // Get the value and adjust sign as required:
   T prefix = spherical_harmonic_prefix(n, m, theta, pol);
   prefix *= cos(m * phi);
   return sign ? T(-prefix) : prefix;
}

template <class T, class Policy>
T spherical_harmonic_i(unsigned n, int m, T theta, T phi, const Policy& pol)
{
   BOOST_MATH_STD_USING  // ADL of std functions

   bool sign = false;
   if(m < 0)
   {
      // Reflect and adjust sign if m < 0:
      sign = !(m&1);
      m = abs(m);
   }
   if(m&1)
   {
      // Check phase if theta is outside [0, PI]:
      T mod = BOOST_MATH_NAMESPACE::tools::fmod_workaround(theta, T(2 * constants::pi<T>()));
      if(mod < 0)
         mod += 2 * constants::pi<T>();
      if(mod > constants::pi<T>())
         sign = !sign;
   }
   // Get the value and adjust sign as required:
   T prefix = spherical_harmonic_prefix(n, m, theta, pol);
   prefix *= sin(m * phi);
   return sign ? T(-prefix) : prefix;
}

template <class T, class U, class Policy>
std::complex<T> spherical_harmonic(unsigned n, int m, U theta, U phi, const Policy& pol)
{
   BOOST_MATH_STD_USING
   //
   // Sort out the signs:
   //
   bool r_sign = false;
   bool i_sign = false;
   if(m < 0)
   {
      // Reflect and adjust sign if m < 0:
      r_sign = m&1;
      i_sign = !(m&1);
      m = abs(m);
   }
   if(m&1)
   {
      // Check phase if theta is outside [0, PI]:
      U mod = BOOST_MATH_NAMESPACE::tools::fmod_workaround(theta, U(2 * constants::pi<U>()));
      if(mod < 0)
         mod += 2 * constants::pi<U>();
      if(mod > constants::pi<U>())
      {
         r_sign = !r_sign;
         i_sign = !i_sign;
      }
   }
   //
   // Calculate the value:
   //
   U prefix = spherical_harmonic_prefix(n, m, theta, pol);
   U r = prefix * cos(m * phi);
   U i = prefix * sin(m * phi);
   //
   // Add in the signs:
   //
   if(r_sign)
      r = -r;
   if(i_sign)
      i = -i;
   static const char* function = "boost::math::spherical_harmonic<%1%>(int, int, %1%, %1%)";
   return std::complex<T>(policies::checked_narrowing_cast<T, Policy>(r, function), policies::checked_narrowing_cast<T, Policy>(i, function));
}

} // namespace detail

BOOST_MATH_EXPORT template <class T1, class T2, class Policy>
inline std::complex<typename tools::promote_args<T1, T2>::type> 
   spherical_harmonic(unsigned n, int m, T1 theta, T2 phi, const Policy& pol)
{
   typedef typename tools::promote_args<T1, T2>::type result_type;
   typedef typename policies::evaluation<result_type, Policy>::type value_type;
   return detail::spherical_harmonic<result_type, value_type>(n, m, static_cast<value_type>(theta), static_cast<value_type>(phi), pol);
}

BOOST_MATH_EXPORT template <class T1, class T2>
inline std::complex<typename tools::promote_args<T1, T2>::type> 
   spherical_harmonic(unsigned n, int m, T1 theta, T2 phi)
{
   return BOOST_MATH_NAMESPACE::spherical_harmonic(n, m, theta, phi, policies::policy<>());
}

BOOST_MATH_EXPORT template <class T1, class T2, class Policy>
inline typename tools::promote_args<T1, T2>::type 
   spherical_harmonic_r(unsigned n, int m, T1 theta, T2 phi, const Policy& pol)
{
   typedef typename tools::promote_args<T1, T2>::type result_type;
   typedef typename policies::evaluation<result_type, Policy>::type value_type;
   return policies::checked_narrowing_cast<result_type, Policy>(detail::spherical_harmonic_r(n, m, static_cast<value_type>(theta), static_cast<value_type>(phi), pol), "boost::math::spherical_harmonic_r<%1%>(unsigned, int, %1%, %1%)");
}

BOOST_MATH_EXPORT template <class T1, class T2>
inline typename tools::promote_args<T1, T2>::type 
   spherical_harmonic_r(unsigned n, int m, T1 theta, T2 phi)
{
   return BOOST_MATH_NAMESPACE::spherical_harmonic_r(n, m, theta, phi, policies::policy<>());
}

BOOST_MATH_EXPORT template <class T1, class T2, class Policy>
inline typename tools::promote_args<T1, T2>::type 
   spherical_harmonic_i(unsigned n, int m, T1 theta, T2 phi, const Policy& pol)
{
   typedef typename tools::promote_args<T1, T2>::type result_type;
   typedef typename policies::evaluation<result_type, Policy>::type value_type;
   return policies::checked_narrowing_cast<result_type, Policy>(detail::spherical_harmonic_i(n, m, static_cast<value_type>(theta), static_cast<value_type>(phi), pol), "boost::math::spherical_harmonic_i<%1%>(unsigned, int, %1%, %1%)");
}

BOOST_MATH_EXPORT template <class T1, class T2>
inline typename tools::promote_args<T1, T2>::type 
   spherical_harmonic_i(unsigned n, int m, T1 theta, T2 phi)
{
   return BOOST_MATH_NAMESPACE::spherical_harmonic_i(n, m, theta, phi, policies::policy<>());
}

BOOST_MATH_NAMESPACE_END

#endif // BOOST_MATH_SPECIAL_SPHERICAL_HARMONIC_HPP



