
//  (C) Copyright John Maddock 2006.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_SPECIAL_SPHERICAL_HARMONIC_HPP
#define BOOST_MATH_SPECIAL_SPHERICAL_HARMONIC_HPP

#include <boost/math/special_functions/legendre.hpp>
#include <complex>

namespace boost{
namespace math{

namespace detail{

//
// Calculates the prefix term that's common to the real
// and imaginary parts.  Does *not* fix up the sign of the result
// though.
//
template <class T>
T spherical_harmonic_prefix(unsigned n, unsigned m, T theta)
{
   using namespace std;

   if(m > n)
      return 0;

   T sin_theta = sin(theta);
   T x = cos(theta);

   T leg = detail::legendre_p_imp(n, m, x, pow(fabs(sin_theta), T(m)));
   
   T prefix = tgamma_delta_ratio(static_cast<T>(n - m + 1), static_cast<T>(2 * m));
   prefix *= (2 * n + 1) / (4 * constants::pi<T>());
   prefix = sqrt(prefix);
   return prefix * leg;
}
//
// Real Part:
//
template <class T>
T spherical_harmonic_r(unsigned n, int m, T theta, T phi)
{
   using namespace std;  // ADL of std functions

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
      T mod = fmod(theta, 2 * constants::pi<T>());
      if(mod < 0)
         mod += 2 * constants::pi<T>();
      if(mod > constants::pi<T>())
         sign = !sign;
   }
   // Get the value and adjust sign as required:
   T prefix = spherical_harmonic_prefix(n, m, theta);
   prefix *= cos(m * phi);
   return sign ? -prefix : prefix;
}

template <class T>
T spherical_harmonic_i(unsigned n, int m, T theta, T phi)
{
   using namespace std;  // ADL of std functions

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
      T mod = fmod(theta, 2 * constants::pi<T>());
      if(mod < 0)
         mod += 2 * constants::pi<T>();
      if(mod > constants::pi<T>())
         sign = !sign;
   }
   // Get the value and adjust sign as required:
   T prefix = spherical_harmonic_prefix(n, m, theta);
   prefix *= sin(m * phi);
   return sign ? -prefix : prefix;
}

template <class T, class U>
std::complex<T> spherical_harmonic(unsigned n, int m, U theta, U phi)
{
   using namespace std;
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
      U mod = fmod(theta, 2 * constants::pi<U>());
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
   U prefix = spherical_harmonic_prefix(n, m, theta);
   U r = prefix * cos(m * phi);
   U i = prefix * sin(m * phi);
   //
   // Add in the signs:
   //
   if(r_sign)
      r = -r;
   if(i_sign)
      i = -i;
   return std::complex<T>(tools::checked_narrowing_cast<T>(r, BOOST_CURRENT_FUNCTION), tools::checked_narrowing_cast<T>(i, BOOST_CURRENT_FUNCTION));
}

} // namespace detail

template <class T>
inline std::complex<T> spherical_harmonic(unsigned n, int m, T theta, T phi)
{
   typedef typename tools::evaluation<typename remove_cv<T>::type>::type value_type;
   return detail::spherical_harmonic<T, value_type>(n, m, static_cast<value_type>(theta), static_cast<value_type>(phi));
}

template <class T>
inline T spherical_harmonic_r(unsigned n, int m, T theta, T phi)
{
   typedef typename tools::evaluation<typename remove_cv<T>::type>::type value_type;
   return tools::checked_narrowing_cast<typename remove_cv<T>::type>(detail::spherical_harmonic_r(n, m, static_cast<value_type>(theta), static_cast<value_type>(phi)), BOOST_CURRENT_FUNCTION);
}

template <class T>
inline T spherical_harmonic_i(unsigned n, int m, T theta, T phi)
{
   typedef typename tools::evaluation<typename remove_cv<T>::type>::type value_type;
   return tools::checked_narrowing_cast<typename remove_cv<T>::type>(detail::spherical_harmonic_i(n, m, static_cast<value_type>(theta), static_cast<value_type>(phi)), BOOST_CURRENT_FUNCTION);
}

} // namespace math
} // namespace boost

#endif // BOOST_MATH_SPECIAL_SPHERICAL_HARMONIC_HPP

