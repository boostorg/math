//  Copyright (c) 2006 Xiaogang Zhang
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
//  History:
//  XZ wrote the original of this file as part of the Google
//  Summer of Code 2006.  JM modified it to fit into the
//  Boost.Math conceptual framework better, and to correctly
//  handle the various corner cases.
//

#ifndef BOOST_MATH_ELLINT_3_HPP
#define BOOST_MATH_ELLINT_3_HPP

#include <boost/math/special_functions/ellint_rf.hpp>
#include <boost/math/special_functions/ellint_rj.hpp>
#include <boost/math/special_functions/ellint_1.hpp>
#include <boost/math/special_functions/ellint_2.hpp>
#include <boost/math/special_functions/log1p.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/math/tools/error_handling.hpp>

// Elliptic integrals (complete and incomplete) of the third kind
// Carlson, Numerische Mathematik, vol 33, 1 (1979)

namespace boost { namespace math { namespace detail{

template <typename T>
T ellint_pi_imp(T v, T k, T vc);

// Elliptic integral (Legendre form) of the third kind
template <typename T>
T ellint_pi_imp(T v, T phi, T k, T vc)
{
    // Note vc = 1-v presumably without cancellation error.
    T value, x, y, z, p, t;

    using namespace std;
    using namespace boost::math::tools;
    using namespace boost::math::constants;

    if (abs(k) > 1)
    {
        return domain_error<T>(BOOST_CURRENT_FUNCTION,
            "Got k = %1%, function requires |k| <= 1", k);
    }

    T sphi = sin(fabs(phi));
    if(v > 1 / (sphi * sphi))
    {
        // Complex result is a domain error:
        return domain_error<T>(BOOST_CURRENT_FUNCTION,
            "Got v = %1%, but result is complex for v > 1 / sin^2(phi)", v);
    }

    // Special cases first:
    if(v == 0)
    {
       // A&S 17.7.18 & 19
       return (k == 0) ? phi : ellint_f_imp(phi, k);
    }
    if(phi == constants::pi<T>() / 2)
    {
       // Have to filter this case out before the next
       // special case, otherwise we might get an infinity from
       // tan(phi).
       // Also note that since we can't represent PI/2 exactly
       // in a T, this is a bit of a guess as to the users true
       // intent...
       //
       return ellint_pi_imp(v, k, vc);
    }
    if(k == 0)
    {
       // A&S 17.7.20:
       if(v < 1)
       {
          T vcr = sqrt(vc);
          return atan(vcr * tan(phi)) / vcr;
       }
       else if(v == 1)
       {
          return tan(phi);
       }
       else
       {
          // v > 1:
          T vcr = sqrt(-vc);
          T arg = vcr * tan(phi);
          return (boost::math::log1p(arg) - boost::math::log1p(-arg)) / (2 * vcr);
       }
    }

    if(v < 0)
    {
       //
       // If we don't shift to 0 <= v <= 1 we get
       // cancellation errors later on.  Use
       // A&S 17.7.15/16 to shift to v > 0:
       //
       T k2 = k * k;
       T N = (k2 - v) / (1 - v);
       T Nm1 = (1 - k2) / (1 - v);
       T p2 = sqrt(-v * (k2 - v) / (1 - v));
       T delta = sqrt(1 - k2 * sphi * sphi);
       T result = ellint_pi_imp(N, phi, k, Nm1);

       result *= sqrt(Nm1 * (1 - k2 / N));
       result += ellint_f_imp(phi, k) * k2 / p2;
       result += atan((p2/2) * sin(2 * phi) / delta);
       result /= sqrt((1 - v) * (1 - k2 / v));
       return result;
    }
#if 0  // disabled but retained for future reference: see below.
    if(v > 1)
    {
       //
       // If v > 1 we can use the identity in A&S 17.7.7/8
       // to shift to 0 <= v <= 1.  Unfortunately this
       // identity appears only to function correctly when
       // 0 <= phi <= PI/2, but it's when phi is outside that
       // range that we really need it: That's when
       // Carlson's formula fails, and the periodicity
       // reduction used below on phi doesn't work when v > 1.
       // So we're stuck... the code is archived here in case
       // some bright spart can figure out the fix.
       //
       T k2 = k * k;
       T N = k2 / v;
       T Nm1 = (v - k2) / v;
       T p1 = sqrt((-vc) * (1 - k2 / v));
       T delta = sqrt(1 - k2 * sphi * sphi);
       T result = -ellint_pi_imp(N, phi, k, Nm1);
       result += ellint_f_imp(phi, k);
       result += log((delta + p1 * tan(phi)) / (delta - p1 * tan(phi))) / (2 * p1);
       return result;
    }
#endif

    // Carlson's algorithm works only for |phi| <= pi/2,
    // use the integrand's periodicity to normalize phi
    //
    // Xiaogang's original code used a cast to long long here
    // but that fails if T has more digits than a long long,
    // so rewritten to use fmod instead:
    //
    if(fabs(phi) > 1 / tools::epsilon<T>())
    {
       if(v > 1)
          return tools::domain_error<T>(
            BOOST_CURRENT_FUNCTION,
            "Got v = %1%, but this is only supported for 0 <= phi <= pi/2", v);
       //  
       // Phi is so large that phi%pi is necessarily zero (or garbage),
       // just return the second part of the duplication formula:
       //
       value = 2 * fabs(phi) * ellint_pi_imp(v, k, vc) / constants::pi<T>();
    }
    else
    {
       T rphi = fmod(fabs(phi), constants::pi<T>() / 2);
       T m = 2 * (fabs(phi) - rphi) / constants::pi<T>();
       int sign = 1;
       if(fmod(m, T(2)) > 0.5)
       {
          m += 1;
          sign = -1;
          rphi = constants::pi<T>() / 2 - rphi;
       }

       if((m > 0) && (v > 1))
          return tools::domain_error<T>(
            BOOST_CURRENT_FUNCTION,
            "Got v = %1%, but this is only supported for 0 <= phi <= pi/2", v);

       T sinp = sin(rphi);
       T cosp = cos(rphi);
       x = cosp * cosp;
       t = sinp * sinp;
       y = 1 - k * k * t;
       z = 1;
       if(v * t < 0.5)
           p = 1 - v * t;
       else
           p = x + vc * t;
       value = sign * sinp * (ellint_rf_imp(x, y, z) + v * t * ellint_rj_imp(x, y, z, p) / 3);
       if(m > 0)
         value += m * ellint_pi_imp(v, k, vc);
    }

    if (phi < 0)
    {
        value = -value;    // odd function
    }
    return value;
}

// Complete elliptic integral (Legendre form) of the third kind
template <typename T>
T ellint_pi_imp(T v, T k, T vc)
{
    // Note arg vc = 1-v, possibly without cancellation errors
    using namespace std;
    using namespace boost::math::tools;

    if (abs(k) >= 1)
    {
        return domain_error<T>(BOOST_CURRENT_FUNCTION,
            "Got k = %1%, function requires |k| <= 1", k);
    }
    if(vc <= 0)
    {
       // Result is complex:
        return domain_error<T>(BOOST_CURRENT_FUNCTION,
            "Got v = %1%, function requires v < 1", v);
    }

    if(v == 0)
    {
       return (k == 0) ? boost::math::constants::pi<T>() / 2 : ellint_k_imp(k);
    }

    if(v < 0)
    {
       T k2 = k * k;
       T N = (k2 - v) / (1 - v);
       T Nm1 = (1 - k2) / (1 - v);
       T p2 = sqrt(-v * (k2 - v) / (1 - v));

       T result = boost::math::detail::ellint_pi_imp(N, k, Nm1);

       result *= sqrt(Nm1 * (1 - k2 / N));
       result += ellint_k_imp(k) * k2 / p2;
       result /= sqrt((1 - v) * (1 - k2 / v));
       return result;
    }

    T x = 0;
    T y = 1 - k * k;
    T z = 1;
    T p = vc;
    T value = ellint_rf_imp(x, y, z) + v * ellint_rj_imp(x, y, z, p) / 3;

    return value;
}

} // namespace detail

template <typename T>
inline T ellint_3(T v, T phi, T k)
{
   typedef typename tools::evaluation<typename remove_cv<T>::type>::type value_type;
   return tools::checked_narrowing_cast<typename remove_cv<T>::type>(
      detail::ellint_pi_imp(
         static_cast<value_type>(v), 
         static_cast<value_type>(phi), 
         static_cast<value_type>(k),
         static_cast<value_type>(1-v)), BOOST_CURRENT_FUNCTION);
}

template <typename T>
inline T ellint_3(T v, T k)
{
   typedef typename tools::evaluation<typename remove_cv<T>::type>::type value_type;
   return tools::checked_narrowing_cast<typename remove_cv<T>::type>(
      detail::ellint_pi_imp(
         static_cast<value_type>(v), 
         static_cast<value_type>(k),
         static_cast<value_type>(1-v)), BOOST_CURRENT_FUNCTION);
}

}} // namespaces

#endif // BOOST_MATH_ELLINT_3_HPP
