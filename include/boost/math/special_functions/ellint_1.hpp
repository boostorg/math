//  Copyright (c) 2006 Xiaogang Zhang
//  Copyright (c) 2006 John Maddock
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
//  History:
//  XZ wrote the original of this file as part of the Google
//  Summer of Code 2006.  JM modified it to fit into the
//  Boost.Math conceptual framework better, and to ensure
//  that the code continues to work no matter how many digits
//  type T has.

#ifndef BOOST_MATH_ELLINT_1_HPP
#define BOOST_MATH_ELLINT_1_HPP

#include <boost/math/special_functions/ellint_rf.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/math/tools/error_handling.hpp>
#include <boost/math/tools/evaluation_type.hpp>

// Elliptic integrals (complete and incomplete) of the first kind
// Carlson, Numerische Mathematik, vol 33, 1 (1979)

namespace boost { namespace math {

namespace detail{

// Elliptic integral (Legendre form) of the first kind
template <typename T>
T ellint_f_imp(T phi, T k)
{
    using namespace std;
    using namespace boost::math::tools;
    using namespace boost::math::constants;

    if (abs(k) > 1)
    {
        return domain_error<T>(BOOST_CURRENT_FUNCTION,
            "Got k = %1%, function requires |k| <= 1", k);
    }

    bool invert = false;
    if(phi < 0)
    {
       phi = fabs(phi);
       invert = true;
    }

    T result;

    if(phi >= tools::max_value<T>())
    {
       // Need to handle infinity as a special case:
       result = tools::overflow_error<T>(BOOST_CURRENT_FUNCTION, 0);
    }
    else if(phi > 1 / tools::epsilon<T>())
    {
       // Phi is so large that phi%pi is necessarily zero (or garbage),
       // just return the second part of the duplication formula:
       result = 2 * phi * ellint_k_imp(k) / constants::pi<T>();
    }
    else
    {
       // Carlson's algorithm works only for |phi| <= pi/2,
       // use the integrand's periodicity to normalize phi
       //
       // Xiaogang's original code used a cast to long long here
       // but that fails if T has more digits than a long long,
       // so rewritten to use fmod instead:
       //
       T rphi = fmod(phi, constants::pi<T>() / 2);
       T m = 2 * (phi - rphi) / constants::pi<T>();
       int s = 1;
       if(fmod(m, T(2)) > 0.5)
       {
          m += 1;
          s = -1;
          rphi = constants::pi<T>() / 2 - rphi;
       }
       T sinp = sin(rphi);
       T cosp = cos(rphi);
       result = s * sinp * ellint_rf_imp(cosp * cosp, 1 - k * k * sinp * sinp, T(1));
       if(m != 0)
          result += m * ellint_k_imp(k);
    }
    return invert ? -result : result;
}

// Complete elliptic integral (Legendre form) of the first kind
template <typename T>
T ellint_k_imp(T k)
{
    using namespace std;
    using namespace boost::math::tools;

    if (abs(k) > 1)
    {
        return domain_error<T>(BOOST_CURRENT_FUNCTION,
            "Got k = %1%, function requires |k| <= 1", k);
    }
    if (abs(k) == 1)
    {
        return overflow_error<T>(BOOST_CURRENT_FUNCTION, 0);
    }

    T x = 0;
    T y = 1 - k * k;
    T z = 1;
    T value = ellint_rf_imp(x, y, z);

    return value;
}

}

// Complete elliptic integral (Legendre form) of the first kind
template <typename T>
inline T ellint_1(T k)
{
   typedef typename tools::evaluation<typename remove_cv<T>::type>::type value_type;
   return tools::checked_narrowing_cast<typename remove_cv<T>::type>(detail::ellint_k_imp(static_cast<value_type>(k)), BOOST_CURRENT_FUNCTION);
}

// Elliptic integral (Legendre form) of the first kind
template <typename T>
inline T ellint_1(T k, T phi)
{
   typedef typename tools::evaluation<typename remove_cv<T>::type>::type value_type;
   return tools::checked_narrowing_cast<typename remove_cv<T>::type>(detail::ellint_f_imp(static_cast<value_type>(phi), static_cast<value_type>(k)), BOOST_CURRENT_FUNCTION);
}

}} // namespaces

#endif // BOOST_MATH_ELLINT_1_HPP
