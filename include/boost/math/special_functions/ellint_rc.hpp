//  Copyright (c) 2006 Xiaogang Zhang
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
//  History:
//  XZ wrote the original of this file as part of the Google
//  Summer of Code 2006.  JM modified it to fit into the
//  Boost.Math conceptual framework better, and to correctly
//  handle the y < 0 case.
//

#ifndef BOOST_MATH_ELLINT_RC_HPP
#define BOOST_MATH_ELLINT_RC_HPP

#include <boost/math/tools/error_handling.hpp>
#include <boost/math/tools/config.hpp>
#include <boost/math/tools/evaluation_type.hpp>

// Carlson's degenerate elliptic integral
// R_C(x, y) = R_F(x, y, y) = 0.5 * \int_{0}^{\infty} (t+x)^{-1/2} (t+y)^{-1} dt
// Carlson, Numerische Mathematik, vol 33, 1 (1979)

namespace boost { namespace math { namespace detail{

template <typename T>
T ellint_rc_imp(T x, T y)
{
    T value, S, u, lambda, tolerance, prefix;
    int k;

    using namespace std;
    using namespace boost::math::tools;

    if(x < 0)
    {
        return domain_error<T>(BOOST_CURRENT_FUNCTION,
            "Argument x must be non-negative but got %1%", x);
    }
    if(y == 0)
    {
        return domain_error<T>(BOOST_CURRENT_FUNCTION,
            "Argument y must not be zero but got %1%", y);
    }

    // error scales as the 6th power of tolerance
    tolerance = pow(4 * tools::epsilon<T>(), T(1) / 6);

    // for y < 0, the integral is singular, return Cauchy principal value
    if (y < 0)
    {
        prefix = sqrt(x / (x - y));
        x = x - y;
        y = -y;
    }
    else
       prefix = 1;

    // duplication
    for (k = 1; k < BOOST_MATH_MAX_ITER; k++)
    {
        u = (x + y + y) / 3;
        S = y / u - 1;               // 1 - x / u = 2 * S

        if (2 * abs(S) < tolerance) 
           break;

        T sx = sqrt(x);
        T sy = sqrt(y);
        lambda = 2 * sx * sy + y;
        x = (x + lambda) / 4;
        y = (y + lambda) / 4;
    }
    // Check to see if we gave up too soon:
    tools::check_series_iterations(BOOST_CURRENT_FUNCTION, k);

    // Taylor series expansion to the 5th order
    value = (1 + S * S * (T(3) / 10 + S * (T(1) / 7 + S * (T(3) / 8 + S * T(9) / 22)))) / sqrt(u);

    return value * prefix;
}

} // namespace detail

template <typename T>
inline T ellint_rc(T x, T y)
{
   typedef typename tools::evaluation<typename remove_cv<T>::type>::type value_type;
   return tools::checked_narrowing_cast<typename remove_cv<T>::type>(
      detail::ellint_rc_imp(
         static_cast<value_type>(x),
         static_cast<value_type>(y)), BOOST_CURRENT_FUNCTION);
}

}} // namespaces

#endif // BOOST_MATH_ELLINT_RC_HPP
