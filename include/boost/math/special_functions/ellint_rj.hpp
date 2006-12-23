//  Copyright (c) 2006 Xiaogang Zhang
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
//  History:
//  XZ wrote the original of this file as part of the Google
//  Summer of Code 2006.  JM modified it to fit into the
//  Boost.Math conceptual framework better, and to correctly
//  handle the p < 0 case.
//

#ifndef BOOST_MATH_ELLINT_RJ_HPP
#define BOOST_MATH_ELLINT_RJ_HPP

#include <boost/math/special_functions/math_fwd.hpp>
#include <boost/math/tools/config.hpp>
#include <boost/math/tools/error_handling.hpp>
#include <boost/math/tools/evaluation_type.hpp>
#include <boost/math/special_functions/ellint_rc.hpp>

// Carlson's elliptic integral of the third kind
// R_J(x, y, z, p) = 1.5 * \int_{0}^{\infty} (t+p)^{-1} [(t+x)(t+y)(t+z)]^{-1/2} dt
// Carlson, Numerische Mathematik, vol 33, 1 (1979)

namespace boost { namespace math { namespace detail{

template <typename T>
T ellint_rj_imp(T x, T y, T z, T p)
{
    T value, u, lambda, alpha, beta, sigma, factor, tolerance;
    T X, Y, Z, P, EA, EB, EC, E2, E3, S1, S2, S3;
    int k;

    using namespace std;
    using namespace boost::math::tools;

    if (x < 0)
    {
        return domain_error<T>(BOOST_CURRENT_FUNCTION,
            "Argument x must be non-negative, but got x = %1%", x);
    }
    if(y < 0)
    {
        return domain_error<T>(BOOST_CURRENT_FUNCTION,
            "Argument y must be non-negative, but got y = %1%", y);
    }
    if(z < 0)
    {
        return domain_error<T>(BOOST_CURRENT_FUNCTION,
            "Argument z must be non-negative, but got z = %1%", z);
    }
    if(p == 0)
    {
        return domain_error<T>(BOOST_CURRENT_FUNCTION,
            "Argument p must not be zero, but got p = %1%", p);
    }
    if (x + y == 0 || y + z == 0 || z + x == 0)
    {
        return domain_error<T>(BOOST_CURRENT_FUNCTION,
            "At most one argument can be zero, "
            "only possible result is %1%.", std::numeric_limits<T>::quiet_NaN());
    }

    // error scales as the 6th power of tolerance
    tolerance = pow(T(1) * tools::epsilon<T>() / 3, T(1) / 6);

    // for p < 0, the integral is singular, return Cauchy principal value
    if (p < 0)
    {
       T q = -p;
       T pmy = (z - y) * (y - x) / (y + q);  // p - y
       if(pmy < 0)
       {
          //
          // TODO, FIXME if you can.... (JM Dec. 2006)
          //
          // The logic breaks down here, we can go into
          // an infinite recursion unless we bail out right away!!
          //
          // This may be fixable by permuting x, y, and z, but may not
          // be worth the hassle, fix this if you care about this use case!
          //
          return tools::domain_error<T>(
             BOOST_CURRENT_FUNCTION,
             "Unable to compute Cauchy principle value, p had the value %1% "
             "and (z-y)(y-x) was negative.  Identity formula could not be applied!",
             p);
       }
       T p = pmy + y;
       value = ellint_rj(x, y, z, p);
       value *= pmy;
       value -= 3 * ellint_rf(x, y, z);
       value += 3 * sqrt((x * y * z) / (x * z + p * q)) * ellint_rc(x * z + p * q, p * q);
       value /= (y + q);
       return value;
    }

    // duplication
    sigma = 0;
    factor = 1;
    for (k = 1; k < BOOST_MATH_MAX_ITER; k++)
    {
        u = (x + y + z + p + p) / 5;
        X = (u - x) / u;
        Y = (u - y) / u;
        Z = (u - z) / u;
        P = (u - p) / u;
        
        if ((tools::max)(abs(X), abs(Y), abs(Z), abs(P)) < tolerance) 
           break;

        T sx = sqrt(x);
        T sy = sqrt(y);
        T sz = sqrt(z);
        
        lambda = sy * (sx + sz) + sz * sx;
        alpha = p * (sx + sy + sz) + sx * sy * sz;
        alpha *= alpha;
        beta = p * (p + lambda) * (p + lambda);
        sigma += factor * ellint_rc(alpha, beta);
        factor /= 4;
        x = (x + lambda) / 4;
        y = (y + lambda) / 4;
        z = (z + lambda) / 4;
        p = (p + lambda) / 4;
    }
    // Check to see if we gave up too soon:
    tools::check_series_iterations(BOOST_CURRENT_FUNCTION, k);

    // Taylor series expansion to the 5th order
    EA = X * Y + Y * Z + Z * X;
    EB = X * Y * Z;
    EC = P * P;
    E2 = EA - 3 * EC;
    E3 = EB + 2 * P * (EA - EC);
    S1 = 1 + E2 * (E2 * T(9) / 88 - E3 * T(9) / 52 - T(3) / 14);
    S2 = EB * (T(1) / 6 + P * (T(-6) / 22 + P * T(3) / 26));
    S3 = P * ((EA - EC) / 3 - P * EA * T(3) / 22);
    value = 3 * sigma + factor * (S1 + S2 + S3) / (u * sqrt(u));

    return value;
}

} // namespace detail

template <typename T>
inline T ellint_rj(T x, T y, T z, T p)
{
   typedef typename tools::evaluation<typename remove_cv<T>::type>::type value_type;
   return tools::checked_narrowing_cast<typename remove_cv<T>::type>(
      detail::ellint_rj_imp(
         static_cast<value_type>(x),
         static_cast<value_type>(y),
         static_cast<value_type>(z),
         static_cast<value_type>(p)), BOOST_CURRENT_FUNCTION);
}

}} // namespaces

#endif // BOOST_MATH_ELLINT_RJ_HPP
