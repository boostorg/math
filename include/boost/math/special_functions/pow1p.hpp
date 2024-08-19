//  (C) Copyright Matt Borland 2024.
//  (C) Copyrigh Fancidev 2024.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_SF_POW1P_HPP
#define BOOST_MATH_SF_POW1P_HPP

#include <boost/math/tools/config.hpp>
#include <boost/math/tools/numeric_limits.hpp>
#include <boost/math/tools/promotion.hpp>
#include <boost/math/special_functions/fpclassify.hpp>
#include <boost/math/policies/error_handling.hpp>

// For cuda we would rather use builtin nextafter than unsupported boost::math::nextafter
// NVRTC does not support the forward declarations header
#ifndef BOOST_MATH_ENABLE_CUDA
#  include <boost/math/special_functions/next.hpp>
#  ifndef BOOST_MATH_HAS_NVRTC
#    include <boost/math/special_functions/math_fwd.hpp>
#  endif // BOOST_MATH_HAS_NVRTC
#endif // BOOST_MATH_ENABLE_CUDA

namespace boost {
namespace math {

namespace detail {

template <typename T, typename Policy>
BOOST_MATH_GPU_ENABLED T pow1p_imp(const T x, const T y, const Policy& pol)
{
    BOOST_MATH_STD_USING
    constexpr auto function = "boost::math::pow1p<%1%>(%1%, %1%)";

    // The special values follow the spec of the pow() function defined in
    // IEEE-754, Section 9.2.1.  The order of the `if` statements matters.
    if (y == T(0)) 
    { 
        // pow(x, +/-0)
        return T(1);
    }
    
    if (x == T(-1))
    { 
        // pow(+/-0, y)
        if (boost::math::isfinite(y) && y < 0) 
        {
            return boost::math::policies::raise_domain_error<T>(function, "Division by 0", x, pol);
        }

        // Gets correct special handling
        return pow(T(0), y);
    }

    if (x == T(-2) && boost::math::isinf(y)) 
    { 
        // pow(-1, +/-inf)
        return T(1);
    }

    if (x == T(0)) { // pow(+1, y)
        return T(1);
    }

    if (boost::math::isinf(y))
    {
        if (boost::math::signbit(y))
        {
            // Pow(y, -inf)
            return (x < 0 && x > -2) ? boost::math::numeric_limits<T>::infinity() : T(0);
        }
        else
        {
            // pow(x, +inf)
            return (x < 0 && x > -2) ? T(0) : boost::math::numeric_limits<T>::infinity();
        }
    }

    if (boost::math::isinf(x)) 
    { 
        // pow(+/-inf, y)
        return pow(x, y);
    }

    // Up to this point, (1+x) = {+/-0, +/-1, +/-inf} have been handled, and
    // and y = {+/-0, +/-inf} have been handled.  Next, we handle `nan` and
    // y = {+/-1}.
    if (boost::math::isnan(x)) 
    {
        return x;
    }
    else if (boost::math::isnan(y))
    {
        return y;
    }

    if (y == T(1)) 
    {
        return T(1) + x;
    }
    if (y == T(-1)) 
    {
        return T(1) / (T(1) + x); // guaranteed no overflow
    }

    // Handle (1+x) < 0
    if (x < -1) 
    {
        if (fmod(y, T(1)) != T(0))
        { 
            // y is not an integer so we have an error
            return boost::math::policies::raise_evaluation_error<T>(function, "Y is not an integer", y, pol);
        }

        // TODO(fancidev): maybe use (1+x)^y == [(1+x)^2]^(y/2)
        return pow(1+x, y);
    }

    // Now x, y are finite and not equal to 0 or +/-1, and x > -1.
    // To compute (1+x)^y, write (1+x) == (s+t) where s is equal to (1+x)
    // rounded toward 1, and t is the (exact) rounding error.
    T s, t;
    if (x < 0) 
    {
        s = T(1) + x;
        t = x - (s - T(1));
        if (t > 0) 
        {
            #ifdef BOOST_MATH_ENABLE_CUDA
            s = ::nextafter(s, T(1));
            #else
            s = boost::math::nextafter(s, T(1));
            #endif

            t = x - (s - T(1));
        }
    } 
    else if (x < 1) 
    {
        s = T(1) + x;
        t = x - (s - T(1));
        if (t < 0) 
        {
            #ifdef BOOST_MATH_ENABLE_CUDA
            s = ::nextafter(s, T(0));
            #else
            s = boost::math::nextafter(s, T(0));
            #endif

            t = x - (s - T(1));
        }
    } 
    else
    {
        s = x + T(1);
        t = T(1) - (s - x);
        if (t < 0) 
        {
            #ifdef BOOST_MATH_ENABLE_CUDA
            s = ::nextafter(s, T(0));
            #else
            s = boost::math::nextafter(s, T(0));
            #endif

            t = T(1) - (s - x);
        }
    }

    // Because x > -1 and s is rounded toward 1, s is guaranteed to be > 0.
    // Then (1+x)^y == (s+t)^y == (s^y)*((1+u)^y), where u := t / s.
    // It can be shown that either both terms <= 1 or both >= 1, so
    // if the first term over/underflows, then the result over/underflows.
    T u = t / s;
    T term1 = pow(s, y);
    if (term1 == T(0) || boost::math::isinf(term1))
    {
        return term1;
    }

    // (1+u)^y == exp(y*log(1+u)).  Since u is close to machine epsilon,
    // log(1+u) ~= u.  Let y*u == z+w, where z is the rounded result and
    // w is the rounding error.  This improves accuracy when y is large.
    // Then exp(y*u) == exp(z)*exp(w).
    T z = y * u;
    T w = fma(y, u, -z);
    T term2 = exp(z) * exp(w);
    
    return term1 * term2;
}

} // Namespace detail

template <typename T1, typename T2, typename Policy>
BOOST_MATH_GPU_ENABLED tools::promote_args_t<T1, T2>
pow1p(const T1 x, const T2 y, const Policy& pol)
{
    using result_type = tools::promote_args_t<T1, T2>;
    return detail::pow1p_imp(static_cast<result_type>(x), static_cast<result_type>(y), pol);
}

template <typename T1, typename T2>
BOOST_MATH_GPU_ENABLED tools::promote_args_t<T1, T2>
pow1p(const T1 x, const T2 y)
{
    using result_type = tools::promote_args_t<T1, T2>;
    return detail::pow1p_imp(static_cast<result_type>(x), static_cast<result_type>(y), policies::policy<>());
}

} // Namespace math
} // Namespace boost

#endif // BOOST_MATH_SF_POW1P_HPP
