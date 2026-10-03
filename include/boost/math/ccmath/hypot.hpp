//  (C) Copyright John Maddock 2005-2021.
//  (C) Copyright Matt Borland 2021.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_CCMATH_HYPOT_HPP
#define BOOST_MATH_CCMATH_HYPOT_HPP

#include <boost/math/ccmath/detail/config.hpp>

#ifdef BOOST_MATH_NO_CCMATH
#error "The header <boost/math/hypot.hpp> can only be used in C++17 and later."
#endif

#ifndef BOOST_MATH_BUILD_MODULE
#include <array>
#endif
#include <boost/math/tools/config.hpp>
#include <boost/math/tools/promotion.hpp>
#include <boost/math/ccmath/sqrt.hpp>
#include <boost/math/ccmath/abs.hpp>
#include <boost/math/ccmath/isinf.hpp>
#include <boost/math/ccmath/isnan.hpp>
#include <boost/math/ccmath/fmax.hpp>
#include <boost/math/ccmath/detail/swap.hpp>

BOOST_MATH_NAMESPACE_BEGIN namespace ccmath {

namespace detail {

template <typename T>
constexpr T hypot_impl(T x, T y) noexcept
{
    x = BOOST_MATH_NAMESPACE::ccmath::abs(x);
    y = BOOST_MATH_NAMESPACE::ccmath::abs(y);

    if (y > x)
    {
        BOOST_MATH_NAMESPACE::ccmath::detail::swap(x, y);
    }

    if(x * std::numeric_limits<T>::epsilon() >= y)
    {
        return x;
    }

    T rat = y / x;
    return x * BOOST_MATH_NAMESPACE::ccmath::sqrt(1 + rat * rat);
}

template <typename T>
constexpr T hypot_impl(T x, T y, T z) noexcept
{
    x = BOOST_MATH_NAMESPACE::ccmath::abs(x);
    y = BOOST_MATH_NAMESPACE::ccmath::abs(y);
    z = BOOST_MATH_NAMESPACE::ccmath::abs(z);

    T a = BOOST_MATH_NAMESPACE::ccmath::fmax(BOOST_MATH_NAMESPACE::ccmath::fmax(x, y), z);
    if (a == 0)
    {
        return a;
    }

    return a * BOOST_MATH_NAMESPACE::ccmath::sqrt((x / a) * (x / a) 
                                       + (y / a) * (y / a) 
                                       + (z / a) * (z / a));
}

} // Namespace detail

BOOST_MATH_EXPORT template <typename Real, std::enable_if_t<!std::is_integral_v<Real>, bool> = true>
constexpr Real hypot(Real x, Real y) noexcept
{
    if(BOOST_MATH_IS_CONSTANT_EVALUATED(x))
    {
        if (BOOST_MATH_NAMESPACE::ccmath::abs(x) == static_cast<Real>(0))
        {
            return BOOST_MATH_NAMESPACE::ccmath::abs(y);
        }
        else if (BOOST_MATH_NAMESPACE::ccmath::abs(y) == static_cast<Real>(0))
        {
            return BOOST_MATH_NAMESPACE::ccmath::abs(x);
        }
        // Return +inf even if the other argument is NaN
        else if (BOOST_MATH_NAMESPACE::ccmath::isinf(x) || BOOST_MATH_NAMESPACE::ccmath::isinf(y))
        {
            return std::numeric_limits<Real>::infinity();
        }
        else if (BOOST_MATH_NAMESPACE::ccmath::isnan(x))
        {
            return x;
        }
        else if (BOOST_MATH_NAMESPACE::ccmath::isnan(y))
        {
            return y;
        }
        
        return BOOST_MATH_NAMESPACE::ccmath::detail::hypot_impl(x, y);
    }
    else
    {
        using std::hypot;
        return hypot(x, y);
    }
}

BOOST_MATH_EXPORT template <typename T1, typename T2>
constexpr auto hypot(T1 x, T2 y) noexcept
{
    if(BOOST_MATH_IS_CONSTANT_EVALUATED(x))
    {
        using promoted_type = BOOST_MATH_NAMESPACE::tools::promote_args_t<T1, T2>;
        return BOOST_MATH_NAMESPACE::ccmath::hypot(static_cast<promoted_type>(x), static_cast<promoted_type>(y));
    }
    else
    {
        using std::hypot;
        return hypot(x, y);
    }
}

constexpr float hypotf(float x, float y) noexcept
{
    return BOOST_MATH_NAMESPACE::ccmath::hypot(x, y);
}

#ifndef BOOST_MATH_NO_LONG_DOUBLE_MATH_FUNCTIONS
constexpr long double hypotl(long double x, long double y) noexcept
{
    return BOOST_MATH_NAMESPACE::ccmath::hypot(x, y);
}
#endif

BOOST_MATH_EXPORT template <typename Real, std::enable_if_t<!std::is_integral_v<Real>, bool> = true>
constexpr Real hypot(Real x, Real y, Real z) noexcept
{
    if (BOOST_MATH_IS_CONSTANT_EVALUATED(x))
    {
        if (BOOST_MATH_NAMESPACE::ccmath::isinf(x) || BOOST_MATH_NAMESPACE::ccmath::isinf(y) || BOOST_MATH_NAMESPACE::ccmath::isinf(z))
        {
            return std::numeric_limits<Real>::infinity();
        }
        else if (BOOST_MATH_NAMESPACE::ccmath::isnan(x))
        {
            return x;
        }
        else if (BOOST_MATH_NAMESPACE::ccmath::isnan(y))
        {
            return y;
        }
        else if (BOOST_MATH_NAMESPACE::ccmath::isnan(z))
        {
            return z;
        }

        return detail::hypot_impl(x, y, z);
    }
    else
    {
        using std::hypot;
        return hypot(x, y, z);
    }
}

BOOST_MATH_EXPORT template <typename T1, typename T2, typename T3>
constexpr auto hypot(T1 x, T2 y, T3 z) noexcept
{
    if (BOOST_MATH_IS_CONSTANT_EVALUATED(x))
    {
        using promoted_type = tools::promote_args_t<T1, T2, T3>;
        return BOOST_MATH_NAMESPACE::ccmath::hypot(static_cast<promoted_type>(x), 
                                          static_cast<promoted_type>(y), 
                                          static_cast<promoted_type>(z));
    }
    else
    {
        using std::hypot;
        return hypot(x, y, z);
    }
}

} BOOST_MATH_NAMESPACE_END

#endif // BOOST_MATH_CCMATH_HYPOT_HPP
