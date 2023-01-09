//  (C) Copyright John Maddock 2005-2021.
//  (C) Copyright Matt Borland 2021.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_CCMATH_HYPOT_HPP
#define BOOST_MATH_CCMATH_HYPOT_HPP

#include <cmath>
#include <array>
#include <limits>
#include <type_traits>
#include <boost/math/tools/config.hpp>
#include <boost/math/tools/is_constant_evaluated.hpp>
#include <boost/math/tools/promotion.hpp>
#include <boost/math/ccmath/sqrt.hpp>
#include <boost/math/ccmath/abs.hpp>
#include <boost/math/ccmath/isinf.hpp>
#include <boost/math/ccmath/isnan.hpp>
#include <boost/math/ccmath/fmax.hpp>
#include <boost/math/ccmath/detail/swap.hpp>

namespace boost::math::ccmath {

namespace detail {

template <typename T>
constexpr T hypot_impl(T x, T y) noexcept
{
    x = boost::math::ccmath::abs(x);
    y = boost::math::ccmath::abs(y);

    if (y > x)
    {
        boost::math::ccmath::detail::swap(x, y);
    }

    if(x * std::numeric_limits<T>::epsilon() >= y)
    {
        return x;
    }

    T rat = y / x;
    return x * boost::math::ccmath::sqrt(1 + rat * rat);
}

template <typename T>
constexpr T hypot_impl(T x, T y, T z) noexcept
{
    x = boost::math::ccmath::abs(x);
    y = boost::math::ccmath::abs(y);
    z = boost::math::ccmath::abs(z);

    T a = boost::math::ccmath::fmax(boost::math::ccmath::fmax(x, y), z);

    return a * boost::math::ccmath::sqrt((x / a) * (x / a) 
                                       + (y / a) * (y / a) 
                                       + (z / a) * (z / a));
}

} // Namespace detail

template <typename Real, std::enable_if_t<!std::is_integral_v<Real>, bool> = true>
constexpr Real hypot(Real x, Real y) noexcept
{
    if(BOOST_MATH_IS_CONSTANT_EVALUATED(x))
    {
        if (boost::math::ccmath::abs(x) == static_cast<Real>(0))
        {
            return boost::math::ccmath::abs(y);
        }
        else if (boost::math::ccmath::abs(y) == static_cast<Real>(0))
        {
            return boost::math::ccmath::abs(x);
        }
        // Return +inf even if the other argument is NaN
        else if (boost::math::ccmath::isinf(x) || boost::math::ccmath::isinf(y))
        {
            return std::numeric_limits<Real>::infinity();
        }
        else if (boost::math::ccmath::isnan(x))
        {
            return x;
        }
        else if (boost::math::ccmath::isnan(y))
        {
            return y;
        }
        
        return boost::math::ccmath::detail::hypot_impl(x, y);
    }
    else
    {
        using std::hypot;
        return hypot(x, y);
    }
}

template <typename T1, typename T2>
constexpr auto hypot(T1 x, T2 y) noexcept
{
    if(BOOST_MATH_IS_CONSTANT_EVALUATED(x))
    {
        using promoted_type = boost::math::tools::promote_args_2_t<T1, T2>;
        return boost::math::ccmath::hypot(static_cast<promoted_type>(x), static_cast<promoted_type>(y));
    }
    else
    {
        using std::hypot;
        return hypot(x, y);
    }
}

constexpr float hypotf(float x, float y) noexcept
{
    return boost::math::ccmath::hypot(x, y);
}

#ifndef BOOST_MATH_NO_LONG_DOUBLE_MATH_FUNCTIONS
constexpr long double hypotl(long double x, long double y) noexcept
{
    return boost::math::ccmath::hypot(x, y);
}
#endif

template <typename Real, std::enable_if_t<!std::is_integral_v<Real>, bool> = true>
constexpr auto hypot(Real x, Real y, Real z) noexcept
{
    if (BOOST_MATH_IS_CONSTANT_EVALUATED(x))
    {
        if (boost::math::ccmath::isinf(x) || boost::math::ccmath::isinf(y) || boost::math::ccmath::isinf(z))
        {
            return std::numeric_limits<Real>::infinity();
        }
        else if (boost::math::ccmath::isnan(x))
        {
            return x;
        }
        else if (boost::math::ccmath::isnan(y))
        {
            return y;
        }
        else if (boost::math::ccmath::isnan(z))
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

template <typename T1, typename T2, typename T3>
constexpr auto hypot(T1 x, T2 y, T3 z) noexcept
{
    if (BOOST_MATH_IS_CONSTANT_EVALUATED(x))
    {
        using promoted_type = tools::promote_args_t<T1, T2, T3>;
        return boost::math::ccmath::hypot(static_cast<promoted_type>(x), 
                                          static_cast<promoted_type>(y), 
                                          static_cast<promoted_type>(z));
    }
    else
    {
        using std::hypot;
        return hypot(x, y, z);
    }
}

} // Namespaces

#endif // BOOST_MATH_CCMATH_HYPOT_HPP
