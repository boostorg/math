//  (C) Copyright Matt Borland 2021 - 2022.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_CCMATH_REMAINDER_HPP
#define BOOST_MATH_CCMATH_REMAINDER_HPP

#include <boost/math/ccmath/detail/config.hpp>

#ifdef BOOST_MATH_NO_CCMATH
#error "The header <boost/math/remainder.hpp> can only be used in C++17 and later."
#endif

#ifndef BOOST_MATH_BUILD_MODULE
#include <cstdint>
#endif
#include <boost/math/tools/promotion.hpp>
#include <boost/math/ccmath/abs.hpp>
#include <boost/math/ccmath/isinf.hpp>
#include <boost/math/ccmath/isnan.hpp>
#include <boost/math/ccmath/isfinite.hpp>
#include <boost/math/ccmath/modf.hpp>

BOOST_MATH_NAMESPACE_BEGIN namespace ccmath {

namespace detail {

template <typename T>
constexpr T remainder_impl(const T x, const T y)
{
    T n = 0;

    if (T fractional_part = BOOST_MATH_NAMESPACE::ccmath::modf((x / y), &n); fractional_part > static_cast<T>(1.0/2))
    {
        ++n;
    }
    else if (fractional_part < static_cast<T>(-1.0/2))
    {
        --n;
    }

    return x - n*y;
}

} // Namespace detail

BOOST_MATH_EXPORT template <typename Real, std::enable_if_t<!std::is_integral_v<Real>, bool> = true>
constexpr Real remainder(Real x, Real y)
{
    if (BOOST_MATH_IS_CONSTANT_EVALUATED(x))
    {
        if (BOOST_MATH_NAMESPACE::ccmath::isinf(x) && !BOOST_MATH_NAMESPACE::ccmath::isnan(y))
        {
            return std::numeric_limits<Real>::quiet_NaN();
        }
        else if (BOOST_MATH_NAMESPACE::ccmath::abs(y) == static_cast<Real>(0) && !BOOST_MATH_NAMESPACE::ccmath::isnan(x))
        {
            return std::numeric_limits<Real>::quiet_NaN();
        }
        else if (BOOST_MATH_NAMESPACE::ccmath::isnan(x))
        {
            return x;
        }
        else if (BOOST_MATH_NAMESPACE::ccmath::isnan(y))
        {
            return y;
        }

        return BOOST_MATH_NAMESPACE::ccmath::detail::remainder_impl(x, y);
    }
    else
    {
        using std::remainder;
        return remainder(x, y);
    }
}

BOOST_MATH_EXPORT template <typename T1, typename T2>
constexpr auto remainder(T1 x, T2 y)
{
    if (BOOST_MATH_IS_CONSTANT_EVALUATED(x))
    {
        using promoted_type = BOOST_MATH_NAMESPACE::tools::promote_args_t<T1, T2>;
        return BOOST_MATH_NAMESPACE::ccmath::remainder(promoted_type(x), promoted_type(y));
    }
    else
    {
        using std::remainder;
        return remainder(x, y);
    }
}

constexpr float remainderf(float x, float y)
{
    return BOOST_MATH_NAMESPACE::ccmath::remainder(x, y);
}

#ifndef BOOST_MATH_NO_LONG_DOUBLE_MATH_FUNCTIONS
constexpr long double remainderl(long double x, long double y)
{
    return BOOST_MATH_NAMESPACE::ccmath::remainder(x, y);
}
#endif

} BOOST_MATH_NAMESPACE_END

#endif // BOOST_MATH_CCMATH_REMAINDER_HPP
