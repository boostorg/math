//  (C) Copyright Matt Borland 2021.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_CCMATH_ILOGB_HPP
#define BOOST_MATH_CCMATH_ILOGB_HPP

#include <boost/math/ccmath/detail/config.hpp>

#ifdef BOOST_MATH_NO_CCMATH
#error "The header <boost/math/ilogb.hpp> can only be used in C++17 and later."
#endif

#include <boost/math/ccmath/logb.hpp>
#include <boost/math/ccmath/isinf.hpp>
#include <boost/math/ccmath/isnan.hpp>
#include <boost/math/ccmath/abs.hpp>

BOOST_MATH_NAMESPACE_BEGIN namespace ccmath {

// If arg is not zero, infinite, or NaN, the value returned is exactly equivalent to static_cast<int>(std::logb(arg))
BOOST_MATH_EXPORT template <typename Real, std::enable_if_t<!std::is_integral_v<Real>, bool> = true>
inline constexpr int ilogb(Real arg) noexcept
{
    if(BOOST_MATH_IS_CONSTANT_EVALUATED(arg))
    {
        return BOOST_MATH_NAMESPACE::ccmath::abs(arg) == Real(0) ? FP_ILOGB0 :
               BOOST_MATH_NAMESPACE::ccmath::isinf(arg) ? INT_MAX :
               BOOST_MATH_NAMESPACE::ccmath::isnan(arg) ? FP_ILOGBNAN :
               static_cast<int>(BOOST_MATH_NAMESPACE::ccmath::logb(arg));
    }
    else
    {
        using std::ilogb;
        return ilogb(arg);
    }
}

BOOST_MATH_EXPORT template <typename Z, std::enable_if_t<std::is_integral_v<Z>, bool> = true>
inline constexpr int ilogb(Z arg) noexcept
{
    return BOOST_MATH_NAMESPACE::ccmath::ilogb(static_cast<double>(arg));
}

inline constexpr int ilogbf(float arg) noexcept
{
    return BOOST_MATH_NAMESPACE::ccmath::ilogb(arg);
}

#ifndef BOOST_MATH_NO_LONG_DOUBLE_MATH_FUNCTIONS
inline constexpr int ilogbl(long double arg) noexcept
{
    return BOOST_MATH_NAMESPACE::ccmath::ilogb(arg);
}
#endif

} BOOST_MATH_NAMESPACE_END

#endif // BOOST_MATH_CCMATH_ILOGB_HPP
