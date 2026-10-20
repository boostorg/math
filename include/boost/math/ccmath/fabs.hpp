//  (C) Copyright Matt Borland 2021.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
//  Constepxr implementation of fabs (see c.math.abs secion 26.8.2 of the ISO standard)

#ifndef BOOST_MATH_CCMATH_FABS
#define BOOST_MATH_CCMATH_FABS

#include <boost/math/ccmath/detail/config.hpp>

#ifdef BOOST_MATH_NO_CCMATH
#error "The header <boost/math/fabs.hpp> can only be used in C++17 and later."
#endif

#include <boost/math/ccmath/abs.hpp>

BOOST_MATH_NAMESPACE_BEGIN namespace ccmath {

BOOST_MATH_EXPORT template <typename T>
inline constexpr auto fabs(T x) noexcept
{
    return BOOST_MATH_NAMESPACE::ccmath::abs(x);
}

inline constexpr float fabsf(float x) noexcept
{
    return BOOST_MATH_NAMESPACE::ccmath::abs(x);
}

#ifndef BOOST_MATH_NO_LONG_DOUBLE_MATH_FUNCTIONS
inline constexpr long double fabsl(long double x) noexcept
{
    return BOOST_MATH_NAMESPACE::ccmath::abs(x);
}
#endif

} BOOST_MATH_NAMESPACE_END

#endif // BOOST_MATH_CCMATH_FABS
