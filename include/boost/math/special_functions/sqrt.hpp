//  (C) Copyright Matt Borland 2021.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
//  Constexpr implementation of sqrt function

#ifndef BOOST_MATH_SPECIAL_FUNCTIONS_SQRT
#define BOOST_MATH_SPECIAL_FUNCTIONS_SQRT

#include <limits>
#include <type_traits>

namespace boost { namespace math { 

namespace detail {

template <typename Real>
inline constexpr Real sqrt_impl_2(Real x, Real s, Real s2)
{
    return !(s < s2) ? s2 : sqrt_impl_2(x, (x / s + s) / 2, s);
}

template <typename Real>
inline constexpr Real sqrt_impl_1(Real x, Real s)
{
    return sqrt_impl_2(x, (x / s + s) / 2, s);
}

template <typename Real>
inline constexpr Real sqrt_impl(Real x)
{
    return sqrt_impl_1(x, x > 1 ? x : Real(1));
}

// std::isnan is not constexpr according to the standard
template <typename Real>
inline constexpr bool is_nan(Real x)
{
    return x != x;
}

template <typename Real>
inline constexpr bool is_inf(Real x)
{
    return x == std::numeric_limits<Real>::infinity() || -x == std::numeric_limits<Real>::infinity();
}

} // namespace detail

template <typename Real, typename std::enable_if<std::is_floating_point<Real>::value, bool>::type = true>
inline constexpr Real sqrt(Real x)
{
    return detail::is_nan(x) ? std::numeric_limits<Real>::quiet_NaN() : 
           detail::is_inf(x) ? std::numeric_limits<Real>::infinity() : 
           detail::sqrt_impl<Real>(x);
}

template <typename Z, typename std::enable_if<std::is_integral<Z>::value, bool>::type = true>
inline constexpr double sqrt(Z x)
{
    return detail::sqrt_impl<double>(static_cast<double>(x));
}

} // namespace math
} // namespace boost

#endif // BOOST_MATH_SPECIAL_FUNCTIONS_SQRT
