//  (C) Copyright Matt Borland 2021.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
//  Constexpr implementation of sqrt function

#include <type_traits>

namespace boost { namespace math { 

namespace detail {

template <typename T>
inline constexpr T sqrt_impl_2(T x, T s, T s2)
{
    return !(s < s2) ? s2 : sqrt_impl_2(x, (x / s + s) / 2, s);
}

template <typename T>
inline constexpr T sqrt_impl_1(T x, T s)
{
    return sqrt_impl_2(x, (x / s + s) / 2, s);
}

template <typename T>
inline constexpr T sqrt_impl(T x)
{
    return sqrt_impl_1(x, x > 1 ? x : T(1));
}
} // namespace detail

template <typename Real, typename std::enable_if<std::is_floating_point<Real>::value, bool>::type = true>
inline constexpr Real sqrt(Real x)
{
    return detail::sqrt_impl<Real>(x);
}

template <typename Z, typename std::enable_if<std::is_integral<Z>::value, bool>::type = true>
inline constexpr double sqrt(Z x)
{
    return detail::sqrt_impl<double>(static_cast<double>(x));
}

} // namespace math
} // namespace boost
