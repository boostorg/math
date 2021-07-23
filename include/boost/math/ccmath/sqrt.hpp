//  (C) Copyright Matt Borland 2021.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
//  Constexpr implementation of sqrt function

#ifndef BOOST_MATH_CCMATH_SQRT
#define BOOST_MATH_CCMATH_SQRT

#include <limits>
#include <type_traits>
#include <boost/math/ccmath/isnan.hpp>
#include <boost/math/ccmath/isinf.hpp>

namespace boost::math::ccmath { 

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

} // namespace detail

template <typename Real, typename std::enable_if<std::is_floating_point<Real>::value, bool>::type = true>
inline constexpr Real sqrt(Real x)
{
    return boost::math::ccmath::isnan(x) ? std::numeric_limits<Real>::quiet_NaN() : 
           boost::math::ccmath::isinf(x) ? std::numeric_limits<Real>::infinity() : 
           detail::sqrt_impl<Real>(x);
}

template <typename Z, typename std::enable_if<std::is_integral<Z>::value, bool>::type = true>
inline constexpr double sqrt(Z x)
{
    return detail::sqrt_impl<double>(static_cast<double>(x));
}

} // Namespaces

#endif // BOOST_MATH_CCMATH_SQRT
