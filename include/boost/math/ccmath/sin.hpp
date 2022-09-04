//  (C) Copyright Matt Borland 2022.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_CCMATH_TRIG_HPP
#define BOOST_MATH_CCMATH_TRIG_HPP

#include <cmath>
#include <type_traits>
#include <boost/math/constants/constants.hpp>
#include <boost/math/tools/is_constant_evaluated.hpp>
#include <boost/math/tools/config.hpp>
#include <boost/math/ccmath/abs.hpp>
#include <boost/math/ccmath/isinf.hpp>
#include <boost/math/ccmath/isnan.hpp>

namespace boost::math::ccmath {

namespace detail {

template <typename T>
[[nodiscard]] constexpr T sin_impl(T x) noexcept
{
    // Precalculated coeffs from the minmax polynomial for up to 128 bit long doubles
    // error not to exceed 2e-16 in the range [-pi/2, pi/2]
    constexpr auto a0 = static_cast<T>(+1L);
    constexpr auto a1 = static_cast<T>(-1.666666666666580809419428987894207e-1L);
    constexpr auto a2 = static_cast<T>(+8.333333333262716094425037738346873e-3L);
    constexpr auto a3 = static_cast<T>(-1.984126982005911439283646346964929e-4L);
    constexpr auto a4 = static_cast<T>(+2.755731607338689220657382272783309e-6L);
    constexpr auto a5 = static_cast<T>(-2.505185130214293595900283001271652e-8L);
    constexpr auto a6 = static_cast<T>(+1.604729591825977403374012010065495e-10L);
    constexpr auto a7 = static_cast<T>(-7.364589573262279913270651228486670e-13L);
    
    const T x2 = x * x;
    const T x3 = x2 * x;
    const T x4 = x2 * x2;
    const T x8 = x4 * x4;
    const T x9 = x8 * x;

    const T A = x3 * (a1 + x2 * (a2 + x2 * a3));
    const T B = a4 + x2 * (a5 + x2 * (a6 + x2 * a7));
    const T C = a0 * x;

    return A + C + x9 * B;
}

}

template <typename Real, std::enable_if_t<!std::is_integral_v<Real>, bool> = true>
[[nodiscard]] constexpr Real sin(Real x) noexcept
{
    if (BOOST_MATH_IS_CONSTANT_EVALUATED(x))
    {
        if (boost::math::ccmath::abs(x) == 0)
        {
            return x;
        }
        else if (boost::math::ccmath::isinf(x))
        {
            return x;
        }
        else if (boost::math::ccmath::isnan(x))
        {
            return x;
        }

        return boost::math::ccmath::detail::sin_impl(x);
    }
    else
    {
        using std::sin;
        return sin(x);
    }
}

template <typename Z, std::enable_if_t<std::is_integral_v<Z>, bool> = true>
[[nodiscard]] constexpr double sin(Z x) noexcept
{
    return boost::math::ccmath::sin(static_cast<double>(x));
}

[[nodiscard]] constexpr float sinf(float x) noexcept
{
    return boost::math::ccmath::sin(x);
}

#ifndef BOOST_MATH_NO_LONG_DOUBLE_MATH_FUNCTIONS
[[nodiscard]] constexpr long double sinl(long double x) noexcept
{
    return boost::math::ccmath::sin(x);
}
#endif

}

#endif // BOOST_MATH_CCMATH_TRIG_HPP
