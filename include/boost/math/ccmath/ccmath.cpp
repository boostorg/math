//  (C) Copyright Matt Borland 2021.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
//  A module containing all of the constexpr implementation of <cmath>

// Global module fragment required for non-module preprocessing
module;

#include <boost/math/tools/config.hpp>
#include <boost/math/tools/is_constant_evaluated.hpp>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cinttypes>
#include <type_traits>
#include <limits>
#include <stdexcept>

export module boost_math_ccmath;

// Forward Declarations
export namespace boost::math::ccmath {

template <typename T>
inline constexpr bool isnan(T x);

template <typename T>
inline constexpr bool isinf(T x);

template <typename T, std::enable_if_t<!std::is_unsigned_v<T>, bool> = true>
inline constexpr T abs(T x) noexcept;

template <typename T, std::enable_if_t<std::is_unsigned_v<T>, bool> = true>
inline constexpr T abs(T x) noexcept;

}

namespace boost::math::ccmath::detail {

template <typename T> 
inline constexpr T abs_impl(T x) noexcept
{
    return boost::math::ccmath::isnan(x) ? std::numeric_limits<T>::quiet_NaN() : 
           boost::math::ccmath::isinf(x) ? std::numeric_limits<T>::infinity() : 
           x == -0 ? T(0) :
           x == (std::numeric_limits<T>::min)() ? std::numeric_limits<T>::quiet_NaN() : 
           x > 0 ? x : -x;
}

}

// Useable Functions

export namespace boost::math::ccmath {

template <typename T>
inline constexpr bool isnan(T x)
{
    if(BOOST_MATH_IS_CONSTANT_EVALUATED(x))
    {
        return x != x;
    }
    else
    {
        using std::isnan;

        if constexpr (!std::is_integral_v<T>)
        {
            return isnan(x);
        }
        else
        {
            return isnan(static_cast<double>(x));
        }
    }
}

template <typename T>
inline constexpr bool isinf(T x)
{
    if(BOOST_MATH_IS_CONSTANT_EVALUATED(x))
    {
        return x == std::numeric_limits<T>::infinity() || -x == std::numeric_limits<T>::infinity();
    }
    else
    {
        using std::isinf;
        
        if constexpr (!std::is_integral_v<T>)
        {
            return isinf(x);
        }
        else
        {
            return isinf(static_cast<double>(x));
        }
    }
}

template <typename T, std::enable_if_t<!std::is_unsigned_v<T>, bool> = true>
inline constexpr T abs(T x) noexcept
{
    if(BOOST_MATH_IS_CONSTANT_EVALUATED(x))
    {
        return detail::abs_impl<T>(x);
    }
    else
    {
        using std::abs;
        return abs(x);
    }
}

// If abs() is called with an argument of type X for which is_unsigned_v<X> is true and if X
// cannot be converted to int by integral promotion (7.3.7), the program is ill-formed.
template <typename T, std::enable_if_t<std::is_unsigned_v<T>, bool> = true>
inline constexpr T abs(T x) noexcept
{
    if constexpr (std::is_convertible_v<T, int>)
    {
        return detail::abs_impl<int>(static_cast<int>(x));
    }
    else
    {
        static_assert(sizeof(T) == 0, "Taking the absolute value of an unsigned value not covertible to int is UB.");
        return T(0); // Unreachable, but suppresses warnings
    }
}

inline constexpr long int labs(long int j) noexcept
{
    return boost::math::ccmath::abs(j);
}

inline constexpr long long int llabs(long long int j) noexcept
{
    return boost::math::ccmath::abs(j);
}

}
