//  (C) Copyright John Maddock 2008 - 2022.
//  (C) Copyright Matt Borland 2022.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_CCMATH_NEXT_HPP
#define BOOST_MATH_CCMATH_NEXT_HPP

#include <cmath>
#include <cfloat>
#include <cstdint>
#include <limits>
#include <type_traits>
#include <boost/math/tools/precision.hpp>
#include <boost/math/tools/traits.hpp>
#include <boost/math/ccmath/ilogb.hpp>
#include <boost/math/ccmath/ldexp.hpp>
#include <boost/math/ccmath/scalbln.hpp>
#include <boost/math/ccmath/round.hpp>
#include <boost/math/ccmath/fabs.hpp>
#include <boost/math/ccmath/fpclassify.hpp>

namespace boost::math::ccmath {

// Forward declarations
template <typename T>
constexpr auto float_next(const T& val);

namespace detail {

template <typename T>
constexpr auto normalize_value(const T& val)
{
    if constexpr (std::numeric_limits<T>::is_specialized && std::numeric_limits<T>::radix != 2)
    {
        std::intmax_t shift = static_cast<std::intmax_t>(std::numeric_limits<T>::digits) - 
                              static_cast<std::intmax_t>(boost::math::ccmath::ilogb(val) - 1);

        T result = boost::math::ccmath::scalbn(val, shift);
        result = boost::math::ccmath::round(result);
        return boost::math::ccmath::scalbn(result, -shift);
    }
    else
    {
        return val;
    }
}

template <typename T>
constexpr T get_smallest_value()
{
    if constexpr (std::numeric_limits<T>::is_specialized && (std::numeric_limits<T>::has_denorm == std::denorm_present))
    {
        if (tools::min_value<T>() / 2 == 0)
        {
            return tools::min_value<T>();
        } 
        else 
        {
            return std::numeric_limits<T>::denorm_min();
        }
    }
    else
    {
        return tools::min_value<T>();
    }
}

template <typename T>
constexpr T get_min_shifted_value()
{
    if constexpr (std::numeric_limits<T>::is_specialized && std::numeric_limits<T>::radix != 2)
    {
        return boost::math::ccmath::scalbln(tools::min_value<T>(), std::numeric_limits<T>::digits + 1);
    }
    else
    {
        return boost::math::ccmath::ldexp(tools::min_value<T>(), tools::digits<T>() + 1);
    }
}

template <typename T>
struct exponent_type
{
    using type = std::conditional_t<boost::math::tools::detail::has_backend_type<T>::value, 
                                    typename T::backend_type::exponent_type, int>;
};

template <typename T>
constexpr T float_next_impl(const T& val)
{
    using exponent_type = typename exponent_type<T>::type;

    exponent_type expon {};

    const int fpclass = (boost::math::ccmath::fpclassify)(val);

    if (fpclass == FP_NAN || fpclass == FP_INFINITE)
    {
        if (val < 0)
        {
            return -tools::max_value<T>();
        }
        else
        {
            static_assert(sizeof(T) == 0, "Argument to float next must be finite");
        }
    }

    if (val >= tools::max_value<T>())
    {
        static_assert(sizeof(T) == 0, "Overflow in constexpr float_next");
    }

    if (val == 0)
    {
        return get_smallest_value<T>();
    }

    T diff {};
    if constexpr (std::numeric_limits<T>::is_specialized && std::numeric_limits<T>::radix != 2)
    {
        if (fpclass != FP_SUBNORMAL && (fpclass != FP_ZERO) && (boost::math::ccmath::fabs(val) < detail::get_min_shifted_value<T>()) && (val != -tools::min_value<T>()))
        {
            // Special case: if the value of the least significant bit is a denorm, and the result
            // would not be a denorm, then shift the input, increment, and shift back.
            // This avoids issues with the Intel SSE2 registers when the FTZ or DAZ flags are set.
            //
            return boost::math::ccmath::scalbn(boost::math::ccmath::float_next(T(boost::math::ccmath::scalbn(val, 2 * std::numeric_limits<T>::digits))), -2 * std::numeric_limits<T>::digits);
        }

        expon = 1 + boost::math::ccmath::ilogb(val);

        if (boost::math::ccmath::scalbn(val, -expon) * std::numeric_limits<T>::radix)
        {
            // reduce expon when val is a power of base and negative
            --expon;
        }
        
        diff = boost::math::ccmath::scalbn(T(1), expon - std::numeric_limits<T>::digits);
    }
    else
    {
        if (fpclass != FP_SUBNORMAL && (fpclass != FP_ZERO) && (boost::math::ccmath::fabs(val) < detail::get_min_shifted_value<T>()) && (val != -tools::min_value<T>()))
        {
            // Special case: if the value of the least significant bit is a denorm, and the result
            // would not be a denorm, then shift the input, increment, and shift back.
            // This avoids issues with the Intel SSE2 registers when the FTZ or DAZ flags are set.
            //
            return boost::math::ccmath::ldexp(boost::math::ccmath::float_next(T(boost::math::ccmath::ldexp(val, 2 * tools::digits<T>()))), -2 * tools::digits<T>());
        }

        if (boost::math::ccmath::frexp(val, &expon) == -0.5)
        {
            // reduce exponent when val is a power of two, and negative.
            --expon;
        }

        diff = boost::math::ccmath::ldexp(T(1), expon - tools::digits<T>());
    }

    if (diff == 0)
    {
        diff = get_smallest_value<T>();
    }

    return val + diff;
}

} // Namespace detail

template <typename T>
constexpr auto float_next(const T& val)
{
    return detail::float_next_impl(detail::normalize_value(val));
}

} // Namespace boost::math::ccmath

#endif // BOOST_MATH_CCMATH_NEXT_HPP
