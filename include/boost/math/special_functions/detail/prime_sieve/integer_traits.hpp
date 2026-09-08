//  (C) Copyright Matt Borland 2026.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
//  Traits that let the sieve accept builtin integers and Boost.Multiprecision integers
//  without including Boost.Multiprecision.

#ifndef BOOST_MATH_SF_DETAIL_PRIME_SIEVE_INTEGER_TRAITS_HPP
#define BOOST_MATH_SF_DETAIL_PRIME_SIEVE_INTEGER_TRAITS_HPP

#include <boost/math/tools/config.hpp>
#include <boost/math/tools/traits.hpp>

#ifndef BOOST_MATH_BUILD_MODULE
#include <cstdint>
#include <limits>
#include <type_traits>
#endif

namespace boost::math::detail::prime_sieve {

template <class T>
inline constexpr bool is_multiprecision_v = boost::math::tools::detail::has_backend_type<T>::value;

template <class T>
inline constexpr bool is_integer_like_v = std::is_integral<T>::value || (std::numeric_limits<T>::is_specialized && std::numeric_limits<T>::is_integer);

// Common type of two bounds: the usual promotion for builtins, otherwise the class type.
template <class Lower, class Upper>
using common_integer_t = std::conditional_t<std::is_integral<Lower>::value && std::is_integral<Upper>::value,
                                            std::common_type_t<Lower, Upper>,
                                            std::conditional_t<std::is_integral<Lower>::value, Upper, Lower>>;

// True when the non-negative value x fits in 64 bits.
template <class Integer>
inline bool fits_u64(const Integer& x)
{
    if constexpr (std::is_integral<Integer>::value)
    {
        (void)x;
        return true;
    }
    else
    {
        return x <= Integer(std::numeric_limits<std::uint64_t>::max());
    }
}

template <class Integer>
inline std::uint64_t to_u64(const Integer& x)
{
    return static_cast<std::uint64_t>(x);
}

template <class Integer>
inline Integer from_u64(std::uint64_t x)
{
    return static_cast<Integer>(x);
}

// Clamps negative inputs of signed types to zero.
template <class Integer>
inline Integer clamp_non_negative(const Integer& x)
{
    if constexpr (std::is_integral<Integer>::value && std::is_unsigned<Integer>::value)
    {
        return x;
    }
    else
    {
        return x < Integer(0) ? Integer(0) : x;
    }
}

} // namespace boost::math::detail::prime_sieve

#endif // BOOST_MATH_SF_DETAIL_PRIME_SIEVE_INTEGER_TRAITS_HPP
