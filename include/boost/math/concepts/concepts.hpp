//  (C) Copyright Matt Borland 2022.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_CONCEPTS_CONCEPTS_HPP
#define BOOST_MATH_CONCEPTS_CONCEPTS_HPP

#if __cplusplus >= 202002L || _MSVC_LANG >= 202002L
#if __has_include(<concepts>)
#include <concepts>

#include <type_traits>
#include <boost/math/tools/config.hpp>

namespace boost::math::concepts {

template <typename T>
concept Integral = std::is_integral_v<T>
                   #ifdef __SIZEOF_INT128__
                   || std::is_same_v<__int128_t, T>
                   || std::is_same_v<__uint128_t, T>
                   #endif
                   ;


template <typename T>
concept Signed_integral = std::is_integral_v<T> && std::is_signed_v<T>
                          #ifdef __SIZEOF_INT128__
                          || std::is_same_v<__int128_t, T>
                          #endif
                          ;

template <typename T>
concept Unsigned_integral = std::is_integral_v<T> && std::is_unsigned_v<T>
                            #ifdef __SIZEOF_INT128__
                            || std::is_same_v<__uint128_t, T>
                            #endif
                            ;

template <typename T>
concept Real = std::is_floating_point_v<T>
               #ifdef BOOST_HAS_FLOAT128
               || std::is_same_v<__float128, T>
               #endif
               ;

template <typename T>
concept Arithmetic = Integral<T> || Real<T>;

template <typename T>
concept Signed_arithmetic = Arithmetic<T> && (std::is_signed_v<T>
                            #ifdef __SIZEOF_INT128__
                            || std::is_same_v<__int128_t, T>
                            #endif
                            );

template <typename T>
concept Unsigned_arithmetic = Arithmetic<T> && (std::is_unsigned_v<T>
                            #ifdef __SIZEOF_INT128__
                            || std::is_same_v<__uint128_t, T>
                            #endif
                            );

}

#define BOOST_MATH_INTEGRAL boost::math::concepts::Integral<T>
#define BOOST_MATH_SIGNED_INTEGRAL boost::math::concepts::Signed_integral<T>
#define BOOST_MATH_UNSIGNED_INTEGRAL boost::math::concepts::Unsigned_integral<T>
#define BOOST_MATH_REAL boost::math::concepts::Real
#define BOOST_MATH_ARITHMETIC boost::math::concepts::Arithmetic
#define BOOST_MATH_SIGNED_ARITHMETIC boost::math::concepts::Signed_arithmetic
#define BOOST_MATH_UNSIGNED_ARITHMETIC boost::math::concepts::Unsigned_arithmetic
#define BOOST_MATH_REQUIRES(X, T) requires X<T>

#endif
#endif

#ifndef BOOST_MATH_INTEGRAL
#  define BOOST_MATH_INTEGRAL typename
#endif

#ifndef BOOST_MATH_SIGNED_INTEGRAL
#  define BOOST_MATH_SIGNED_INTEGRAL typename
#endif

#ifndef BOOST_MATH_UNSIGNED_INTEGRAL
#  define BOOST_MATH_UNSIGNED_INTEGRAL typename
#endif

#ifndef BOOST_MATH_REAL
#  define BOOST_MATH_REAL typename
#endif

#ifndef BOOST_MATH_ARITHMETIC
#  define BOOST_MATH_ARITHMETIC typename
#endif

#ifndef BOOST_MATH_SIGNED_ARITHMETIC
#  define BOOST_MATH_SIGNED_ARITHMETIC typename
#endif

#ifndef BOOST_MATH_UNSIGNED_ARITHMETIC
#  define BOOST_MATH_UNSIGNED_ARITHMETIC typename
#endif

#ifndef BOOST_MATH_REQUIRES
#  define BOOST_MATH_REQUIRES(X, T)
#endif

#endif // BOOST_MATH_CONCEPTS_CONCEPTS_HPP
