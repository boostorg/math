//  (C) Copyright Matt Borland 2022.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_CONCEPTS_CONCEPTS_HPP
#define BOOST_MATH_CONCEPTS_CONCEPTS_HPP

#if (__cplusplus >= 202002L || _MSVC_LANG >= 202002L) && !defined(BOOST_MATH_DISABLE_CONCEPTS)
#if __has_include(<concepts>)

#include <utility>
#include <algorithm>
#include <concepts>
#include <functional>
#include <type_traits>
#include <limits>
#include <boost/math/tools/config.hpp>

namespace boost::math::concepts {

template <typename X, typename Y, typename Op>
struct op_valid_impl
{
    template <typename U, typename L, typename R>
    static constexpr auto test(int) -> decltype(std::declval<U>()(std::declval<L>(), std::declval<R>()),
                                                void(), std::true_type());

    template <typename U, typename L, typename R>
    static constexpr auto test(...) -> std::false_type;

    using type = decltype(test<Op, X, Y>(0));
};

template <typename X, typename Y, typename Op> 
using op_valid = typename op_valid_impl<X, Y, Op>::type;

template <typename X, typename Y, typename Op>
inline constexpr bool op_valid_v = op_valid<X, Y, Op>::value;

template <typename T>
concept Integral = std::is_integral_v<T>
                   #ifdef __SIZEOF_INT128__
                   || std::is_same_v<__int128_t, T>
                   || std::is_same_v<__uint128_t, T>
                   #endif
                   ;

template <typename T>
concept Signed_integral = std::is_integral_v<T> && (std::is_signed_v<T>
                          #ifdef __SIZEOF_INT128__
                          || std::is_same_v<__int128_t, T>
                          #endif
                          );

template <typename T>
concept Unsigned_integral = std::is_integral_v<T> && (std::is_unsigned_v<T>
                            #ifdef __SIZEOF_INT128__
                            || std::is_same_v<__uint128_t, T>
                            #endif
                            );

template <typename T>
concept Real = std::is_floating_point_v<T>
               #ifdef BOOST_HAS_FLOAT128
               #if defined(__INTEL_LLVM_COMPILER) || defined(__INTEL_COMPILER)
               || std::is_same_v<_Quad, T>
               #else
               || std::is_same_v<__float128, T>
               #endif
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

template <typename T>
concept Arbitrary_unsigned_arithmetic_type = Unsigned_arithmetic<T> ||
                                             (op_valid_v<T, T, std::equal_to<>> &&
                                              op_valid_v<T, T, std::not_equal_to<>> &&
                                              op_valid_v<T, T, std::greater<>> &&
                                              op_valid_v<T, T, std::less<>> &&
                                              op_valid_v<T, T, std::greater_equal<>> &&
                                              op_valid_v<T, T, std::less_equal<>> &&
                                              op_valid_v<T, T, std::plus<>> &&
                                              op_valid_v<T, T, std::minus<>> &&
                                              op_valid_v<T, T, std::multiplies<>> &&
                                              op_valid_v<T, T, std::divides<>>);

template <typename T>
concept Arbitrary_signed_arithmetic_type = Signed_arithmetic<T> ||
                                           (Arbitrary_unsigned_arithmetic_type<T> &&
                                            (op_valid_v<T, T, std::negate<>> ||
                                             std::numeric_limits<T>::is_signed));

template <typename T>
concept Arbitrary_arithmetic_type = Arbitrary_unsigned_arithmetic_type<T> ||
                                    Arbitrary_signed_arithmetic_type<T>;

template <typename T>
concept Aribitrary_unsigned_integer_type = Arbitrary_unsigned_arithmetic_type<T> &&
                                           std::numeric_limits<T>::is_integer;

template <typename T>
concept Aribitrary_signed_integer_type = Arbitrary_signed_arithmetic_type<T> &&
                                         std::numeric_limits<T>::is_integer;

template <typename T>
concept Aribitrary_integer_type = Aribitrary_unsigned_integer_type<T> ||
                                  Aribitrary_signed_integer_type<T>;

template <typename T>
concept Aribitrary_real_type = Arbitrary_arithmetic_type<T> &&
                               !std::numeric_limits<T>::is_integer;

}

#define BOOST_MATH_INTEGRAL boost::math::concepts::Integral
#define BOOST_MATH_SIGNED_INTEGRAL boost::math::concepts::Signed_integral
#define BOOST_MATH_UNSIGNED_INTEGRAL boost::math::concepts::Unsigned_integral
#define BOOST_MATH_REAL boost::math::concepts::Real
#define BOOST_MATH_ARITHMETIC boost::math::concepts::Arithmetic
#define BOOST_MATH_SIGNED_ARITHMETIC boost::math::concepts::Signed_arithmetic
#define BOOST_MATH_UNSIGNED_ARITHMETIC boost::math::concepts::Unsigned_arithmetic
#define BOOST_MATH_ARBITRARY_UNSIGNED_ARITHMETIC boost::math::concepts::Arbitrary_unsigned_arithmetic_type
#define BOOST_MATH_ARBITRARY_SIGNED_ARITHMETIC boost::math::concepts::Arbitrary_signed_arithmetic_type
#define BOOST_MATH_ARBITRARY_ARITHMETIC boost::math::concepts::Arbitrary_arithmetic_type
#define BOOST_MATH_ARBITRARY_UNSIGNED_INTEGER boost::math::concepts::Aribitrary_unsigned_integer_type
#define BOOST_MATH_ARBITRARY_SIGNED_INTEGER boost::math::concepts::Aribitrary_signed_integer_type
#define BOOST_MATH_ARBITRARY_INTEGER boost::math::concepts::Aribitrary_integer_type
#define BOOST_MATH_ARBITRARY_REAL boost::math::concepts::Aribitrary_real_type
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

#ifndef BOOST_MATH_ARBITRARY_UNSIGNED_ARITHMETIC
#  define BOOST_MATH_ARBITRARY_UNSIGNED_ARITHMETIC typename
#endif

#ifndef BOOST_MATH_ARBITRARY_SIGNED_ARITHMETIC
#  define BOOST_MATH_ARBITRARY_SIGNED_ARITHMETIC typename
#endif

#ifndef BOOST_MATH_ARBITRARY_ARITHMETIC
#  define BOOST_MATH_ARBITRARY_ARITHMETIC typename
#endif

#ifndef BOOST_MATH_ARBITRARY_UNSIGNED_INTEGER
#  define BOOST_MATH_ARBITRARY_UNSIGNED_INTEGER typename
#endif

#ifndef BOOST_MATH_ARBITRARY_SIGNED_INTEGER
#  define BOOST_MATH_ARBITRARY_SIGNED_INTEGER typename
#endif

#ifndef BOOST_MATH_ARBITRARY_INTEGER
#  define BOOST_MATH_ARBITRARY_INTEGER typename
#endif

#ifndef BOOST_MATH_ARBITRARY_REAL
#  define BOOST_MATH_ARBITRARY_REAL typename
#endif

#ifndef BOOST_MATH_REQUIRES
#  define BOOST_MATH_REQUIRES(X, T)
#endif

#endif // BOOST_MATH_CONCEPTS_CONCEPTS_HPP
