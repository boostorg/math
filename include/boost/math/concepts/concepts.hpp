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
#include <iterator>
#include <complex>
#include <boost/math/tools/config.hpp>
#include <boost/math/policies/policy.hpp>

namespace boost::math::concepts {

namespace detail {

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
using op_valid_t = typename op_valid_impl<X, Y, Op>::type;

template <typename X, typename Y, typename Op>
inline constexpr bool op_valid_v = op_valid_t<X, Y, Op>::value;

// Detector for class member functions
struct nonesuch 
{
    nonesuch(const nonesuch&) = delete;
    ~nonesuch() = delete;
    void operator=(const nonesuch&) = delete;
};

template <typename Default, typename AlwaysVoid, template<typename...> typename Op, typename... Args>
struct detector 
{
    using value_t = std::false_type;
    using type = Default;
};
 
template <typename Default, template<typename...> typename Op, typename... Args>
struct detector<Default, std::void_t<Op<Args...>>, Op, Args...> 
{
    using value_t = std::true_type;
    using type = Op<Args...>;
};

template <template<typename...> typename Op, typename... Args>
using is_detected = typename detector<nonesuch, void, Op, Args...>::value_t;
 
template <template<typename...> typename Op, typename... Args>
using detected_t = typename detector<nonesuch, void, Op, Args...>::type;

#define BOOST_MATH_HAS_MEMBER_FUNCTION(member)                                      \
template <typename T>                                                               \
using has_##member##_t = decltype(std::declval<T&>().member());                     \
                                                                                    \
template <typename T>                                                               \
inline constexpr bool has_##member##_v = is_detected<has_##member##_t, T>::value;       

BOOST_MATH_HAS_MEMBER_FUNCTION(begin)
BOOST_MATH_HAS_MEMBER_FUNCTION(end)
BOOST_MATH_HAS_MEMBER_FUNCTION(real)
BOOST_MATH_HAS_MEMBER_FUNCTION(imag)

} // Namespace detail

template <typename T>
concept integral = std::is_integral_v<T>;

template <typename T>
concept signed_integral = integral<T> && std::is_signed_v<T>;

template <typename T>
concept unsigned_integral = integral<T> && std::is_unsigned_v<T>;

template <typename T>
concept real = std::is_floating_point_v<T>;

template <typename T>
concept complex = std::is_same_v<T, std::complex<float>>
                  || std::is_same_v<T, std::complex<double>>
                  #ifndef BOOST_MATH_NO_LONG_DOUBLE_MATH_FUNCTIONS
                  || std::is_same_v<T, std::complex<long double>>
                  #endif
                  ;

template <typename T>
concept real_or_complex = real<T> || complex<T>;

template <typename T>
concept arithmetic = std::is_arithmetic_v<T>;

template <typename T>
concept numerical = arithmetic<T> || real_or_complex<T>;

template <typename T>
concept signed_arithmetic = arithmetic<T> && std::is_signed_v<T>;

template <typename T>
concept unsigned_arithmetic = arithmetic<T> && std::is_unsigned_v<T>;

template <typename T>
concept arbitrary_unsigned_arithmetic_type = unsigned_arithmetic<T> ||
                                             (detail::op_valid_v<T, T, std::equal_to<>> &&
                                              detail::op_valid_v<T, T, std::not_equal_to<>> &&
                                              detail::op_valid_v<T, T, std::greater<>> &&
                                              detail::op_valid_v<T, T, std::less<>> &&
                                              detail::op_valid_v<T, T, std::greater_equal<>> &&
                                              detail::op_valid_v<T, T, std::less_equal<>> &&
                                              detail::op_valid_v<T, T, std::plus<>> &&
                                              detail::op_valid_v<T, T, std::minus<>> &&
                                              detail::op_valid_v<T, T, std::multiplies<>> &&
                                              detail::op_valid_v<T, T, std::divides<>>);

template <typename T>
concept arbitrary_signed_arithmetic_type = signed_arithmetic<T> ||
                                           (arbitrary_unsigned_arithmetic_type<T> &&
                                            (detail::op_valid_v<T, T, std::negate<>> ||
                                             std::numeric_limits<T>::is_signed));

template <typename T>
concept arbitrary_arithmetic_type = arbitrary_unsigned_arithmetic_type<T> ||
                                    arbitrary_signed_arithmetic_type<T>;

template <typename T>
concept arbitrary_unsigned_integer_type = arbitrary_unsigned_arithmetic_type<T> &&
                                           std::numeric_limits<T>::is_integer;

template <typename T>
concept arbitrary_signed_integer_type = arbitrary_signed_arithmetic_type<T> &&
                                         std::numeric_limits<T>::is_integer;

template <typename T>
concept arbitrary_integer_type = arbitrary_unsigned_integer_type<T> ||
                                  arbitrary_signed_integer_type<T>;

template <typename T>
concept arbitrary_real_type = arbitrary_arithmetic_type<T> &&
                               !std::numeric_limits<T>::is_integer;

template <typename T>
concept arbitrary_complex_type = complex<T> ||
                                 (detail::has_real_v<T> &&
                                  detail::has_imag_v<T>);

template <typename T>
concept arbitrary_real_or_complex_type = arbitrary_real_type<T> ||
                                         arbitrary_complex_type<T>;

template <typename T>
concept arbitrary_numerical_type = arbitrary_real_or_complex_type<T> ||
                                   arbitrary_arithmetic_type<T>;

template <typename T>
concept policy = boost::math::policies::is_policy<T>::value;

template <typename T>
concept forward_iterator = std::derived_from<typename std::iterator_traits<T>::iterator_category, std::forward_iterator_tag>;

template <typename T>
concept bidirectional_iterator = std::derived_from<typename std::iterator_traits<T>::iterator_category, std::bidirectional_iterator_tag>;

template <typename T>
concept random_access_iterator = std::derived_from<typename std::iterator_traits<T>::iterator_category, std::random_access_iterator_tag>;

template <typename I, typename T>
concept output_iterator = std::derived_from<typename std::iterator_traits<I>::iterator_category, std::input_iterator_tag> &&
                          std::derived_from<typename std::iterator_traits<T>::iterator_category, std::output_iterator_tag>;

template <typename T>
concept is_container = detail::has_begin_v<T> &&
                       detail::has_end_v<T>;

template <typename T>
concept random_access_container = is_container<T> &&
                                  boost::math::concepts::random_access_iterator<typename T::iterator>;

} // boost::math::concepts

#define BOOST_MATH_INTEGRAL boost::math::concepts::integral
#define BOOST_MATH_SIGNED_INTEGRAL boost::math::concepts::signed_integral
#define BOOST_MATH_UNSIGNED_INTEGRAL boost::math::concepts::unsigned_integral
#define BOOST_MATH_REAL boost::math::concepts::real
#define BOOST_MATH_COMPLEX boost::math::concepts::complex
#define BOOST_MATH_REAL_OR_COMPLEX boost::math::concepts::real_or_complex
#define BOOST_MATH_ARITHMETIC boost::math::concepts::arithmetic
#define BOOST_MATH_NUMERICAL boost::math::concepts::numerical
#define BOOST_MATH_SIGNED_ARITHMETIC boost::math::concepts::signed_arithmetic
#define BOOST_MATH_UNSIGNED_ARITHMETIC boost::math::concepts::unsigned_arithmetic
#define BOOST_MATH_ARBITRARY_UNSIGNED_ARITHMETIC boost::math::concepts::arbitrary_unsigned_arithmetic_type
#define BOOST_MATH_ARBITRARY_SIGNED_ARITHMETIC boost::math::concepts::arbitrary_signed_arithmetic_type
#define BOOST_MATH_ARBITRARY_ARITHMETIC boost::math::concepts::arbitrary_arithmetic_type
#define BOOST_MATH_ARBITRARY_UNSIGNED_INTEGER boost::math::concepts::arbitrary_unsigned_integer_type
#define BOOST_MATH_ARBITRARY_SIGNED_INTEGER boost::math::concepts::arbitrary_signed_integer_type
#define BOOST_MATH_ARBITRARY_INTEGER boost::math::concepts::arbitrary_integer_type
#define BOOST_MATH_ARBITRARY_REAL boost::math::concepts::arbitrary_real_type
#define BOOST_MATH_ARBITRARY_COMPLEX boost::math::concepts::arbitrary_complex_type
#define BOOST_MATH_ARBITRARY_REAL_OR_COMPLEX boost::math::concepts::arbitrary_real_or_complex_type
#define BOOST_MATH_ARBITRARY_NUMERICAL boost::math::concepts::arbitrary_numerical_type

#define BOOST_MATH_POLICY boost::math::concepts::policy

#define BOOST_MATH_CONTAINER boost::math::concepts::is_container
#define BOOST_MATH_RANDOM_ACCESS_CONTAINER boost::math::concepts::random_access_container

#define BOOST_MATH_FORWARD_ITER boost::math::concepts::forward_iterator
#define BOOST_MATH_BIDIRECTIONAL_ITER boost::math::concepts::bidirectional_iterator
#define BOOST_MATH_RANDOM_ACCESS_ITER boost::math::concepts::random_access_iterator
#define BOOST_MATH_OUTPUT_ITER(I, T) boost::math::concepts::output_iterator<I, T>
#define BOOST_MATH_REQUIRES_ITER(X) requires X

#define BOOST_MATH_REQUIRES(X, T) requires X<T>

#ifdef BOOST_MATH_EXEC_COMPATIBLE
#include <execution>

namespace boost::math::concepts {

template <typename T>
concept execution_policy = std::is_execution_policy_v<std::remove_cvref_t<T>>;

} // Namespace boost::math::concepts

#define BOOST_MATH_EXECUTION_POLICY boost::math::concepts::execution_policy

#endif // Has <execution>

#endif // Has <concepts>
#endif // C++20

// If concepts are unavailable replace them with typename for compatibility

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

#ifndef BOOST_MATH_COMPLEX
#  define BOOST_MATH_COMPLEX typename
#endif

#ifndef BOOST_MATH_REAL_OR_COMPLEX
#  define BOOST_MATH_REAL_OR_COMPLEX typename
#endif

#ifndef BOOST_MATH_ARITHMETIC
#  define BOOST_MATH_ARITHMETIC typename
#endif

#ifndef BOOST_MATH_NUMERICAL
#  define BOOST_MATH_NUMERICAL typename
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

#ifndef BOOST_MATH_ARBITRARY_COMPLEX
#  define BOOST_MATH_ARBITRARY_COMPLEX typename
#endif

#ifndef BOOST_MATH_ARBITRARY_REAL_OR_COMPLEX
#  define BOOST_MATH_ARBITRARY_REAL_OR_COMPLEX typename
#endif

#ifndef BOOST_MATH_ARBITRARY_NUMERICAL
#  define BOOST_MATH_ARBITRARY_NUMERICAL typename
#endif

#ifndef BOOST_MATH_POLICY
#  define BOOST_MATH_POLICY typename
#endif

#ifndef BOOST_MATH_FORWARD_ITER
#  define BOOST_MATH_FORWARD_ITER typename
#endif

#ifndef BOOST_MATH_BIDIRECTIONAL_ITER
#  define BOOST_MATH_BIDIRECTIONAL_ITER typename
#endif

#ifndef BOOST_MATH_RANDOM_ACCESS_ITER
#  define BOOST_MATH_RANDOM_ACCESS_ITER typename
#endif

#ifndef BOOST_MATH_OUTPUT_ITER
#  define BOOST_MATH_OUTPUT_ITER(I, T)
#endif

#ifndef BOOST_MATH_REQUIRES_ITER
#  define BOOST_MATH_REQUIRES_ITER(X)
#endif

#ifndef BOOST_MATH_CONTAINER
#  define BOOST_MATH_CONTAINER typename
#endif

#ifndef BOOST_MATH_RANDOM_ACCESS_CONTAINER
#  define BOOST_MATH_RANDOM_ACCESS_CONTAINER typename
#endif

#ifndef BOOST_MATH_EXECUTION_POLICY
#  define BOOST_MATH_EXECUTION_POLICY typename
#endif

#ifndef BOOST_MATH_REQUIRES
#  define BOOST_MATH_REQUIRES(X, T)
#endif

#endif // BOOST_MATH_CONCEPTS_CONCEPTS_HPP
