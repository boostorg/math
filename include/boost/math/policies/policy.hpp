//  Copyright John Maddock 2007.
//  Copyright Matt Borland 2021.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_POLICY_HPP
#define BOOST_MATH_POLICY_HPP

#include <boost/mp11.hpp>
#include <boost/mp11/mpl.hpp>
#include <boost/mpl/contains.hpp>
#include <boost/mpl/find_if.hpp>
#include <boost/mpl/comparison.hpp>
#include <boost/math/tools/config.hpp>
#include <boost/math/tools/mp.hpp>
#include <limits>
#include <type_traits>
#include <cmath>
#include <cstdint>

using namespace boost::mp11;

namespace boost{ namespace math{ 

namespace tools{

template <class T>
constexpr int digits(BOOST_MATH_EXPLICIT_TEMPLATE_TYPE(T)) noexcept;
template <class T>
constexpr T epsilon(BOOST_MATH_EXPLICIT_TEMPLATE_TYPE(T)) noexcept(std::is_floating_point<T>::value);

}

namespace policies{

//
// Define macros for our default policies, if they're not defined already:
//
// Special cases for exceptions disabled first:
//
#ifdef BOOST_NO_EXCEPTIONS
#  ifndef BOOST_MATH_DOMAIN_ERROR_POLICY
#    define BOOST_MATH_DOMAIN_ERROR_POLICY errno_on_error
#  endif
#  ifndef BOOST_MATH_POLE_ERROR_POLICY
#     define BOOST_MATH_POLE_ERROR_POLICY errno_on_error
#  endif
#  ifndef BOOST_MATH_OVERFLOW_ERROR_POLICY
#     define BOOST_MATH_OVERFLOW_ERROR_POLICY errno_on_error
#  endif
#  ifndef BOOST_MATH_EVALUATION_ERROR_POLICY
#     define BOOST_MATH_EVALUATION_ERROR_POLICY errno_on_error
#  endif
#  ifndef BOOST_MATH_ROUNDING_ERROR_POLICY
#     define BOOST_MATH_ROUNDING_ERROR_POLICY errno_on_error
#  endif
#endif
//
// Then the regular cases:
//
#ifndef BOOST_MATH_DOMAIN_ERROR_POLICY
#define BOOST_MATH_DOMAIN_ERROR_POLICY throw_on_error
#endif
#ifndef BOOST_MATH_POLE_ERROR_POLICY
#define BOOST_MATH_POLE_ERROR_POLICY throw_on_error
#endif
#ifndef BOOST_MATH_OVERFLOW_ERROR_POLICY
#define BOOST_MATH_OVERFLOW_ERROR_POLICY throw_on_error
#endif
#ifndef BOOST_MATH_EVALUATION_ERROR_POLICY
#define BOOST_MATH_EVALUATION_ERROR_POLICY throw_on_error
#endif
#ifndef BOOST_MATH_ROUNDING_ERROR_POLICY
#define BOOST_MATH_ROUNDING_ERROR_POLICY throw_on_error
#endif
#ifndef BOOST_MATH_UNDERFLOW_ERROR_POLICY
#define BOOST_MATH_UNDERFLOW_ERROR_POLICY ignore_error
#endif
#ifndef BOOST_MATH_DENORM_ERROR_POLICY
#define BOOST_MATH_DENORM_ERROR_POLICY ignore_error
#endif
#ifndef BOOST_MATH_INDETERMINATE_RESULT_ERROR_POLICY
#define BOOST_MATH_INDETERMINATE_RESULT_ERROR_POLICY ignore_error
#endif
#ifndef BOOST_MATH_DIGITS10_POLICY
#define BOOST_MATH_DIGITS10_POLICY 0
#endif
#ifndef BOOST_MATH_PROMOTE_FLOAT_POLICY
#define BOOST_MATH_PROMOTE_FLOAT_POLICY true
#endif
#ifndef BOOST_MATH_PROMOTE_DOUBLE_POLICY
#ifdef BOOST_MATH_NO_LONG_DOUBLE_MATH_FUNCTIONS
#define BOOST_MATH_PROMOTE_DOUBLE_POLICY false
#else
#define BOOST_MATH_PROMOTE_DOUBLE_POLICY true
#endif
#endif
#ifndef BOOST_MATH_DISCRETE_QUANTILE_POLICY
#define BOOST_MATH_DISCRETE_QUANTILE_POLICY integer_round_outwards
#endif
#ifndef BOOST_MATH_ASSERT_UNDEFINED_POLICY
#define BOOST_MATH_ASSERT_UNDEFINED_POLICY true
#endif
#ifndef BOOST_MATH_MAX_SERIES_ITERATION_POLICY
#define BOOST_MATH_MAX_SERIES_ITERATION_POLICY 1000000
#endif
#ifndef BOOST_MATH_MAX_ROOT_ITERATION_POLICY
#define BOOST_MATH_MAX_ROOT_ITERATION_POLICY 200
#endif

#define BOOST_MATH_META_INT(Type, name, Default)                                               \
   template <Type N = Default>                                                                 \
   class name : public std::integral_constant<int, N>{};                                        \
                                                                                                \
   namespace detail{                                                                            \
   template <Type N>                                                                           \
   char test_is_valid_arg(const name<N>* = nullptr);                                            \
   char test_is_default_arg(const name<Default>* = nullptr);                                    \
                                                                                                \
   template <typename T>                                                                        \
   class is_##name##_imp                                                                        \
   {                                                                                            \
   private:                                                                                     \
      template <Type N>                                                                        \
      static char test(const name<N>* = nullptr);                                               \
      static double test(...);                                                                  \
   public:                                                                                      \
      static constexpr bool value = sizeof(test(static_cast<T*>(0))) == sizeof(char);           \
   };                                                                                           \
   }                                                                                            \
                                                                                                \
   template <typename T>                                                                        \
   class is_##name                                                                              \
   {                                                                                            \
   public:                                                                                      \
      static constexpr bool value = boost::math::policies::detail::is_##name##_imp<T>::value;   \
      using type = std::integral_constant<bool, value>;                                         \
   };

#define BOOST_MATH_META_BOOL(name, Default)                                                     \
   template <bool N = Default>                                                                  \
   class name : public std::integral_constant<bool, N>{};                                       \
                                                                                                \
   namespace detail{                                                                            \
   template <bool N>                                                                            \
   char test_is_valid_arg(const name<N>* = nullptr);                                            \
   char test_is_default_arg(const name<Default>* = nullptr);                                    \
                                                                                                \
   template <typename T>                                                                        \
   class is_##name##_imp                                                                        \
   {                                                                                            \
   private:                                                                                     \
      template <bool N>                                                                         \
      static char test(const name<N>* = nullptr);                                               \
      static double test(...);                                                                  \
   public:                                                                                      \
      static constexpr bool value = sizeof(test(static_cast<T*>(0))) == sizeof(char);           \
   };                                                                                           \
   }                                                                                            \
                                                                                                \
   template <typename T>                                                                        \
   class is_##name                                                                              \
   {                                                                                            \
   public:                                                                                      \
      static constexpr bool value = boost::math::policies::detail::is_##name##_imp<T>::value;   \
      using type = std::integral_constant<bool, value>;                                         \
   };

//
// Begin by defining policy types for error handling:
//
enum error_policy_type
{
   throw_on_error = 0,
   errno_on_error = 1,
   ignore_error = 2,
   user_error = 3
};

BOOST_MATH_META_INT(error_policy_type, domain_error, BOOST_MATH_DOMAIN_ERROR_POLICY)
BOOST_MATH_META_INT(error_policy_type, pole_error, BOOST_MATH_POLE_ERROR_POLICY)
BOOST_MATH_META_INT(error_policy_type, overflow_error, BOOST_MATH_OVERFLOW_ERROR_POLICY)
BOOST_MATH_META_INT(error_policy_type, underflow_error, BOOST_MATH_UNDERFLOW_ERROR_POLICY)
BOOST_MATH_META_INT(error_policy_type, denorm_error, BOOST_MATH_DENORM_ERROR_POLICY)
BOOST_MATH_META_INT(error_policy_type, evaluation_error, BOOST_MATH_EVALUATION_ERROR_POLICY)
BOOST_MATH_META_INT(error_policy_type, rounding_error, BOOST_MATH_ROUNDING_ERROR_POLICY)
BOOST_MATH_META_INT(error_policy_type, indeterminate_result_error, BOOST_MATH_INDETERMINATE_RESULT_ERROR_POLICY)

//
// Policy types for internal promotion:
//
BOOST_MATH_META_BOOL(promote_float, BOOST_MATH_PROMOTE_FLOAT_POLICY)
BOOST_MATH_META_BOOL(promote_double, BOOST_MATH_PROMOTE_DOUBLE_POLICY)
BOOST_MATH_META_BOOL(assert_undefined, BOOST_MATH_ASSERT_UNDEFINED_POLICY)
//
// Policy types for discrete quantiles:
//
enum discrete_quantile_policy_type
{
   real,
   integer_round_outwards,
   integer_round_inwards,
   integer_round_down,
   integer_round_up,
   integer_round_nearest
};

BOOST_MATH_META_INT(discrete_quantile_policy_type, discrete_quantile, BOOST_MATH_DISCRETE_QUANTILE_POLICY)
//
// Precision:
//
BOOST_MATH_META_INT(int, digits10, BOOST_MATH_DIGITS10_POLICY)
BOOST_MATH_META_INT(int, digits2, 0)
//
// Iterations:
//
BOOST_MATH_META_INT(unsigned long, max_series_iterations, BOOST_MATH_MAX_SERIES_ITERATION_POLICY)
BOOST_MATH_META_INT(unsigned long, max_root_iterations, BOOST_MATH_MAX_ROOT_ITERATION_POLICY)
//
// Define the names for each possible policy:
//
#define BOOST_MATH_PARAMETER(name)\
   BOOST_PARAMETER_TEMPLATE_KEYWORD(name##_name)\
   BOOST_PARAMETER_NAME(name##_name)

struct default_policy{};

namespace detail{
//
// Trait to work out bits precision from digits10 and digits2:
//
template <class Digits10, class Digits2>
struct precision
{
   //
   // Now work out the precision:
   //
   using digits2_type = typename std::conditional<
      (Digits10::value == 0),
      digits2<0>,
      digits2<((Digits10::value + 1) * 1000L) / 301L>
   >::type;
public:
#ifdef BOOST_BORLANDC
   using type = typename std::conditional<
      (Digits2::value > ::boost::math::policies::detail::precision<Digits10,Digits2>::digits2_type::value),
      Digits2, digits2_type>::type;
#else
   using type = typename std::conditional<
      (Digits2::value > digits2_type::value),
      Digits2, digits2_type>::type;
#endif
};

template <class A, class B, bool b>
struct select_result
{
   using type = A;
};

// Type that underlying the iterator position (e.g. std::iterator_traits<B>::value_type)
template <class A, class B>
struct select_result<A, B, false>
{
   using type = typename mpl::deref<B>::type;
};

// Select the argument based on if the found iterator is the end
// If true: use the default type
// If false: use the underlying type from the iterator position
template <class Seq, class Pred, class DefaultType>
struct find_arg
{
private:
   using iter = typename mpl::find_if<Seq, Pred>::type;
   using end_type = typename mpl::end<Seq>::type;
public:
   using type = typename select_result<DefaultType, iter, std::is_same<iter, end_type>::value>::type;
};

double test_is_valid_arg(...);
double test_is_default_arg(...);
char test_is_valid_arg(const default_policy*);
char test_is_default_arg(const default_policy*);

template <typename T>
class is_valid_policy_imp 
{
public:
   static constexpr bool value = sizeof(boost::math::policies::detail::test_is_valid_arg(static_cast<T*>(0))) == sizeof(char);
};

template <typename T> 
class is_valid_policy 
{
public:
   static constexpr bool value = boost::math::policies::detail::is_valid_policy_imp<T>::value;
};

template <typename T>
class is_default_policy_imp
{
public:
   static constexpr bool value = sizeof(boost::math::policies::detail::test_is_default_arg(static_cast<T*>(0))) == sizeof(char);
};

template <typename T>
class is_default_policy
{
public:
   static constexpr bool value = boost::math::policies::detail::is_default_policy_imp<T>::value;
   using type = std::integral_constant<bool, value>;

   template <typename U>
   struct apply
   {
      using type = is_default_policy<U>;
   };
};

template <class Seq, class T, std::size_t N>
struct append_N
{
   using type = typename append_N<mp_push_back<Seq, T>, T, N-1>::type;
};

template <class Seq, class T>
struct append_N<Seq, T, 0>
{
   using type = Seq;
};

//
// Traits class to work out what template parameters our default
// policy<> class will have when modified for forwarding:
//
template <bool f, bool d>
struct default_args
{
   typedef promote_float<false> arg1;
   typedef promote_double<false> arg2;
};

template <>
struct default_args<false, false>
{
   typedef default_policy arg1;
   typedef default_policy arg2;
};

template <>
struct default_args<true, false>
{
   typedef promote_float<false> arg1;
   typedef default_policy arg2;
};

template <>
struct default_args<false, true>
{
   typedef promote_double<false> arg1;
   typedef default_policy arg2;
};

typedef default_args<BOOST_MATH_PROMOTE_FLOAT_POLICY, BOOST_MATH_PROMOTE_DOUBLE_POLICY>::arg1 forwarding_arg1;
typedef default_args<BOOST_MATH_PROMOTE_FLOAT_POLICY, BOOST_MATH_PROMOTE_DOUBLE_POLICY>::arg2 forwarding_arg2;

} // detail
//
// Now define the policy type with enough arguments to handle all
// the policies:
//
template <class A1 = default_policy, 
          class A2 = default_policy, 
          class A3 = default_policy,
          class A4 = default_policy,
          class A5 = default_policy,
          class A6 = default_policy,
          class A7 = default_policy,
          class A8 = default_policy,
          class A9 = default_policy,
          class A10 = default_policy,
          class A11 = default_policy,
          class A12 = default_policy,
          class A13 = default_policy>
struct policy
{
private:
   //
   // Validate all our arguments:
   //
   static_assert(::boost::math::policies::detail::is_valid_policy<A1>::value, "::boost::math::policies::detail::is_valid_policy<A1>::value");
   static_assert(::boost::math::policies::detail::is_valid_policy<A2>::value, "::boost::math::policies::detail::is_valid_policy<A2>::value");
   static_assert(::boost::math::policies::detail::is_valid_policy<A3>::value, "::boost::math::policies::detail::is_valid_policy<A3>::value");
   static_assert(::boost::math::policies::detail::is_valid_policy<A4>::value, "::boost::math::policies::detail::is_valid_policy<A4>::value");
   static_assert(::boost::math::policies::detail::is_valid_policy<A5>::value, "::boost::math::policies::detail::is_valid_policy<A5>::value");
   static_assert(::boost::math::policies::detail::is_valid_policy<A6>::value, "::boost::math::policies::detail::is_valid_policy<A6>::value");
   static_assert(::boost::math::policies::detail::is_valid_policy<A7>::value, "::boost::math::policies::detail::is_valid_policy<A7>::value");
   static_assert(::boost::math::policies::detail::is_valid_policy<A8>::value, "::boost::math::policies::detail::is_valid_policy<A8>::value");
   static_assert(::boost::math::policies::detail::is_valid_policy<A9>::value, "::boost::math::policies::detail::is_valid_policy<A9>::value");
   static_assert(::boost::math::policies::detail::is_valid_policy<A10>::value, "::boost::math::policies::detail::is_valid_policy<A10>::value");
   static_assert(::boost::math::policies::detail::is_valid_policy<A11>::value, "::boost::math::policies::detail::is_valid_policy<A11>::value");
   static_assert(::boost::math::policies::detail::is_valid_policy<A12>::value, "::boost::math::policies::detail::is_valid_policy<A12>::value");
   static_assert(::boost::math::policies::detail::is_valid_policy<A13>::value, "::boost::math::policies::detail::is_valid_policy<A13>::value");
   //
   // Typelist of the arguments:
   //
   using arg_list = mp_list<A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11,A12,A13>;

public:
   typedef typename detail::find_arg<arg_list, is_domain_error<mpl::_1>, domain_error<> >::type domain_error_type;
   typedef typename detail::find_arg<arg_list, is_pole_error<mpl::_1>, pole_error<> >::type pole_error_type;
   typedef typename detail::find_arg<arg_list, is_overflow_error<mpl::_1>, overflow_error<> >::type overflow_error_type;
   typedef typename detail::find_arg<arg_list, is_underflow_error<mpl::_1>, underflow_error<> >::type underflow_error_type;
   typedef typename detail::find_arg<arg_list, is_denorm_error<mpl::_1>, denorm_error<> >::type denorm_error_type;
   typedef typename detail::find_arg<arg_list, is_evaluation_error<mpl::_1>, evaluation_error<> >::type evaluation_error_type;
   typedef typename detail::find_arg<arg_list, is_rounding_error<mpl::_1>, rounding_error<> >::type rounding_error_type;
   typedef typename detail::find_arg<arg_list, is_indeterminate_result_error<mpl::_1>, indeterminate_result_error<> >::type indeterminate_result_error_type;
private:
   //
   // Now work out the precision:
   //
   typedef typename detail::find_arg<arg_list, is_digits10<mpl::_1>, digits10<> >::type digits10_type;
   typedef typename detail::find_arg<arg_list, is_digits2<mpl::_1>, digits2<> >::type bits_precision_type;
public:
   typedef typename detail::precision<digits10_type, bits_precision_type>::type precision_type;
   //
   // Internal promotion:
   //
   typedef typename detail::find_arg<arg_list, is_promote_float<mpl::_1>, promote_float<> >::type promote_float_type;
   typedef typename detail::find_arg<arg_list, is_promote_double<mpl::_1>, promote_double<> >::type promote_double_type;
   //
   // Discrete quantiles:
   //
   typedef typename detail::find_arg<arg_list, is_discrete_quantile<mpl::_1>, discrete_quantile<> >::type discrete_quantile_type;
   //
   // Mathematically undefined properties:
   //
   typedef typename detail::find_arg<arg_list, is_assert_undefined<mpl::_1>, assert_undefined<> >::type assert_undefined_type;
   //
   // Max iterations:
   //
   typedef typename detail::find_arg<arg_list, is_max_series_iterations<mpl::_1>, max_series_iterations<> >::type max_series_iterations_type;
   typedef typename detail::find_arg<arg_list, is_max_root_iterations<mpl::_1>, max_root_iterations<> >::type max_root_iterations_type;
};
//
// These full specializations are defined to reduce the amount of
// template instantiations that have to take place when using the default
// policies, they have quite a large impact on compile times:
//
template <>
struct policy<default_policy, default_policy, default_policy, default_policy, default_policy, default_policy, default_policy, default_policy, default_policy, default_policy, default_policy>
{
public:
   typedef domain_error<> domain_error_type;
   typedef pole_error<> pole_error_type;
   typedef overflow_error<> overflow_error_type;
   typedef underflow_error<> underflow_error_type;
   typedef denorm_error<> denorm_error_type;
   typedef evaluation_error<> evaluation_error_type;
   typedef rounding_error<> rounding_error_type;
   typedef indeterminate_result_error<> indeterminate_result_error_type;
#if BOOST_MATH_DIGITS10_POLICY == 0
   typedef digits2<> precision_type;
#else
   typedef detail::precision<digits10<>, digits2<> >::type precision_type;
#endif
   typedef promote_float<> promote_float_type;
   typedef promote_double<> promote_double_type;
   typedef discrete_quantile<> discrete_quantile_type;
   typedef assert_undefined<> assert_undefined_type;
   typedef max_series_iterations<> max_series_iterations_type;
   typedef max_root_iterations<> max_root_iterations_type;
};

template <>
struct policy<detail::forwarding_arg1, detail::forwarding_arg2, default_policy, default_policy, default_policy, default_policy, default_policy, default_policy, default_policy, default_policy, default_policy>
{
public:
   typedef domain_error<> domain_error_type;
   typedef pole_error<> pole_error_type;
   typedef overflow_error<> overflow_error_type;
   typedef underflow_error<> underflow_error_type;
   typedef denorm_error<> denorm_error_type;
   typedef evaluation_error<> evaluation_error_type;
   typedef rounding_error<> rounding_error_type;
   typedef indeterminate_result_error<> indeterminate_result_error_type;
#if BOOST_MATH_DIGITS10_POLICY == 0
   typedef digits2<> precision_type;
#else
   typedef detail::precision<digits10<>, digits2<> >::type precision_type;
#endif
   typedef promote_float<false> promote_float_type;
   typedef promote_double<false> promote_double_type;
   typedef discrete_quantile<> discrete_quantile_type;
   typedef assert_undefined<> assert_undefined_type;
   typedef max_series_iterations<> max_series_iterations_type;
   typedef max_root_iterations<> max_root_iterations_type;
};

template <class Policy, 
          class A1 = default_policy, 
          class A2 = default_policy, 
          class A3 = default_policy,
          class A4 = default_policy,
          class A5 = default_policy,
          class A6 = default_policy,
          class A7 = default_policy,
          class A8 = default_policy,
          class A9 = default_policy,
          class A10 = default_policy,
          class A11 = default_policy,
          class A12 = default_policy,
          class A13 = default_policy>
struct normalise
{
private:
   using arg_list = mp_list<A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11,A12,A13>;
   typedef typename detail::find_arg<arg_list, is_domain_error<mpl::_1>, typename Policy::domain_error_type >::type domain_error_type;
   typedef typename detail::find_arg<arg_list, is_pole_error<mpl::_1>, typename Policy::pole_error_type >::type pole_error_type;
   typedef typename detail::find_arg<arg_list, is_overflow_error<mpl::_1>, typename Policy::overflow_error_type >::type overflow_error_type;
   typedef typename detail::find_arg<arg_list, is_underflow_error<mpl::_1>, typename Policy::underflow_error_type >::type underflow_error_type;
   typedef typename detail::find_arg<arg_list, is_denorm_error<mpl::_1>, typename Policy::denorm_error_type >::type denorm_error_type;
   typedef typename detail::find_arg<arg_list, is_evaluation_error<mpl::_1>, typename Policy::evaluation_error_type >::type evaluation_error_type;
   typedef typename detail::find_arg<arg_list, is_rounding_error<mpl::_1>, typename Policy::rounding_error_type >::type rounding_error_type;
   typedef typename detail::find_arg<arg_list, is_indeterminate_result_error<mpl::_1>, typename Policy::indeterminate_result_error_type >::type indeterminate_result_error_type;
   //
   // Now work out the precision:
   //
   typedef typename detail::find_arg<arg_list, is_digits10<mpl::_1>, digits10<> >::type digits10_type;
   typedef typename detail::find_arg<arg_list, is_digits2<mpl::_1>, typename Policy::precision_type >::type bits_precision_type;
   typedef typename detail::precision<digits10_type, bits_precision_type>::type precision_type;
   //
   // Internal promotion:
   //
   typedef typename detail::find_arg<arg_list, is_promote_float<mpl::_1>, typename Policy::promote_float_type >::type promote_float_type;
   typedef typename detail::find_arg<arg_list, is_promote_double<mpl::_1>, typename Policy::promote_double_type >::type promote_double_type;
   //
   // Discrete quantiles:
   //
   typedef typename detail::find_arg<arg_list, is_discrete_quantile<mpl::_1>, typename Policy::discrete_quantile_type >::type discrete_quantile_type;
   //
   // Mathematically undefined properties:
   //
   typedef typename detail::find_arg<arg_list, is_assert_undefined<mpl::_1>, typename Policy::assert_undefined_type >::type assert_undefined_type;
   //
   // Max iterations:
   //
   typedef typename detail::find_arg<arg_list, is_max_series_iterations<mpl::_1>, typename Policy::max_series_iterations_type>::type max_series_iterations_type;
   typedef typename detail::find_arg<arg_list, is_max_root_iterations<mpl::_1>, typename Policy::max_root_iterations_type>::type max_root_iterations_type;
   //
   // Define a typelist of the policies:
   //
   using result_list = mp_list<
      domain_error_type,
      pole_error_type,
      overflow_error_type,
      underflow_error_type,
      denorm_error_type,
      evaluation_error_type,
      rounding_error_type,
      indeterminate_result_error_type,
      precision_type,
      promote_float_type,
      promote_double_type,
      discrete_quantile_type,
      assert_undefined_type,
      max_series_iterations_type,
      max_root_iterations_type>;
   //
   // Remove all the policies that are the same as the default:
   //
   using fn = mp_quote_trait<detail::is_default_policy>;
   using reduced_list = mp_remove_if_q<result_list, fn>;
   //
   // Pad out the list with defaults:
   //
   using result_type = typename detail::append_N<reduced_list, default_policy, (15UL - mp_size<reduced_list>::value)>::type;
public:
   using type = policy<
      mp_at_c<result_type, 0>,
      mp_at_c<result_type, 1>,
      mp_at_c<result_type, 2>,
      mp_at_c<result_type, 3>,
      mp_at_c<result_type, 4>,
      mp_at_c<result_type, 5>,
      mp_at_c<result_type, 6>,
      mp_at_c<result_type, 7>,
      mp_at_c<result_type, 8>,
      mp_at_c<result_type, 9>,
      mp_at_c<result_type, 10>,
      mp_at_c<result_type, 11>,
      mp_at_c<result_type, 12>
      >;
};
//
// Full specialisation to speed up compilation of the common case:
//
template <>
struct normalise<policy<>, 
          promote_float<false>, 
          promote_double<false>, 
          discrete_quantile<>,
          assert_undefined<>,
          default_policy,
          default_policy,
          default_policy,
          default_policy,
          default_policy,
          default_policy,
          default_policy>
{
   typedef policy<detail::forwarding_arg1, detail::forwarding_arg2> type;
};

template <>
struct normalise<policy<detail::forwarding_arg1, detail::forwarding_arg2>,
          promote_float<false>,
          promote_double<false>,
          discrete_quantile<>,
          assert_undefined<>,
          default_policy,
          default_policy,
          default_policy,
          default_policy,
          default_policy,
          default_policy,
          default_policy>
{
   typedef policy<detail::forwarding_arg1, detail::forwarding_arg2> type;
};

inline constexpr policy<> make_policy() noexcept
{ return policy<>(); }

template <class A1>
inline constexpr typename normalise<policy<>, A1>::type make_policy(const A1&) noexcept
{ 
   typedef typename normalise<policy<>, A1>::type result_type;
   return result_type(); 
}

template <class A1, class A2>
inline constexpr typename normalise<policy<>, A1, A2>::type make_policy(const A1&, const A2&) noexcept
{ 
   typedef typename normalise<policy<>, A1, A2>::type result_type;
   return result_type(); 
}

template <class A1, class A2, class A3>
inline constexpr typename normalise<policy<>, A1, A2, A3>::type make_policy(const A1&, const A2&, const A3&) noexcept
{ 
   typedef typename normalise<policy<>, A1, A2, A3>::type result_type;
   return result_type(); 
}

template <class A1, class A2, class A3, class A4>
inline constexpr typename normalise<policy<>, A1, A2, A3, A4>::type make_policy(const A1&, const A2&, const A3&, const A4&) noexcept
{ 
   typedef typename normalise<policy<>, A1, A2, A3, A4>::type result_type;
   return result_type(); 
}

template <class A1, class A2, class A3, class A4, class A5>
inline constexpr typename normalise<policy<>, A1, A2, A3, A4, A5>::type make_policy(const A1&, const A2&, const A3&, const A4&, const A5&) noexcept
{ 
   typedef typename normalise<policy<>, A1, A2, A3, A4, A5>::type result_type;
   return result_type(); 
}

template <class A1, class A2, class A3, class A4, class A5, class A6>
inline constexpr typename normalise<policy<>, A1, A2, A3, A4, A5, A6>::type make_policy(const A1&, const A2&, const A3&, const A4&, const A5&, const A6&) noexcept
{ 
   typedef typename normalise<policy<>, A1, A2, A3, A4, A5, A6>::type result_type;
   return result_type(); 
}

template <class A1, class A2, class A3, class A4, class A5, class A6, class A7>
inline constexpr typename normalise<policy<>, A1, A2, A3, A4, A5, A6, A7>::type make_policy(const A1&, const A2&, const A3&, const A4&, const A5&, const A6&, const A7&) noexcept
{ 
   typedef typename normalise<policy<>, A1, A2, A3, A4, A5, A6, A7>::type result_type;
   return result_type(); 
}

template <class A1, class A2, class A3, class A4, class A5, class A6, class A7, class A8>
inline constexpr typename normalise<policy<>, A1, A2, A3, A4, A5, A6, A7, A8>::type make_policy(const A1&, const A2&, const A3&, const A4&, const A5&, const A6&, const A7&, const A8&) noexcept
{ 
   typedef typename normalise<policy<>, A1, A2, A3, A4, A5, A6, A7, A8>::type result_type;
   return result_type(); 
}

template <class A1, class A2, class A3, class A4, class A5, class A6, class A7, class A8, class A9>
inline constexpr typename normalise<policy<>, A1, A2, A3, A4, A5, A6, A7, A8, A9>::type make_policy(const A1&, const A2&, const A3&, const A4&, const A5&, const A6&, const A7&, const A8&, const A9&) noexcept
{ 
   typedef typename normalise<policy<>, A1, A2, A3, A4, A5, A6, A7, A8, A9>::type result_type;
   return result_type(); 
}

template <class A1, class A2, class A3, class A4, class A5, class A6, class A7, class A8, class A9, class A10>
inline constexpr typename normalise<policy<>, A1, A2, A3, A4, A5, A6, A7, A8, A9, A10>::type make_policy(const A1&, const A2&, const A3&, const A4&, const A5&, const A6&, const A7&, const A8&, const A9&, const A10&) noexcept
{ 
   typedef typename normalise<policy<>, A1, A2, A3, A4, A5, A6, A7, A8, A9, A10>::type result_type;
   return result_type(); 
}

template <class A1, class A2, class A3, class A4, class A5, class A6, class A7, class A8, class A9, class A10, class A11>
inline constexpr typename normalise<policy<>, A1, A2, A3, A4, A5, A6, A7, A8, A9, A10, A11>::type make_policy(const A1&, const A2&, const A3&, const A4&, const A5&, const A6&, const A7&, const A8&, const A9&, const A10&, const A11&) noexcept
{
   typedef typename normalise<policy<>, A1, A2, A3, A4, A5, A6, A7, A8, A9, A10, A11>::type result_type;
   return result_type();
}

//
// Traits class to handle internal promotion:
//
template <class Real, class Policy>
struct evaluation
{
   typedef Real type;
};

template <class Policy>
struct evaluation<float, Policy>
{
   using type = typename std::conditional<Policy::promote_float_type::value, double, float>::type;
};

template <class Policy>
struct evaluation<double, Policy>
{
   using type = typename std::conditional<Policy::promote_double_type::value, long double, double>::type;
};

template <class Real, class Policy>
struct precision
{
   static_assert((std::numeric_limits<Real>::radix == 2) || ((std::numeric_limits<Real>::is_specialized == 0) || (std::numeric_limits<Real>::digits == 0)),
   "(std::numeric_limits<Real>::radix == 2) || ((std::numeric_limits<Real>::is_specialized == 0) || (std::numeric_limits<Real>::digits == 0))");
#ifndef BOOST_BORLANDC
   using precision_type = typename Policy::precision_type;
   using type = typename std::conditional<
      ((std::numeric_limits<Real>::is_specialized == 0) || (std::numeric_limits<Real>::digits == 0)),
      // Possibly unknown precision:
      precision_type,
      typename std::conditional<
         ((std::numeric_limits<Real>::digits <= precision_type::value) 
         || (Policy::precision_type::value <= 0)),
         // Default case, full precision for RealType:
         digits2< std::numeric_limits<Real>::digits>,
         // User customised precision:
         precision_type
      >::type
   >::type;
#else
   using precision_type = typename Policy::precision_type;
   using digits_t = std::integral_constant<int, std::numeric_limits<Real>::digits>;
   using spec_t = std::integral_constant<bool, std::numeric_limits<Real>::is_specialized>;
   using type = typename std::conditional<
      (spec_t::value == true std::true_type || digits_t::value == 0),
      // Possibly unknown precision:
      precision_type,
      typename std::conditional<
         (digits_t::value <= precision_type::value || precision_type::value <= 0),
         // Default case, full precision for RealType:
         digits2< std::numeric_limits<Real>::digits>,
         // User customised precision:
         precision_type
      >::type
   >::type;
#endif
};

#ifdef BOOST_MATH_USE_FLOAT128

template <class Policy>
struct precision<BOOST_MATH_FLOAT128_TYPE, Policy>
{
   typedef std::integral_constant<int, 113> type;
};

#endif

namespace detail{

template <class T, class Policy>
inline constexpr int digits_imp(std::true_type const&) noexcept
{
   static_assert( std::numeric_limits<T>::is_specialized, "std::numeric_limits<T>::is_specialized");
   typedef typename boost::math::policies::precision<T, Policy>::type p_t;
   return p_t::value;
}

template <class T, class Policy>
inline constexpr int digits_imp(std::false_type const&) noexcept
{
   return tools::digits<T>();
}

} // namespace detail

template <class T, class Policy>
inline constexpr int digits(BOOST_MATH_EXPLICIT_TEMPLATE_TYPE(T)) noexcept
{
   typedef std::integral_constant<bool, std::numeric_limits<T>::is_specialized > tag_type;
   return detail::digits_imp<T, Policy>(tag_type());
}
template <class T, class Policy>
inline constexpr int digits_base10(BOOST_MATH_EXPLICIT_TEMPLATE_TYPE(T)) noexcept
{
   return boost::math::policies::digits<T, Policy>() * 301 / 1000L;
}

template <class Policy>
inline constexpr unsigned long get_max_series_iterations() noexcept
{
   typedef typename Policy::max_series_iterations_type iter_type;
   return iter_type::value;
}

template <class Policy>
inline constexpr unsigned long get_max_root_iterations() noexcept
{
   typedef typename Policy::max_root_iterations_type iter_type;
   return iter_type::value;
}

namespace detail{

template <class T, class Digits, class Small, class Default>
struct series_factor_calc
{
   static T get() noexcept(std::is_floating_point<T>::value)
   {
      return ldexp(T(1.0), 1 - Digits::value);
   }
};

template <class T, class Digits>
struct series_factor_calc<T, Digits, std::true_type, std::true_type>
{
   static constexpr T get() noexcept(std::is_floating_point<T>::value)
   {
      return boost::math::tools::epsilon<T>();
   }
};
template <class T, class Digits>
struct series_factor_calc<T, Digits, std::true_type, std::false_type>
{
   static constexpr T get() noexcept(std::is_floating_point<T>::value)
   {
      return 1 / static_cast<T>(static_cast<std::uintmax_t>(1u) << (Digits::value - 1));
   }
};
template <class T, class Digits>
struct series_factor_calc<T, Digits, std::false_type, std::true_type>
{
   static constexpr T get() noexcept(std::is_floating_point<T>::value)
   {
      return boost::math::tools::epsilon<T>();
   }
};

template <class T, class Policy>
inline constexpr T get_epsilon_imp(std::true_type const&) noexcept(std::is_floating_point<T>::value)
{
   static_assert(std::numeric_limits<T>::is_specialized, "std::numeric_limits<T>::is_specialized");
   static_assert(std::numeric_limits<T>::radix == 2, "std::numeric_limits<T>::radix == 2");

   typedef typename boost::math::policies::precision<T, Policy>::type p_t;
   typedef std::integral_constant<bool, p_t::value <= std::numeric_limits<std::uintmax_t>::digits> is_small_int;
   typedef std::integral_constant<bool, p_t::value >= std::numeric_limits<T>::digits> is_default_value;
   return series_factor_calc<T, p_t, is_small_int, is_default_value>::get();
}

template <class T, class Policy>
inline constexpr T get_epsilon_imp(std::false_type const&) noexcept(std::is_floating_point<T>::value)
{
   return tools::epsilon<T>();
}

} // namespace detail

template <class T, class Policy>
inline constexpr T get_epsilon(BOOST_MATH_EXPLICIT_TEMPLATE_TYPE(T)) noexcept(std::is_floating_point<T>::value)
{
   typedef std::integral_constant<bool, (std::numeric_limits<T>::is_specialized && (std::numeric_limits<T>::radix == 2)) > tag_type;
   return detail::get_epsilon_imp<T, Policy>(tag_type());
}

namespace detail{

template <class A1, 
          class A2, 
          class A3,
          class A4,
          class A5,
          class A6,
          class A7,
          class A8,
          class A9,
          class A10,
          class A11>
char test_is_policy(const policy<A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11>*);
double test_is_policy(...);

template <typename P>
class is_policy_imp
{
public:
   static constexpr bool value = (sizeof(::boost::math::policies::detail::test_is_policy(static_cast<P*>(0))) == sizeof(char));
};

}

template <typename P>
class is_policy
{
public:
   static constexpr bool value = boost::math::policies::detail::is_policy_imp<P>::value;
   using type = std::integral_constant<bool, value>;
};

//
// Helper traits class for distribution error handling:
//
template <class Policy>
struct constructor_error_check
{
   using domain_error_type = typename Policy::domain_error_type;
   using type = typename std::conditional<
      (domain_error_type::value == throw_on_error) || (domain_error_type::value == user_error) || (domain_error_type::value == errno_on_error),
      std::true_type,
      std::false_type>::type;
};

template <class Policy>
struct method_error_check
{
   using domain_error_type = typename Policy::domain_error_type;
   using type = typename std::conditional<
      (domain_error_type::value == throw_on_error) && (domain_error_type::value != user_error),
      std::false_type,
      std::true_type>::type;
};
//
// Does the Policy ever throw on error?
//
template <class Policy>
struct is_noexcept_error_policy
{
   typedef typename Policy::domain_error_type               t1;
   typedef typename Policy::pole_error_type                 t2;
   typedef typename Policy::overflow_error_type             t3;
   typedef typename Policy::underflow_error_type            t4;
   typedef typename Policy::denorm_error_type               t5;
   typedef typename Policy::evaluation_error_type           t6;
   typedef typename Policy::rounding_error_type             t7;
   typedef typename Policy::indeterminate_result_error_type t8;

   static constexpr bool value = 
      ((t1::value != throw_on_error) && (t1::value != user_error)
      && (t2::value != throw_on_error) && (t2::value != user_error)
      && (t3::value != throw_on_error) && (t3::value != user_error)
      && (t4::value != throw_on_error) && (t4::value != user_error)
      && (t5::value != throw_on_error) && (t5::value != user_error)
      && (t6::value != throw_on_error) && (t6::value != user_error)
      && (t7::value != throw_on_error) && (t7::value != user_error)
      && (t8::value != throw_on_error) && (t8::value != user_error));
};

}}} // namespaces

#endif // BOOST_MATH_POLICY_HPP

