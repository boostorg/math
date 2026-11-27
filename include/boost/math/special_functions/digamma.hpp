//  (C) Copyright John Maddock 2006.
//  (C) Copyright Matt Borland 2024.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_SF_DIGAMMA_HPP
#define BOOST_MATH_SF_DIGAMMA_HPP

#ifdef _MSC_VER
#pragma once
#pragma warning(push)
#pragma warning(disable:4702) // Unreachable code (release mode only warning)
#endif

#include <boost/math/tools/config.hpp>
#include <boost/math/tools/type_traits.hpp>
#include <boost/math/tools/rational.hpp>
#include <boost/math/tools/promotion.hpp>
#include <boost/math/policies/policy.hpp>
#include <boost/math/policies/error_handling.hpp>
#include <boost/math/constants/constants.hpp>

#ifndef BOOST_MATH_HAS_NVRTC
#include <boost/math/special_functions/math_fwd.hpp>
#include <boost/math/special_functions/bernoulli.hpp>
#include <boost/math/tools/series.hpp>
#include <boost/math/policies/error_handling.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/math/tools/big_constant.hpp>
#endif

#ifndef BOOST_MATH_HAS_GPU_SUPPORT
#include <boost/math/tools/compensated_arithmetic.hpp>
#include <boost/math/tools/numeric_limits.hpp>
#include <boost/math/tools/precision.hpp>
#endif

#if defined(__GNUC__) && defined(BOOST_MATH_USE_FLOAT128)
//
// This is the only way we can avoid
// warning: non-standard suffix on floating constant [-Wpedantic]
// when building with -Wall -pedantic.  Neither __extension__
// nor #pragma diagnostic ignored work :(
//
#pragma GCC system_header
#endif

BOOST_MATH_NAMESPACE_BEGIN
namespace detail{
//
// Begin by defining the smallest value for which it is safe to
// use the asymptotic expansion for digamma:
//
BOOST_MATH_GPU_ENABLED inline unsigned digamma_large_lim(const BOOST_MATH_NAMESPACE::integral_constant<int, 0>*)
{  return 20;  }
BOOST_MATH_GPU_ENABLED inline unsigned digamma_large_lim(const BOOST_MATH_NAMESPACE::integral_constant<int, 113>*)
{  return 20;  }
BOOST_MATH_GPU_ENABLED inline unsigned digamma_large_lim(const void*)
{  return 10;  }
//
// Implementations of the asymptotic expansion come next,
// the coefficients of the series have been evaluated
// in advance at high precision, and the series truncated
// at the first term that's too small to effect the result.
// Note that the series becomes divergent after a while
// so truncation is very important.
//
// This first one gives 34-digit precision for x >= 20:
//

#ifndef BOOST_MATH_HAS_NVRTC
template <class T>
inline T digamma_imp_large(T x, const BOOST_MATH_NAMESPACE::integral_constant<int, 113>*)
{
   BOOST_MATH_STD_USING // ADL of std functions.
   static const T P[] = {
      BOOST_MATH_BIG_CONSTANT(T, 113, 0.083333333333333333333333333333333333333333333333333),
      BOOST_MATH_BIG_CONSTANT(T, 113, -0.0083333333333333333333333333333333333333333333333333),
      BOOST_MATH_BIG_CONSTANT(T, 113, 0.003968253968253968253968253968253968253968253968254),
      BOOST_MATH_BIG_CONSTANT(T, 113, -0.0041666666666666666666666666666666666666666666666667),
      BOOST_MATH_BIG_CONSTANT(T, 113, 0.0075757575757575757575757575757575757575757575757576),
      BOOST_MATH_BIG_CONSTANT(T, 113, -0.021092796092796092796092796092796092796092796092796),
      BOOST_MATH_BIG_CONSTANT(T, 113, 0.083333333333333333333333333333333333333333333333333),
      BOOST_MATH_BIG_CONSTANT(T, 113, -0.44325980392156862745098039215686274509803921568627),
      BOOST_MATH_BIG_CONSTANT(T, 113, 3.0539543302701197438039543302701197438039543302701),
      BOOST_MATH_BIG_CONSTANT(T, 113, -26.456212121212121212121212121212121212121212121212),
      BOOST_MATH_BIG_CONSTANT(T, 113, 281.4601449275362318840579710144927536231884057971),
      BOOST_MATH_BIG_CONSTANT(T, 113, -3607.510546398046398046398046398046398046398046398),
      BOOST_MATH_BIG_CONSTANT(T, 113, 54827.583333333333333333333333333333333333333333333),
      BOOST_MATH_BIG_CONSTANT(T, 113, -974936.82385057471264367816091954022988505747126437),
      BOOST_MATH_BIG_CONSTANT(T, 113, 20052695.796688078946143462272494530559046688078946),
      BOOST_MATH_BIG_CONSTANT(T, 113, -472384867.72162990196078431372549019607843137254902),
      BOOST_MATH_BIG_CONSTANT(T, 113, 12635724795.916666666666666666666666666666666666667)
   };
   x -= 1;
   T result = log(x);
   result += 1 / (2 * x);
   T z = 1 / (x*x);
   result -= z * tools::evaluate_polynomial(P, z);
   return result;
}
//
// 19-digit precision for x >= 10:
//
template <class T>
inline T digamma_imp_large(T x, const BOOST_MATH_NAMESPACE::integral_constant<int, 64>*)
{
   BOOST_MATH_STD_USING // ADL of std functions.
   static const T P[] = {
      BOOST_MATH_BIG_CONSTANT(T, 64, 0.083333333333333333333333333333333333333333333333333),
      BOOST_MATH_BIG_CONSTANT(T, 64, -0.0083333333333333333333333333333333333333333333333333),
      BOOST_MATH_BIG_CONSTANT(T, 64, 0.003968253968253968253968253968253968253968253968254),
      BOOST_MATH_BIG_CONSTANT(T, 64, -0.0041666666666666666666666666666666666666666666666667),
      BOOST_MATH_BIG_CONSTANT(T, 64, 0.0075757575757575757575757575757575757575757575757576),
      BOOST_MATH_BIG_CONSTANT(T, 64, -0.021092796092796092796092796092796092796092796092796),
      BOOST_MATH_BIG_CONSTANT(T, 64, 0.083333333333333333333333333333333333333333333333333),
      BOOST_MATH_BIG_CONSTANT(T, 64, -0.44325980392156862745098039215686274509803921568627),
      BOOST_MATH_BIG_CONSTANT(T, 64, 3.0539543302701197438039543302701197438039543302701),
      BOOST_MATH_BIG_CONSTANT(T, 64, -26.456212121212121212121212121212121212121212121212),
      BOOST_MATH_BIG_CONSTANT(T, 64, 281.4601449275362318840579710144927536231884057971),
   };
   x -= 1;
   T result = log(x);
   result += 1 / (2 * x);
   T z = 1 / (x*x);
   result -= z * tools::evaluate_polynomial(P, z);
   return result;
}
#endif
//
// 17-digit precision for x >= 10:
//
template <class T>
BOOST_MATH_GPU_ENABLED inline T digamma_imp_large(T x, const BOOST_MATH_NAMESPACE::integral_constant<int, 53>*)
{
   BOOST_MATH_STD_USING // ADL of std functions.
   BOOST_MATH_STATIC const T P[] = {
      0.083333333333333333333333333333333333333333333333333,
      -0.0083333333333333333333333333333333333333333333333333,
      0.003968253968253968253968253968253968253968253968254,
      -0.0041666666666666666666666666666666666666666666666667,
      0.0075757575757575757575757575757575757575757575757576,
      -0.021092796092796092796092796092796092796092796092796,
      0.083333333333333333333333333333333333333333333333333,
      -0.44325980392156862745098039215686274509803921568627
   };
   x -= 1;
   T result = log(x);
   result += 1 / (2 * x);
   T z = 1 / (x*x);
   result -= z * tools::evaluate_polynomial(P, z);
   return result;
}
//
// 9-digit precision for x >= 10:
//
template <class T>
BOOST_MATH_GPU_ENABLED inline T digamma_imp_large(T x, const BOOST_MATH_NAMESPACE::integral_constant<int, 24>*)
{
   BOOST_MATH_STD_USING // ADL of std functions.
   BOOST_MATH_STATIC const T P[] = {
      0.083333333333333333333333333333333333333333333333333f,
      -0.0083333333333333333333333333333333333333333333333333f,
      0.003968253968253968253968253968253968253968253968254f
   };
   x -= 1;
   T result = log(x);
   result += 1 / (2 * x);
   T z = 1 / (x*x);
   result -= z * tools::evaluate_polynomial(P, z);
   return result;
}

#ifndef BOOST_MATH_HAS_NVRTC
//
// Fully generic asymptotic expansion in terms of Bernoulli numbers, see:
// http://functions.wolfram.com/06.14.06.0012.01
//
// LCOV_EXCL_START muliprecision only.
template <class T>
struct digamma_series_func
{
private:
   int k;
   T xx;
   T term;
public:
   digamma_series_func(T x) : k(1), xx(x * x), term(1 / (x * x)) {}
   T operator()()
   {
      T result = term * BOOST_MATH_NAMESPACE::bernoulli_b2n<T>(k) / (2 * k);
      term /= xx;
      ++k;
      return result;
   }
   typedef T result_type;
};

template <class T, class Policy>
inline T digamma_imp_large(T x, const Policy& pol, const BOOST_MATH_NAMESPACE::integral_constant<int, 0>*)
{
   BOOST_MATH_STD_USING
   digamma_series_func<T> s(x);
   T result = log(x) - 1 / (2 * x);
   std::uintmax_t max_iter = policies::get_max_series_iterations<Policy>();
   result = BOOST_MATH_NAMESPACE::tools::sum_series(s, BOOST_MATH_NAMESPACE::policies::get_epsilon<T, Policy>(), max_iter, -result);
   result = -result;
   policies::check_series_iterations<T>("boost::math::digamma<%1%>(%1%)", max_iter, pol);
   return result;
}
// LCOV_EXCL_STOP
//
// Now follow rational approximations over the range [1,2].
//
// 35-digit precision:
//
template <class T>
T digamma_imp_1_2(T x, const BOOST_MATH_NAMESPACE::integral_constant<int, 113>*)
{
   //
   // Now the approximation, we use the form:
   //
   // digamma(x) = (x - root) * (Y + R(x-1))
   //
   // Where root is the location of the positive root of digamma,
   // Y is a constant, and R is optimised for low absolute error
   // compared to Y.
   //
   // Max error found at 128-bit long double precision:  5.541e-35
   // Maximum Deviation Found (approximation error):     1.965e-35
   //
   // LCOV_EXCL_START
   static const float Y = 0.99558162689208984375F;

   static const T root1 = T(1569415565) / 1073741824uL;
   static const T root2 = (T(381566830) / 1073741824uL) / 1073741824uL;
   static const T root3 = ((T(111616537) / 1073741824uL) / 1073741824uL) / 1073741824uL;
   static const T root4 = (((T(503992070) / 1073741824uL) / 1073741824uL) / 1073741824uL) / 1073741824uL;
   static const T root5 = BOOST_MATH_BIG_CONSTANT(T, 113, 0.52112228569249997894452490385577338504019838794544e-36);

   static const T P[] = {
      BOOST_MATH_BIG_CONSTANT(T, 113, 0.25479851061131551526977464225335883769),
      BOOST_MATH_BIG_CONSTANT(T, 113, -0.18684290534374944114622235683619897417),
      BOOST_MATH_BIG_CONSTANT(T, 113, -0.80360876047931768958995775910991929922),
      BOOST_MATH_BIG_CONSTANT(T, 113, -0.67227342794829064330498117008564270136),
      BOOST_MATH_BIG_CONSTANT(T, 113, -0.26569010991230617151285010695543858005),
      BOOST_MATH_BIG_CONSTANT(T, 113, -0.05775672694575986971640757748003553385),
      BOOST_MATH_BIG_CONSTANT(T, 113, -0.0071432147823164975485922555833274240665),
      BOOST_MATH_BIG_CONSTANT(T, 113, -0.00048740753910766168912364555706064993274),
      BOOST_MATH_BIG_CONSTANT(T, 113, -0.16454996865214115723416538844975174761e-4),
      BOOST_MATH_BIG_CONSTANT(T, 113, -0.20327832297631728077731148515093164955e-6)
   };
   static const T Q[] = {
      BOOST_MATH_BIG_CONSTANT(T, 113, 1.0),
      BOOST_MATH_BIG_CONSTANT(T, 113, 2.6210924610812025425088411043163287646),
      BOOST_MATH_BIG_CONSTANT(T, 113, 2.6850757078559596612621337395886392594),
      BOOST_MATH_BIG_CONSTANT(T, 113, 1.4320913706209965531250495490639289418),
      BOOST_MATH_BIG_CONSTANT(T, 113, 0.4410872083455009362557012239501953402),
      BOOST_MATH_BIG_CONSTANT(T, 113, 0.081385727399251729505165509278152487225),
      BOOST_MATH_BIG_CONSTANT(T, 113, 0.0089478633066857163432104815183858149496),
      BOOST_MATH_BIG_CONSTANT(T, 113, 0.00055861622855066424871506755481997374154),
      BOOST_MATH_BIG_CONSTANT(T, 113, 0.1760168552357342401304462967950178554e-4),
      BOOST_MATH_BIG_CONSTANT(T, 113, 0.20585454493572473724556649516040874384e-6),
      BOOST_MATH_BIG_CONSTANT(T, 113, -0.90745971844439990284514121823069162795e-11),
      BOOST_MATH_BIG_CONSTANT(T, 113, 0.48857673606545846774761343500033283272e-13),
   };
   // LCOV_EXCL_STOP
   T g = x - root1;
   g -= root2;
   g -= root3;
   g -= root4;
   g -= root5;
   T r = tools::evaluate_polynomial(P, T(x-1)) / tools::evaluate_polynomial(Q, T(x-1));
   T result = g * Y + g * r;

   return result;
}
//
// 19-digit precision:
//
template <class T>
T digamma_imp_1_2(T x, const BOOST_MATH_NAMESPACE::integral_constant<int, 64>*)
{
   //
   // Now the approximation, we use the form:
   //
   // digamma(x) = (x - root) * (Y + R(x-1))
   //
   // Where root is the location of the positive root of digamma,
   // Y is a constant, and R is optimised for low absolute error
   // compared to Y.
   //
   // Max error found at 80-bit long double precision:   5.016e-20
   // Maximum Deviation Found (approximation error):     3.575e-20
   //
   // LCOV_EXCL_START
   static const float Y = 0.99558162689208984375F;

   static const T root1 = T(1569415565) / 1073741824uL;
   static const T root2 = (T(381566830) / 1073741824uL) / 1073741824uL;
   static const T root3 = BOOST_MATH_BIG_CONSTANT(T, 64, 0.9016312093258695918615325266959189453125e-19);

   static const T P[] = {
      BOOST_MATH_BIG_CONSTANT(T, 64, 0.254798510611315515235),
      BOOST_MATH_BIG_CONSTANT(T, 64, -0.314628554532916496608),
      BOOST_MATH_BIG_CONSTANT(T, 64, -0.665836341559876230295),
      BOOST_MATH_BIG_CONSTANT(T, 64, -0.314767657147375752913),
      BOOST_MATH_BIG_CONSTANT(T, 64, -0.0541156266153505273939),
      BOOST_MATH_BIG_CONSTANT(T, 64, -0.00289268368333918761452)
   };
   static const T Q[] = {
      BOOST_MATH_BIG_CONSTANT(T, 64, 1.0),
      BOOST_MATH_BIG_CONSTANT(T, 64, 2.1195759927055347547),
      BOOST_MATH_BIG_CONSTANT(T, 64, 1.54350554664961128724),
      BOOST_MATH_BIG_CONSTANT(T, 64, 0.486986018231042975162),
      BOOST_MATH_BIG_CONSTANT(T, 64, 0.0660481487173569812846),
      BOOST_MATH_BIG_CONSTANT(T, 64, 0.00298999662592323990972),
      BOOST_MATH_BIG_CONSTANT(T, 64, -0.165079794012604905639e-5),
      BOOST_MATH_BIG_CONSTANT(T, 64, 0.317940243105952177571e-7)
   };
   // LCOV_EXCL_STOP
   T g = x - root1;
   g -= root2;
   g -= root3;
   T r = tools::evaluate_polynomial(P, T(x-1)) / tools::evaluate_polynomial(Q, T(x-1));
   T result = g * Y + g * r;

   return result;
}

#endif
//
// 18-digit precision:
//
template <class T>
BOOST_MATH_GPU_ENABLED T digamma_imp_1_2(T x, const BOOST_MATH_NAMESPACE::integral_constant<int, 53>*)
{
   //
   // Now the approximation, we use the form:
   //
   // digamma(x) = (x - root) * (Y + R(x-1))
   //
   // Where root is the location of the positive root of digamma,
   // Y is a constant, and R is optimised for low absolute error
   // compared to Y.
   //
   // Maximum Deviation Found:               1.466e-18
   // At double precision, max error found:  2.452e-17
   //
   // LCOV_EXCL_START
   BOOST_MATH_STATIC const float Y = 0.99558162689208984F;

   BOOST_MATH_STATIC const T root1 = T(1569415565) / 1073741824uL;
   BOOST_MATH_STATIC const T root2 = (T(381566830) / 1073741824uL) / 1073741824uL;
   BOOST_MATH_STATIC const T root3 = BOOST_MATH_BIG_CONSTANT(T, 53, 0.9016312093258695918615325266959189453125e-19);

   BOOST_MATH_STATIC const T P[] = {
      BOOST_MATH_BIG_CONSTANT(T, 53, 0.25479851061131551),
      BOOST_MATH_BIG_CONSTANT(T, 53, -0.32555031186804491),
      BOOST_MATH_BIG_CONSTANT(T, 53, -0.65031853770896507),
      BOOST_MATH_BIG_CONSTANT(T, 53, -0.28919126444774784),
      BOOST_MATH_BIG_CONSTANT(T, 53, -0.045251321448739056),
      BOOST_MATH_BIG_CONSTANT(T, 53, -0.0020713321167745952)
   };
   BOOST_MATH_STATIC const T Q[] = {
      BOOST_MATH_BIG_CONSTANT(T, 53, 1.0),
      BOOST_MATH_BIG_CONSTANT(T, 53, 2.0767117023730469),
      BOOST_MATH_BIG_CONSTANT(T, 53, 1.4606242909763515),
      BOOST_MATH_BIG_CONSTANT(T, 53, 0.43593529692665969),
      BOOST_MATH_BIG_CONSTANT(T, 53, 0.054151797245674225),
      BOOST_MATH_BIG_CONSTANT(T, 53, 0.0021284987017821144),
      BOOST_MATH_BIG_CONSTANT(T, 53, -0.55789841321675513e-6)
   };
   // LCOV_EXCL_STOP
   T g = x - root1;
   g -= root2;
   g -= root3;
   T r = tools::evaluate_polynomial(P, T(x-1)) / tools::evaluate_polynomial(Q, T(x-1));
   T result = g * Y + g * r;

   return result;
}
//
// 9-digit precision:
//
template <class T>
BOOST_MATH_GPU_ENABLED inline T digamma_imp_1_2(T x, const BOOST_MATH_NAMESPACE::integral_constant<int, 24>*)
{
   //
   // Now the approximation, we use the form:
   //
   // digamma(x) = (x - root) * (Y + R(x-1))
   //
   // Where root is the location of the positive root of digamma,
   // Y is a constant, and R is optimised for low absolute error
   // compared to Y.
   //
   // Maximum Deviation Found:              3.388e-010
   // At float precision, max error found:  2.008725e-008
   //
   // LCOV_EXCL_START
   BOOST_MATH_STATIC const float Y = 0.99558162689208984f;
   BOOST_MATH_STATIC const T root = 1532632.0f / 1048576;
   BOOST_MATH_STATIC const T root_minor = static_cast<T>(0.3700660185912626595423257213284682051735604e-6L);
   BOOST_MATH_STATIC const T P[] = {
      0.25479851023250261e0f,
      -0.44981331915268368e0f,
      -0.43916936919946835e0f,
      -0.61041765350579073e-1f
   };
   BOOST_MATH_STATIC const T Q[] = {
      0.1e1f,
      0.15890202430554952e1f,
      0.65341249856146947e0f,
      0.63851690523355715e-1f
   };
   // LCOV_EXCL_STOP
   T g = x - root;
   g -= root_minor;
   T r = tools::evaluate_polynomial(P, T(x-1)) / tools::evaluate_polynomial(Q, T(x-1));
   T result = g * Y + g * r;

   return result;
}

#ifndef BOOST_MATH_HAS_GPU_SUPPORT
//
// digamma at negative x in double-word arithmetic.  Near each of its zeros there, one between every pair of
// poles, digamma(1 - x) and the cotangent of the reflection cancel, and only their rounding errors are left.
//
template <class T>
class digamma_double_word
{
public:
   using arithmetic = double_word_arithmetic<T>;
   using word = double_word<T>;

   digamma_double_word()
   {
      BOOST_MATH_STD_USING
      // Relative to the largest: a term below drop is left out, and one below sloppy needs only working precision.
      roundoff = ldexp(T(1), -tools::digits<T>());
      full = roundoff / 256;
      const T sloppy {full};
      const T drop {sloppy * roundoff};
      const word one {T(1), T(0)};
      // pi - hi = sin(hi), for hi = pi rounded
      pi = word{constants::pi<T>(), sin(constants::pi<T>())};
      // (-1)^k / (2k)! and (-1)^k / (2k + 1)!, for arguments up to pi / 4
      word factorial {one};
      T bound {1};
      for (unsigned k {0}; (bound > drop) && (k < capacity); ++k)
      {
         if (k > 0)
         {
            const T step {static_cast<T>((2 * k - 1) * (2 * k))};
            factorial = arithmetic::negate(a.divide(factorial, step));
            bound *= T(0.62) / step;
         }
         push(cosine_series, factorial);
         push(sine_series, a.divide(factorial, static_cast<T>(2 * k + 1)));
      }
      // 1 / (2k + 1), for squared arguments up to (sqrt(2) - 1)^2 / (sqrt(2) + 1)^2
      bound = 1;
      for (unsigned k {0}; (bound > drop) && (k < capacity); ++k, bound *= T(0.0295))
      {
         push(atanh_series, a.divide(one, static_cast<T>(2 * k + 1)));
      }
      // log(2) = 7 log(16 / 15) + 5 log(25 / 24) + 3 log(81 / 80)
      ln2 = a.add(a.add(a.multiply(log_ratio(T(31), sloppy, drop), T(7)), a.multiply(log_ratio(T(49), sloppy, drop), T(5))),
                  a.multiply(log_ratio(T(161), sloppy, drop), T(3)));
      // B_2n / (2n) for n = 1 to 17: the last that are ratios of integers exact in a double.
      constexpr double numerators[] { 1, -1, 1, -1, 1, -691, 1, -3617, 43867, -174611, 77683, -236364091, 657931, -3392780147.0, 1723168255201.0, -7709321041217.0, 151628697551.0 };
      constexpr double denominators[] { 12, 120, 252, 240, 132, 32760, 12, 8160, 14364, 6600, 276, 65520, 12, 3480, 85932, 16320, 12 };
      for (unsigned n {0}; n < 17; ++n)
      {
         const T hi {static_cast<T>(numerators[n])};
         const word numerator {hi, static_cast<T>(numerators[n] - static_cast<double>(hi))};
         push(bernoulli_series, a.divide(numerator, static_cast<T>(denominators[n])));
      }
      // The first term left out of that series at w: B_36 / (36 w^36)
      for (unsigned w {1}; w < truncations; ++w)
      {
         truncation[w] = T(380879311253.0) * pow(static_cast<T>(w), T(-36));
      }
   }

   // digamma(x) for x < 0, not an integer, where the result is smaller than its largest term by the factor cancellation.
   T operator()(const T& x, const T& cancellation) const
   {
      // No more of the terms are taken in double words than the cancellation calls for.
      const T wanted {1 / (4096 * cancellation)};
      const T sloppy {wanted > full ? wanted : full};
      const T drop {sloppy * roundoff};
      const auto result {a.add(positive(arithmetic::two_sum(T(1), T(-x)), sloppy, drop), reflection(x, sloppy, drop))};
      return result.hi + result.lo;
   }

private:
   enum { capacity = 32, truncations = 40 };

   // The coefficients of a power series, and their sizes relative to the first.
   struct series
   {
      word coefficients[capacity];
      T sizes[capacity];
      unsigned size = 0;
   };

   static void push(series& s, const word& coefficient)
   {
      BOOST_MATH_STD_USING
      s.coefficients[s.size] = coefficient;
      s.sizes[s.size] = fabs(coefficient.hi / s.coefficients[0].hi);
      ++s.size;
   }

   // Horner's rule: in working precision for the terms too small to need more, then in double words.
   word evaluate(const series& s, const word& z, const T& sloppy, const T& drop) const
   {
      BOOST_MATH_STD_USING
      const T size_of_z {fabs(z.hi)};
      unsigned split {s.size};
      unsigned size {s.size};
      T power {1};
      for (unsigned k {1}; k < s.size; ++k)
      {
         power *= size_of_z;
         const T bound {power * s.sizes[k]};
         if ((bound < sloppy) && (split == s.size))
         {
            split = k;
         }
         if (bound < drop)
         {
            size = k;
            break;
         }
      }
      unsigned k {size - 1};
      word result {s.coefficients[k]};
      if (k >= split)
      {
         T tail {result.hi};
         while (k > split)
         {
            --k;
            tail = tail * z.hi + s.coefficients[k].hi;
         }
         result = word{tail, T(0)};
      }
      while (k > 0)
      {
         --k;
         result = a.add(a.multiply(result, z), s.coefficients[k]);
      }
      return result;
   }

   // 2 atanh(s) = log((1 + s) / (1 - s))
   word twice_atanh(const word& s, const T& sloppy, const T& drop) const
   {
      const auto half {a.multiply(s, evaluate(atanh_series, a.multiply(s, s), sloppy, drop))};
      return word{2 * half.hi, 2 * half.lo};
   }

   // log((m + 1) / m) = 2 atanh(1 / (2m + 1)), given odd = 2m + 1.
   word log_ratio(const T& odd, const T& sloppy, const T& drop) const
   {
      return twice_atanh(a.divide(word{T(1), T(0)}, odd), sloppy, drop);
   }

   // log(v) = e log(2) + 2 atanh((m - 1) / (m + 1)), with v = m 2^e and m in [1 / sqrt(2), sqrt(2)).
   word log_of(const word& v, const T& sloppy, const T& drop) const
   {
      BOOST_MATH_STD_USING
      int e {};
      T hi {frexp(v.hi, &e)};
      if (hi < T(0.70710678118654752440L))
      {
         hi *= 2;
         --e;
      }
      // hi / v.hi is the power of two that scales v to m
      const word m {hi, v.lo * (hi / v.hi)};
      return a.add(twice_atanh(a.divide(a.add(m, T(-1)), a.add(m, T(1))), sloppy, drop), a.multiply(ln2, static_cast<T>(e)));
   }

   // pi cot(pi (1 - x)) = -pi cot(pi d), where d = x less the nearest integer is exact.
   word reflection(const T& x, const T& sloppy, const T& drop) const
   {
      BOOST_MATH_STD_USING
      T d {x};
      if (x <= T(-0.5))
      {
         d = x - floor(x);
         if (d > T(0.5))
         {
            d -= 1;
         }
      }
      // cot(pi r) = tan(pi (1/2 - r)) keeps the argument of the series within pi / 4
      const T r {fabs(d)};
      const bool complement {r > T(0.25)};
      const auto t {a.multiply(pi, complement ? T(T(0.5) - r) : r)};
      const auto t2 {a.multiply(t, t)};
      const auto s {a.multiply(t, evaluate(sine_series, t2, sloppy, drop))};
      const auto c {evaluate(cosine_series, t2, sloppy, drop)};
      const auto result {a.multiply(complement ? a.divide(s, c) : a.divide(c, s), pi)};
      return d < 0 ? result : arithmetic::negate(result);
   }

   // digamma(w) for w >= 1: shifted up to where the asymptotic series is good to drop.
   word positive(word w, const T& sloppy, const T& drop) const
   {
      unsigned limit {2};
      while ((truncation[limit] > drop) && (limit < truncations - 1))
      {
         ++limit;
      }
      word shift {T(0), T(0)};
      while (w.hi < static_cast<T>(limit))
      {
         shift = a.add(shift, reciprocal(w));
         w = a.add(w, T(1));
      }
      // log(w) - 1 / (2w) - sum B_2n / (2n w^2n)
      const auto r {reciprocal(w)};
      const auto z {a.multiply(r, r)};
      const auto tail {a.add(a.multiply(evaluate(bernoulli_series, z, sloppy, drop), z), word{r.hi / 2, r.lo / 2})};
      return a.add(a.add(log_of(w, sloppy, drop), arithmetic::negate(tail)), arithmetic::negate(shift));
   }

   // 1 / w from one division: with q = fl(1 / w.hi), 1 / w = q + q (1 - w q) to double-word precision.
   word reciprocal(const word& w) const
   {
      const T q {1 / w.hi};
      const auto residual {a.add(a.multiply(w, T(-q)), T(1))};
      return arithmetic::fast_two_sum(q, residual.hi * q);
   }

   arithmetic a;
   T roundoff;
   T full;
   word pi;
   word ln2;
   T truncation[truncations];
   series sine_series;
   series cosine_series;
   series atanh_series;
   series bernoulli_series;
};

template <class T>
inline T digamma_by_double_word(const T& x, const T& cancellation)
{
   static const digamma_double_word<T> evaluator;
   return evaluator(x, cancellation);
}

// Built-in binary types from float to double: a long double's fma is in software and far too slow, a narrower
// type cannot hold the tables, and a class that only emulates one of these has not been verified here.
template <class T>
struct digamma_has_double_word : public BOOST_MATH_NAMESPACE::integral_constant<bool,
   BOOST_MATH_NAMESPACE::is_floating_point<T>::value && (BOOST_MATH_NAMESPACE::numeric_limits<T>::radix == 2)
   && (BOOST_MATH_NAMESPACE::numeric_limits<T>::digits >= 24) && (BOOST_MATH_NAMESPACE::numeric_limits<T>::digits <= 53)> {};

// Recomputes a result at x < 0 that is small beside largest, the largest of the terms that were added to give it.
template <class T, class Tag>
inline T digamma_cancellation(const T& x, const T& result, const T& largest, const Tag*, const BOOST_MATH_NAMESPACE::true_type&)
{
   BOOST_MATH_STD_USING
   // An evaluation type wider than the result has bits to lose
   const T limit {(tools::digits<T>() > Tag::value) ? T(32) : T(2)};
   return (fabs(result) * limit < largest) ? digamma_by_double_word(x, T(largest / fabs(result))) : result;
}

template <class T, class Tag>
inline T digamma_cancellation(const T&, const T& result, const T&, const Tag*, const BOOST_MATH_NAMESPACE::false_type&)
{
   return result;
}
#endif // BOOST_MATH_HAS_GPU_SUPPORT

template <class T, class Tag, class Policy>
BOOST_MATH_GPU_ENABLED T digamma_imp(T x, const Tag* t, const Policy& pol)
{
   //
   // This handles reflection of negative arguments, and all our
   // error handling, then forwards to the T-specific approximation.
   //
   BOOST_MATH_STD_USING // ADL of std functions.

   T result = 0;
#ifndef BOOST_MATH_HAS_GPU_SUPPORT
   const T original {x};
#endif
   //
   // Check for negative arguments and use reflection:
   //
   if(x <= -1)
   {
      // Argument reduction for tan uses x itself: 1 - x is not always exact, x - floor(x) is.
      T remainder = x - floor(x);
      //
      // check for evaluation at a negative pole:
      //
      if(remainder == 0)
      {
         return policies::raise_pole_error<T>("boost::math::digamma<%1%>(%1%)", nullptr, x, pol);
      }
      // Reflect:
      x = 1 - x;
      if(remainder < T(0.25))
      {
         // pi cot(pi (1 - x)) = -pi cot(pi remainder)
         result = -constants::pi<T>() / tan(constants::pi<T>() * remainder);
      }
      else if(remainder > T(0.75))
      {
         result = constants::pi<T>() / tan(constants::pi<T>() * (1 - remainder));
      }
      else
      {
         // tan loses digits near pi / 2, so here cot(pi remainder) = tan(pi (1/2 - remainder)), whose argument is exact
         result = -constants::pi<T>() * tan(constants::pi<T>() * (T(0.5) - remainder));
      }
   }
   if(x == 0)
      return policies::raise_pole_error<T>("boost::math::digamma<%1%>(%1%)", nullptr, x, pol);
#ifndef BOOST_MATH_HAS_GPU_SUPPORT
   const T reflection {result};
#endif
   //
   // If we're above the lower-limit for the
   // asymptotic expansion then use it:
   //
   #ifndef BOOST_MATH_HAS_NVRTC
   if(x >= digamma_large_lim(t))
   {
      result += digamma_imp_large(x, t);
   }
   else
   #endif
   {
      //
      // If x > 2 reduce to the interval [1,2]:
      //
      while(x > 2)
      {
         x -= 1;
         result += 1/x;
      }
      //
      // If x < 1 use recurrence to shift to > 1:
      //
      while(x < 1)
      {
         result -= 1/x;
         x += 1;
      }
      result += digamma_imp_1_2(x, t);
   }
#ifndef BOOST_MATH_HAS_GPU_SUPPORT
   if(original < 0)
   {
      // The largest term at negative x: the cotangent where x was reflected, 1 / x where it was not.
      const T largest {fabs(original <= -1 ? reflection : T(1 / original))};
      result = digamma_cancellation(original, result, largest, t, digamma_has_double_word<T>());
   }
#endif
   return result;
}

#ifndef BOOST_MATH_HAS_NVRTC

// LCOV_EXCL_START
template <class T, class Policy>
T digamma_imp(T x, const BOOST_MATH_NAMESPACE::integral_constant<int, 0>* t, const Policy& pol)
{
   //
   // This handles reflection of negative arguments, and all our
   // error handling, then forwards to the T-specific approximation.
   //
   // This is covered by our real_concept tests, but these are disabled for
   // code coverage runs for performance reasons.
   //
   BOOST_MATH_STD_USING // ADL of std functions.

   T result = 0;
   //
   // Check for negative arguments and use reflection:
   //
   if(x <= -1)
   {
      // Argument reduction for tan uses x itself: 1 - x is not always exact, x - floor(x) is.
      T remainder = x - floor(x);
      //
      // check for evaluation at a negative pole:
      //
      if(remainder == 0)
      {
         return policies::raise_pole_error<T>("boost::math::digamma<%1%>(%1%)", nullptr, x, pol);
      }
      // Reflect:
      x = 1 - x;
      if(remainder < T(0.25))
      {
         // pi cot(pi (1 - x)) = -pi cot(pi remainder)
         result = -constants::pi<T>() / tan(constants::pi<T>() * remainder);
      }
      else if(remainder > T(0.75))
      {
         result = constants::pi<T>() / tan(constants::pi<T>() * (1 - remainder));
      }
      else
      {
         // tan loses digits near pi / 2, so here cot(pi remainder) = tan(pi (1/2 - remainder)), whose argument is exact
         result = -constants::pi<T>() * tan(constants::pi<T>() * (T(0.5) - remainder));
      }
   }
   if(x == 0)
      return policies::raise_pole_error<T>("boost::math::digamma<%1%>(%1%)", nullptr, x, pol);
   //
   // If we're above the lower-limit for the
   // asymptotic expansion then use it, the
   // limit is a linear interpolation with
   // limit = 10 at 50 bit precision and
   // limit = 250 at 1000 bit precision.
   //
   int lim = 10 + ((tools::digits<T>() - 50) * 240L) / 950;
   T two_x = ldexp(x, 1);
   if(x >= lim)
   {
      result += digamma_imp_large(x, pol, t);
   }
   else if(floor(x) == x)
   {
      //
      // Special case for integer arguments, see
      // http://functions.wolfram.com/06.14.03.0001.01
      //
      // Summed apart from any reflection term: x is also an integer when 1 - x has rounded to one.
      T sum {-constants::euler<T, Policy>()};
      T val = 1;
      while(val < x)
      {
         sum += 1 / val;
         val += 1;
      }
      result += sum;
   }
   else if(floor(two_x) == two_x)
   {
      //
      // Special case for half integer arguments, see:
      // http://functions.wolfram.com/06.14.03.0007.01
      //
      T sum {-2 * constants::ln_two<T, Policy>() - constants::euler<T, Policy>()};
      int n = itrunc(x);
      if(n)
      {
         for(int k = 1; k < n; ++k)
            sum += 1 / T(k);
         for(int k = n; k <= 2 * n - 1; ++k)
            sum += 2 / T(k);
      }
      result += sum;
   }
   else
   {
      //
      // Rescale so we can use the asymptotic expansion:
      //
      while(x < lim)
      {
         result -= 1 / x;
         x += 1;
      }
      result += digamma_imp_large(x, pol, t);
   }
   return result;
}
// LCOV_EXCL_STOP

#endif

} // namespace detail

BOOST_MATH_EXPORT template <class T, class Policy>
BOOST_MATH_GPU_ENABLED inline typename tools::promote_args<T>::type
   digamma(T x, const Policy&)
{
   typedef typename tools::promote_args<T>::type result_type;
   typedef typename policies::evaluation<result_type, Policy>::type value_type;
   typedef typename policies::precision<T, Policy>::type precision_type;
   typedef BOOST_MATH_NAMESPACE::integral_constant<int,
      (precision_type::value <= 0) || (precision_type::value > 113) ? 0 :
      precision_type::value <= 24 ? 24 :
      precision_type::value <= 53 ? 53 :
      precision_type::value <= 64 ? 64 :
      precision_type::value <= 113 ? 113 : 0 > tag_type;
   typedef typename policies::normalise<
      Policy,
      policies::promote_float<false>,
      policies::promote_double<false>,
      policies::discrete_quantile<>,
      policies::assert_undefined<> >::type forwarding_policy;

   return policies::checked_narrowing_cast<result_type, Policy>(detail::digamma_imp(static_cast<value_type>(x), static_cast<const tag_type*>(nullptr), forwarding_policy()), "boost::math::digamma<%1%>(%1%)");
}

BOOST_MATH_EXPORT template <class T>
BOOST_MATH_GPU_ENABLED inline typename tools::promote_args<T>::type
   digamma(T x)
{
   return digamma(x, policies::policy<>());
}

BOOST_MATH_NAMESPACE_END

#ifdef _MSC_VER
#pragma warning(pop)
#endif

#endif

