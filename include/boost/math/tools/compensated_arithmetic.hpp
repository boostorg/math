//  (C) Copyright Nick Thompson 2026.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_TOOLS_COMPENSATED_ARITHMETIC_HPP
#define BOOST_MATH_TOOLS_COMPENSATED_ARITHMETIC_HPP

#ifdef _MSC_VER
#pragma once
#endif

#include <boost/math/tools/config.hpp>
#include <boost/math/tools/precision.hpp>

#if !defined(BOOST_MATH_BUILD_MODULE) && !defined(BOOST_MATH_HAS_NVRTC)
#include <cmath>
#endif

//
// Reassociation (-ffast-math, /fp:fast) simplifies every error term below to zero, so forbid it here.  GCC won't
// inline these into fast-math callers, so only pin it there when it's needed.  Clang has the pragma from 12
// (Apple clang 13).
//
#if defined(__clang__) && (__clang_major__ >= (defined(__apple_build_version__) ? 13 : 12))
#define BOOST_MATH_COMPENSATED_ARITHMETIC_CLANG_FP
#endif
#if defined(BOOST_MATH_COMPENSATED_ARITHMETIC_CLANG_FP)
#pragma float_control(push)
#pragma clang fp reassociate(off)
#elif defined(__GNUC__) && defined(__ASSOCIATIVE_MATH__)
#pragma GCC push_options
#pragma GCC optimize("no-associative-math")
#elif defined(_MSC_VER) && !defined(__clang__) && defined(_M_FP_FAST)
#pragma float_control(precise, on, push)
#endif

BOOST_MATH_NAMESPACE_BEGIN namespace detail {

#ifndef BOOST_MATH_HAS_NVRTC

namespace two_product_fma {
using std::fma;

// result = fma(a, b, c) where the type has an fma, found by argument-dependent lookup for the
// non-built-in types.  Returns false, leaving result alone, where there is none.
template <class T>
inline auto evaluate(const T& a, const T& b, const T& c, T& result, int) -> decltype(static_cast<void>(fma(a, b, c)), bool())
{
   result = fma(a, b, c);
   return true;
}
template <class T>
inline bool evaluate(const T&, const T&, const T&, T&, long)
{
   return false;
}
} // namespace two_product_fma

//
// Splits a product into its rounded value and its rounding error: p + e = a b exactly.  A fused fma gives e
// directly, but having an fma is no guarantee of one: some C libraries (newlib, hence Cygwin) and the
// Boost.Multiprecision backends without a native one (cpp_bin_float) evaluate it as a * b + c with two
// roundings, which would make every e zero.  So probe once per type: (1 + eps)^2 = 1 + 2 eps + eps^2, and only
// a fused fma(x, x, -(1 + 2 eps)) returns the last term.  Otherwise split the factors in half (Veltkamp) so that
// the partial products are exact (Dekker).
// The rounded product comes from fma(a, b, 0), which the compiler can't contract into a later addition.
//
template <class T>
class two_product
{
public:
   two_product() : splitter(0)
   {
      BOOST_MATH_STD_USING
      if (!is_fused())
         splitter = ldexp(T(1), (tools::digits<T>() + 1) / 2) + 1;
   }
   void operator()(const T& a, const T& b, T& p, T& e) const
   {
      if (is_fused())
         two_product_fma::evaluate(a, b, T(0), p, 0);
      else
         p = a * b;
      e = error(a, b, p);
   }
   // The rounding error of ab = fl(a b): a b = ab + error(a, b, ab) exactly.
   T error(const T& a, const T& b, const T& ab) const
   {
      if (is_fused())
      {
         T e;
         two_product_fma::evaluate(a, b, T(-ab), e, 0);
         return e;
      }
      T ah, al, bh, bl;
      split(a, ah, al);
      split(b, bh, bl);
      return (((ah * bh - ab) + ah * bl) + al * bh) + al * bl;
   }
private:
   void split(const T& t, T& hi, T& lo) const
   {
      T g = splitter * t;
      hi = g - (g - t);
      lo = t - hi;
   }
   static bool is_fused()
   {
      static const bool result = []
      {
         // From a volatile, so that the compiler can't fold the probe at compile time:
         volatile int volatile_one = 1;
         T one = static_cast<int>(volatile_one);
         T eps = tools::epsilon<T>();
         T x = one + eps;
         T eps_squared = 0;
         return two_product_fma::evaluate(x, x, T(-(one + 2 * eps)), eps_squared, 0) && (eps_squared != 0);
      }();
      return result;
   }
   T splitter;
};

//
// a b - c d without the cancellation: with p1 + e1 = a b and p2 + e2 = c d exactly, p1 - p2 is exact where they
// cancel (Sterbenz's lemma).  Measured near cancellation it was correctly rounded, where Kahan's
// fma(a, b, -p2) - e2 is out by up to 1.5 ulps.
//
template <class T>
inline T difference_of_products(const T& a, const T& b, const T& c, const T& d)
{
   two_product<T> product;
   T p1, e1, p2, e2;
   product(a, b, p1, e1);
   product(c, d, p2, e2);
   return (p1 - p2) + (e1 - e2);
}

#endif // BOOST_MATH_HAS_NVRTC

//
// The exact sum a + b is s + (return value), where s = fl(a + b):
//
template <class T>
BOOST_MATH_GPU_ENABLED inline T two_sum_error(T a, T b, T s)
{
   T bb = s - a;
   return (a - (s - bb)) + (b - bb);
}

//
// Neumaier's compensated sum: sum + comp carries the total to about one rounding, however many terms are added.
//
template <class T>
struct compensated_sum
{
   T sum;
   T comp;
   BOOST_MATH_GPU_ENABLED explicit compensated_sum(T s0 = 0) : sum(s0), comp(0) {}
   BOOST_MATH_GPU_ENABLED void add(T x)
   {
      BOOST_MATH_STD_USING
      T t = sum + x;
      comp += (fabs(sum) >= fabs(x)) ? T((sum - t) + x) : T((x - t) + sum);
      sum = t;
   }
   BOOST_MATH_GPU_ENABLED T value() const { return sum + comp; }
};

} BOOST_MATH_NAMESPACE_END // namespace detail

#if defined(BOOST_MATH_COMPENSATED_ARITHMETIC_CLANG_FP)
#undef BOOST_MATH_COMPENSATED_ARITHMETIC_CLANG_FP
#pragma float_control(pop)
#elif defined(__GNUC__) && defined(__ASSOCIATIVE_MATH__)
#pragma GCC pop_options
#elif defined(_MSC_VER) && !defined(__clang__) && defined(_M_FP_FAST)
#pragma float_control(pop)
#endif

#endif // BOOST_MATH_TOOLS_COMPENSATED_ARITHMETIC_HPP
