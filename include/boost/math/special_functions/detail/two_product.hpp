//  (C) Copyright Nick Thompson 2026.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_SPECIAL_FUNCTIONS_DETAIL_TWO_PRODUCT_HPP
#define BOOST_MATH_SPECIAL_FUNCTIONS_DETAIL_TWO_PRODUCT_HPP

#ifdef _MSC_VER
#pragma once
#endif

#include <boost/math/tools/config.hpp>
#include <boost/math/tools/precision.hpp>

#ifndef BOOST_MATH_BUILD_MODULE
#include <cmath>
#endif

BOOST_MATH_NAMESPACE_BEGIN namespace detail {

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
// directly, but having an fma is no guarantee of one: the Boost.Multiprecision backends without a native one
// (cpp_bin_float among them) evaluate it as a * b + c with two roundings, which would make every e zero.  So
// probe once per type: (1 + eps)^2 = 1 + 2 eps + eps^2, and only a fused fma(x, x, -(1 + 2 eps)) returns the
// last term.  Otherwise split the factors in half (Veltkamp) so that the partial products are exact (Dekker).
// The rounded product comes from fma(a, b, 0), which the compiler can't contract into a later addition.
//
template <class T>
class two_product
{
public:
   two_product() : fused(is_fused()), splitter(0)
   {
      BOOST_MATH_STD_USING
      if (!fused)
         splitter = ldexp(T(1), (tools::digits<T>() + 1) / 2) + 1;
   }
   void operator()(const T& a, const T& b, T& p, T& e) const
   {
      if (fused)
      {
         two_product_fma::evaluate(a, b, T(0), p, 0);
         two_product_fma::evaluate(a, b, T(-p), e, 0);
         return;
      }
      p = a * b;
      T ah, al, bh, bl;
      split(a, ah, al);
      split(b, bh, bl);
      e = (((ah * bh - p) + ah * bl) + al * bh) + al * bl;
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
   bool fused;
   T splitter;
};

} BOOST_MATH_NAMESPACE_END // namespace detail

#endif // BOOST_MATH_SPECIAL_FUNCTIONS_DETAIL_TWO_PRODUCT_HPP
