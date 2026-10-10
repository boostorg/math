
//  (C) Copyright John Maddock 2006.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_SPECIAL_LAGUERRE_HPP
#define BOOST_MATH_SPECIAL_LAGUERRE_HPP

#ifdef _MSC_VER
#pragma once
#endif

#include <boost/math/special_functions/math_fwd.hpp>
#include <boost/math/tools/config.hpp>
#include <boost/math/policies/error_handling.hpp>

#ifndef BOOST_MATH_BUILD_MODULE
#include <cmath>
#endif

BOOST_MATH_NAMESPACE_BEGIN

// Recurrence relation for Laguerre polynomials:
BOOST_MATH_EXPORT template <class T1, class T2, class T3>
inline typename tools::promote_args<T1, T2, T3>::type  
   laguerre_next(unsigned n, T1 x, T2 Ln, T3 Lnm1)
{
   typedef typename tools::promote_args<T1, T2, T3>::type result_type;
   return ((2 * n + 1 - result_type(x)) * result_type(Ln) - n * result_type(Lnm1)) / (n + 1);
}

namespace detail{

//
// Implement Laguerre polynomials via recurrence.  On its own the recurrence loses up to n^1.5 epsilon or so
// where L_n oscillates, as both of its solutions are of similar size there and the rounding errors of each step
// accumulate.  So we find the rounding error of each step almost exactly, and carry it forward with the same
// recurrence: the result is then within the conditioning of L_n.  The rounded products come from fma(a, b, 0),
// which the compiler can't contract into a later addition, as that would break the error free transformations:
//
template <class T>
T laguerre_imp(unsigned n, T x)
{
   using std::fma;
   T p0 = 1;
   T p1 = 1 - x;

   if(n == 0)
      return p0;

   // Errors in p0 and p1, starting with that of 1 - x:
   T e0 = 0;
   T t = p1 - 1;
   T e1 = (1 - (p1 - t)) + (-x - t);

   for(unsigned c = 1; c < n; ++c)
   {
      // a + ea = 2c + 1 - x exactly:
      T a = (2 * c + 1) - x;
      t = a - (2 * c + 1);
      T ea = ((2 * c + 1) - (a - t)) + (-x - t);
      // s1 + r1 = a p1 and s2 + r2 = c p0 exactly:
      T s1 = fma(a, p1, T(0));
      T r1 = fma(a, p1, -s1);
      T s2 = fma(T(c), p0, T(0));
      T r2 = fma(T(c), p0, -s2);
      // d + ed = s1 - s2 exactly:
      T d = s1 - s2;
      t = d - s1;
      T ed = (s1 - (d - t)) + (-s2 - t);
      T p2 = d / (c + 1);
      // d - (c + 1) p2 exactly:
      T rd = fma(-T(c + 1), p2, d);
      T e2 = (a * e1 - c * e0 + rd + ed + r1 - r2 + ea * p1) / (c + 1);
      p0 = p1;
      p1 = p2;
      e0 = e1;
      e1 = e2;
   }
   return p1 + e1;
}

template <class T, class Policy>
inline typename tools::promote_args<T>::type 
laguerre(unsigned n, T x, const Policy&, const std::true_type&)
{
   typedef typename tools::promote_args<T>::type result_type;
   typedef typename policies::evaluation<result_type, Policy>::type value_type;
   return policies::checked_narrowing_cast<result_type, Policy>(detail::laguerre_imp(n, static_cast<value_type>(x)), "boost::math::laguerre<%1%>(unsigned, %1%)");
}

template <class T>
inline typename tools::promote_args<T>::type 
   laguerre(unsigned n, unsigned m, T x, const std::false_type&)
{
   return BOOST_MATH_NAMESPACE::laguerre(n, m, x, policies::policy<>());
}

} // namespace detail

BOOST_MATH_EXPORT template <class T>
inline typename tools::promote_args<T>::type 
   laguerre(unsigned n, T x)
{
   return laguerre(n, x, policies::policy<>());
}

// Recurrence for associated polynomials:
BOOST_MATH_EXPORT template <class T1, class T2, class T3>
inline typename tools::promote_args<T1, T2, T3>::type  
   laguerre_next(unsigned n, unsigned l, T1 x, T2 Pl, T3 Plm1)
{
   typedef typename tools::promote_args<T1, T2, T3>::type result_type;
   return ((2 * n + l + 1 - result_type(x)) * result_type(Pl) - (n + l) * result_type(Plm1)) / (n+1);
}

namespace detail{
// Laguerre Associated Polynomial:
template <class T, class Policy>
T laguerre_imp(unsigned n, unsigned m, T x, const Policy& pol)
{
   // Special cases:
   if(m == 0)
      return BOOST_MATH_NAMESPACE::laguerre(n, x, pol);

   T p0 = 1;
   
   if(n == 0)
      return p0;

   T p1 = m + 1 - x;

   unsigned c = 1;

   while(c < n)
   {
      std::swap(p0, p1);
      p1 = static_cast<T>(laguerre_next(c, m, x, p0, p1));
      ++c;
   }
   return p1;
}

}

BOOST_MATH_EXPORT template <class T, class Policy>
inline typename tools::promote_args<T>::type 
   laguerre(unsigned n, unsigned m, T x, const Policy& pol)
{
   typedef typename tools::promote_args<T>::type result_type;
   typedef typename policies::evaluation<result_type, Policy>::type value_type;
   return policies::checked_narrowing_cast<result_type, Policy>(detail::laguerre_imp(n, m, static_cast<value_type>(x), pol), "boost::math::laguerre<%1%>(unsigned, unsigned, %1%)");
}

BOOST_MATH_EXPORT template <class T1, class T2>
inline typename laguerre_result<T1, T2>::type 
   laguerre(unsigned n, T1 m, T2 x)
{
   typedef typename policies::is_policy<T2>::type tag_type;
   return detail::laguerre(n, m, x, tag_type());
}

BOOST_MATH_NAMESPACE_END

#endif // BOOST_MATH_SPECIAL_LAGUERRE_HPP



