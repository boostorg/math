
//  (C) Copyright John Maddock 2006.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_SPECIAL_LAGUERRE_HPP
#define BOOST_MATH_SPECIAL_LAGUERRE_HPP

#include <boost/math/special_functions/math_fwd.hpp>
#include <boost/math/tools/config.hpp>
#include <boost/math/tools/evaluation_type.hpp>

namespace boost{
namespace math{

// Recurrance relation for Laguerre polynomials:
template <class T>
inline T laguerre_next(unsigned n, T x, T Ln, T Lnm1)
{
   return ((2 * n + 1 - x) * Ln - n * Lnm1) / (n + 1);
}

namespace detail{

// Implement Laguerre polynomials via recurrance:
template <class T>
T laguerre_imp(unsigned n, T x)
{
   T p0 = 1;
   T p1 = 1 - x;

   if(n == 0)
      return p0;

   unsigned c = 1;

   while(c < n)
   {
      std::swap(p0, p1);
      p1 = laguerre_next(c, x, p0, p1);
      ++c;
   }
   return p1;
}

} // namespace detail

template <class T>
inline T laguerre(unsigned n, T x)
{
   typedef typename tools::evaluation<typename remove_cv<T>::type>::type value_type;
   return tools::checked_narrowing_cast<typename remove_cv<T>::type>(detail::laguerre_imp(n, static_cast<value_type>(x)), BOOST_CURRENT_FUNCTION);
}

// Recurrence for associated polynomials:
template <class T>
inline T laguerre_next(unsigned n, unsigned l, T x, T Pl, T Plm1)
{
   return ((2 * n + l + 1 - x) * Pl - (n + l) * Plm1) / (n+1);
}

namespace detail{
// Laguerre Associated Polynomial:
template <class T>
T laguerre_imp(unsigned n, unsigned m, T x)
{
   // Special cases:
   if(m == 0)
      return laguerre(n, x);

   T p0 = 1;
   
   if(n == 0)
      return p0;

   T p1 = m + 1 - x;

   unsigned c = 1;

   while(c < n)
   {
      std::swap(p0, p1);
      p1 = laguerre_next(c, m, x, p0, p1);
      ++c;
   }
   return p1;
}

}

template <class T>
inline T laguerre(unsigned n, unsigned m, T x)
{
   typedef typename tools::evaluation<typename remove_cv<T>::type>::type value_type;
   return tools::checked_narrowing_cast<typename remove_cv<T>::type>(detail::laguerre_imp(n, m, static_cast<value_type>(x)), BOOST_CURRENT_FUNCTION);
}

} // namespace math
} // namespace boost

#endif // BOOST_MATH_SPECIAL_LAGUERRE_HPP

