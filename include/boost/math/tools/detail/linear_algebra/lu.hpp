//  (C) Copyright Nick Thompson 2026.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_TOOLS_DETAIL_LINEAR_ALGEBRA_LU_HPP
#define BOOST_MATH_TOOLS_DETAIL_LINEAR_ALGEBRA_LU_HPP

#include <boost/math/tools/config.hpp>
#include <boost/math/special_functions/fpclassify.hpp>
#include <boost/math/tools/precision.hpp>
#ifndef BOOST_MATH_BUILD_MODULE
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <utility>
#include <vector>
#endif

BOOST_MATH_NAMESPACE_BEGIN
namespace tools {
namespace detail {

// LU with partial pivoting of the row-major n x n matrix a, in place. False when a pivot is not above n eps times
// the largest entry of its column, so the test does not depend on the scaling of the variables. With eps = 0 only
// a zero or non-finite pivot fails.
template <class Real>
bool lu_factor(std::vector<Real>& a, std::vector<std::size_t>& piv, std::vector<Real>& colmax, std::size_t n,
               Real const eps = tools::epsilon<Real>())
{
   using std::abs;
   for (std::size_t j = 0; j < n; ++j)
   {
      colmax[j] = 0;
      for (std::size_t i = 0; i < n; ++i)
         colmax[j] = (std::max)(colmax[j], Real(abs(a[i * n + j])));
   }
   Real const tiny = n * eps;
   for (std::size_t k = 0; k < n; ++k)
   {
      std::size_t p = k;
      for (std::size_t i = k + 1; i < n; ++i)
         if (abs(a[i * n + k]) > abs(a[p * n + k]))
            p = i;
      piv[k] = p;
      if (p != k)
         for (std::size_t j = 0; j < n; ++j)
            std::swap(a[k * n + j], a[p * n + j]);
      Real const d = a[k * n + k];
      if (!(abs(d) > tiny * colmax[k]) || !(BOOST_MATH_NAMESPACE::isfinite)(d))
         return false;
      for (std::size_t i = k + 1; i < n; ++i)
      {
         Real const l = a[i * n + k] /= d;
         for (std::size_t j = k + 1; j < n; ++j)
            a[i * n + j] -= l * a[k * n + j];
      }
   }
   return true;
}

template <class Real, class V>
void lu_solve(std::vector<Real> const& a, std::vector<std::size_t> const& piv, std::size_t n, V& b)
{
   for (std::size_t k = 0; k < n; ++k)
      std::swap(b[k], b[piv[k]]);
   for (std::size_t i = 1; i < n; ++i)
      for (std::size_t j = 0; j < i; ++j)
         b[i] -= a[i * n + j] * b[j];
   for (std::size_t i = n; i-- > 0;)
   {
      for (std::size_t j = i + 1; j < n; ++j)
         b[i] -= a[i * n + j] * b[j];
      b[i] /= a[i * n + i];
   }
}

} // namespace detail
} // namespace tools
BOOST_MATH_NAMESPACE_END

#endif // BOOST_MATH_TOOLS_DETAIL_LINEAR_ALGEBRA_LU_HPP
