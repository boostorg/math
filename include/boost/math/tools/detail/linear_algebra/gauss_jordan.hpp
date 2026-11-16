//  (C) Copyright Nick Thompson 2026.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_TOOLS_DETAIL_LINEAR_ALGEBRA_GAUSS_JORDAN_HPP
#define BOOST_MATH_TOOLS_DETAIL_LINEAR_ALGEBRA_GAUSS_JORDAN_HPP

#include <boost/math/tools/config.hpp>
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

// x = a^-1 for the row-major n x n matrix a, which is destroyed, by Gauss-Jordan elimination with partial pivoting;
// false if a is singular (a zero pivot).
template <class Real>
bool gauss_jordan_inverse(Real* a, Real* x, std::size_t n)
{
   using std::abs;
   for (std::size_t i = 0; i < n * n; ++i) x[i] = 0;
   for (std::size_t i = 0; i < n; ++i) x[i * n + i] = 1;
   for (std::size_t c = 0; c < n; ++c)
   {
      std::size_t piv = c;
      for (std::size_t r = c + 1; r < n; ++r)
         if (abs(a[r * n + c]) > abs(a[piv * n + c])) piv = r;
      if (a[piv * n + c] == 0) return false;
      if (piv != c)
         for (std::size_t j = 0; j < n; ++j) { std::swap(a[c * n + j], a[piv * n + j]); std::swap(x[c * n + j], x[piv * n + j]); }
      const Real inv = 1 / a[c * n + c];
      for (std::size_t j = 0; j < n; ++j) { a[c * n + j] *= inv; x[c * n + j] *= inv; }
      for (std::size_t r = 0; r < n; ++r)
      {
         if (r == c || a[r * n + c] == 0) continue;
         const Real f = a[r * n + c];
         for (std::size_t j = 0; j < n; ++j) { a[r * n + j] -= f * a[c * n + j]; x[r * n + j] -= f * x[c * n + j]; }
      }
   }
   return true;
}

} // namespace detail
} // namespace tools
BOOST_MATH_NAMESPACE_END

#endif // BOOST_MATH_TOOLS_DETAIL_LINEAR_ALGEBRA_GAUSS_JORDAN_HPP
