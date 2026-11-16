//  (C) Copyright Nick Thompson 2026.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_TOOLS_DETAIL_LINEAR_ALGEBRA_QR_HPP
#define BOOST_MATH_TOOLS_DETAIL_LINEAR_ALGEBRA_QR_HPP

#include <boost/math/tools/config.hpp>
#include <boost/math/tools/assert.hpp>
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

// Householder QR of the row-major m x n matrix A, in place, optionally with column pivoting: min(m, n, r)
// reflections, which leave R in the upper triangle of A. Q, if not null, receives the first r columns of the
// orthogonal factor (m x r, row-major).
template <class Real>
void householder_qr(Real* A, std::size_t m, std::size_t n, std::size_t r, bool pivot, Real* Q)
{
   using std::sqrt;
   BOOST_MATH_ASSERT(!Q || r <= m);
   const std::size_t steps = (std::min)((std::min)(m, n), r);
   std::vector<std::vector<Real>> vs;
   std::vector<Real> beta;
   for (std::size_t k = 0; k < steps; ++k)
   {
      if (pivot)
      {
         std::size_t best = k;
         Real bn = -1;
         for (std::size_t j = k; j < n; ++j)
         {
            Real s = 0;
            for (std::size_t i = k; i < m; ++i) s += A[i * n + j] * A[i * n + j];
            if (s > bn) { bn = s; best = j; }
         }
         if (best != k)
            for (std::size_t i = 0; i < m; ++i) std::swap(A[i * n + k], A[i * n + best]);
      }
      std::vector<Real> v(m, Real(0));
      Real nx = 0;
      for (std::size_t i = k; i < m; ++i) { v[i] = A[i * n + k]; nx += v[i] * v[i]; }
      nx = sqrt(nx);
      Real b = 0;
      if (nx > 0)
      {
         v[k] -= v[k] >= 0 ? -nx : nx;
         Real vv = 0;
         for (std::size_t i = k; i < m; ++i) vv += v[i] * v[i];
         b = vv > 0 ? Real(2) / vv : Real(0);
         for (std::size_t j = k; j < n; ++j)
         {
            Real s = 0;
            for (std::size_t i = k; i < m; ++i) s += v[i] * A[i * n + j];
            s *= b;
            for (std::size_t i = k; i < m; ++i) A[i * n + j] -= s * v[i];
         }
      }
      if (Q)
      {
         vs.push_back(v);
         beta.push_back(b);
      }
   }
   if (!Q) return;
   for (std::size_t i = 0; i < m * r; ++i) Q[i] = 0;
   for (std::size_t j = 0; j < r; ++j) Q[j * r + j] = 1;
   for (std::size_t k = vs.size(); k-- > 0;)
      for (std::size_t j = 0; j < r; ++j)
      {
         Real s = 0;
         for (std::size_t i = k; i < m; ++i) s += vs[k][i] * Q[i * r + j];
         s *= beta[k];
         for (std::size_t i = k; i < m; ++i) Q[i * r + j] -= s * vs[k][i];
      }
}

} // namespace detail
} // namespace tools
BOOST_MATH_NAMESPACE_END

#endif // BOOST_MATH_TOOLS_DETAIL_LINEAR_ALGEBRA_QR_HPP
