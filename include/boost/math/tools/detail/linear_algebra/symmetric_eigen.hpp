//  (C) Copyright Nick Thompson 2026.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_TOOLS_DETAIL_LINEAR_ALGEBRA_SYMMETRIC_EIGEN_HPP
#define BOOST_MATH_TOOLS_DETAIL_LINEAR_ALGEBRA_SYMMETRIC_EIGEN_HPP

#include <boost/math/tools/config.hpp>
#ifndef BOOST_MATH_BUILD_MODULE
#include <cmath>
#include <cstddef>
#include <limits>
#include <utility>
#include <vector>
#endif

BOOST_MATH_NAMESPACE_BEGIN
namespace tools {
namespace detail {

// Computes A = V diag(eigenvalues) V^T for a symmetric n x n matrix A stored row-major,
// by Householder tridiagonalization followed by the implicit QL algorithm;
// see Wilkinson & Reinsch, Handbook for Automatic Computation, Vol. II (tred2 and tql2),
// and Golub & Van Loan, Matrix Computations, 4th ed., Section 8.3.
// V is returned row-major, with the eigenvectors as its columns.
// Returns false if some eigenvalue did not converge within max_iterations QL iterations.
template <typename Real>
bool symmetric_eigendecomposition(std::vector<Real> A, size_t n, std::vector<Real> &V, std::vector<Real> &eigenvalues, size_t max_iterations = 30) {
  using std::abs;
  using std::hypot;
  using std::sqrt;
  using std::max;
  V = std::move(A);
  std::vector<Real> &d = eigenvalues;
  d.resize(n);
  std::vector<Real> e(n, Real(0));
  // Householder reduction to tridiagonal form; d is the diagonal and e the subdiagonal:
  for (size_t j = 0; j < n; ++j) {
    d[j] = V[(n - 1)*n + j];
  }
  for (size_t i = n - 1; i > 0; --i) {
    Real scale = 0;
    Real h = 0;
    for (size_t k = 0; k < i; ++k) {
      scale += abs(d[k]);
    }
    if (scale == 0) {
      e[i] = d[i - 1];
      for (size_t j = 0; j < i; ++j) {
        d[j] = V[(i - 1)*n + j];
        V[i*n + j] = 0;
        V[j*n + i] = 0;
      }
    }
    else {
      for (size_t k = 0; k < i; ++k) {
        d[k] /= scale;
        h += d[k]*d[k];
      }
      Real f = d[i - 1];
      Real g = sqrt(h);
      if (f > 0) {
        g = -g;
      }
      e[i] = scale*g;
      h -= f*g;
      d[i - 1] = f - g;
      for (size_t j = 0; j < i; ++j) {
        e[j] = 0;
      }
      for (size_t j = 0; j < i; ++j) {
        f = d[j];
        V[j*n + i] = f;
        g = e[j] + V[j*n + j]*f;
        for (size_t k = j + 1; k < i; ++k) {
          g += V[k*n + j]*d[k];
          e[k] += V[k*n + j]*f;
        }
        e[j] = g;
      }
      f = 0;
      for (size_t j = 0; j < i; ++j) {
        e[j] /= h;
        f += e[j]*d[j];
      }
      const Real hh = f/(h + h);
      for (size_t j = 0; j < i; ++j) {
        e[j] -= hh*d[j];
      }
      for (size_t j = 0; j < i; ++j) {
        f = d[j];
        g = e[j];
        for (size_t k = j; k < i; ++k) {
          V[k*n + j] -= f*e[k] + g*d[k];
        }
        d[j] = V[(i - 1)*n + j];
        V[i*n + j] = 0;
      }
    }
    d[i] = h;
  }
  // Accumulate the transformations:
  for (size_t i = 0; i + 1 < n; ++i) {
    V[(n - 1)*n + i] = V[i*n + i];
    V[i*n + i] = 1;
    const Real h = d[i + 1];
    if (h != 0) {
      for (size_t k = 0; k <= i; ++k) {
        d[k] = V[k*n + i + 1]/h;
      }
      for (size_t j = 0; j <= i; ++j) {
        Real g = 0;
        for (size_t k = 0; k <= i; ++k) {
          g += V[k*n + i + 1]*V[k*n + j];
        }
        for (size_t k = 0; k <= i; ++k) {
          V[k*n + j] -= g*d[k];
        }
      }
    }
    for (size_t k = 0; k <= i; ++k) {
      V[k*n + i + 1] = 0;
    }
  }
  for (size_t j = 0; j < n; ++j) {
    d[j] = V[(n - 1)*n + j];
    V[(n - 1)*n + j] = 0;
  }
  V[(n - 1)*n + n - 1] = 1;
  // The QL rotations act on columns of V, so apply them to the rows of V^T, which are contiguous:
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < i; ++j) {
      std::swap(V[i*n + j], V[j*n + i]);
    }
  }
  // Implicit QL iterations on the tridiagonal matrix:
  for (size_t i = 1; i < n; ++i) {
    e[i - 1] = e[i];
  }
  e[n - 1] = 0;
  Real f = 0;
  Real tst1 = 0;
  const Real eps = std::numeric_limits<Real>::epsilon();
  for (size_t l = 0; l < n; ++l) {
    // Find a small subdiagonal element:
    tst1 = (max)(tst1, abs(d[l]) + abs(e[l]));
    size_t m = l;
    while (m < n - 1 && !(abs(e[m]) <= eps*tst1)) {
      ++m;
    }
    // If m == l, d[l] is already an eigenvalue; otherwise, iterate:
    size_t iterations = 0;
    while (m > l && abs(e[l]) > eps*tst1) {
      if (iterations++ == max_iterations) {
        break;
      }
      // Compute the implicit shift:
      Real g = d[l];
      Real p = (d[l + 1] - g)/(2*e[l]);
      Real r = hypot(p, Real(1));
      if (p < 0) {
        r = -r;
      }
      d[l] = e[l]/(p + r);
      d[l + 1] = e[l]*(p + r);
      const Real dl1 = d[l + 1];
      Real h = g - d[l];
      for (size_t i = l + 2; i < n; ++i) {
        d[i] -= h;
      }
      f += h;
      p = d[m];
      Real c = 1;
      Real c2 = c;
      Real c3 = c;
      const Real el1 = e[l + 1];
      Real s = 0;
      Real s2 = 0;
      for (size_t i = m; i-- > l;) {
        c3 = c2;
        c2 = c;
        s2 = s;
        g = c*e[i];
        h = c*p;
        r = hypot(p, e[i]);
        e[i + 1] = s*r;
        s = e[i]/r;
        c = p/r;
        p = c*d[i] - s*g;
        d[i + 1] = h + s*(c*g + s*d[i]);
        for (size_t k = 0; k < n; ++k) {
          h = V[(i + 1)*n + k];
          V[(i + 1)*n + k] = s*V[i*n + k] + c*h;
          V[i*n + k] = c*V[i*n + k] - s*h;
        }
      }
      p = -s*s2*c3*el1*e[l]/dl1;
      e[l] = s*p;
      d[l] = c*p;
    }
    if (iterations > max_iterations) {
      return false;
    }
    d[l] += f;
    e[l] = 0;
  }
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < i; ++j) {
      std::swap(V[i*n + j], V[j*n + i]);
    }
  }
  return true;
}

} // namespace detail
} // namespace tools
BOOST_MATH_NAMESPACE_END

#endif // BOOST_MATH_TOOLS_DETAIL_LINEAR_ALGEBRA_SYMMETRIC_EIGEN_HPP
