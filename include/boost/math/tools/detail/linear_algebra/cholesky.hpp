//  (C) Copyright Nick Thompson 2026.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_TOOLS_DETAIL_LINEAR_ALGEBRA_CHOLESKY_HPP
#define BOOST_MATH_TOOLS_DETAIL_LINEAR_ALGEBRA_CHOLESKY_HPP

#include <boost/math/tools/config.hpp>
#ifndef BOOST_MATH_BUILD_MODULE
#include <algorithm>
#include <cmath>
#include <cstddef>
#endif

BOOST_MATH_NAMESPACE_BEGIN
namespace tools {
namespace detail {

// Hand-rolled dense linear algebra on row-major arrays. The kernels with long inner loops are written so that
// the compiler can vectorize them without reassociating floating point sums: either as axpy updates of contiguous
// rows, or as several rows' dot products with one shared operand, which gives independent accumulators.

// s[r] = sum_{k < m} a[r][k] b[k] for four rows a[0..3].
template <class Real> void dot4(Real const *const *a, Real const *b, size_t m, Real *s) {
  Real s0 = 0, s1 = 0, s2 = 0, s3 = 0, t0 = 0, t1 = 0, t2 = 0, t3 = 0;
  size_t k = 0;
  for (; k + 1 < m; k += 2) {
    Real const x = b[k];
    Real const y = b[k + 1];
    s0 += a[0][k] * x;
    t0 += a[0][k + 1] * y;
    s1 += a[1][k] * x;
    t1 += a[1][k + 1] * y;
    s2 += a[2][k] * x;
    t2 += a[2][k + 1] * y;
    s3 += a[3][k] * x;
    t3 += a[3][k + 1] * y;
  }
  if (k < m) {
    Real const x = b[k];
    s0 += a[0][k] * x;
    s1 += a[1][k] * x;
    s2 += a[2][k] * x;
    s3 += a[3][k] * x;
  }
  s[0] = s0 + t0;
  s[1] = s1 + t1;
  s[2] = s2 + t2;
  s[3] = s3 + t3;
}

template <class Real> Real dot(Real const *a, Real const *b, size_t m) {
  Real s0 = 0, s1 = 0, s2 = 0, s3 = 0;
  size_t k = 0;
  for (; k + 3 < m; k += 4) {
    s0 += a[k] * b[k];
    s1 += a[k + 1] * b[k + 1];
    s2 += a[k + 2] * b[k + 2];
    s3 += a[k + 3] * b[k + 3];
  }
  for (; k < m; ++k) {
    s0 += a[k] * b[k];
  }
  return (s0 + s1) + (s2 + s3);
}

// C[i][j] -= sum_{k < len} X[i][k] X[j][k] for lo <= j <= i < hi, where row i of X starts at X + i ldx and C has
// leading dimension ldc. Tiles of 64 rows j stay in cache while every block of four rows i below them goes by.
// Rows i < first of C are left alone; the others get what they would with first = 0.
template <class Real>
void subtract_gram(Real *C, size_t ldc, Real const *X, size_t ldx, size_t len, size_t lo, size_t hi, size_t first = 0) {
  constexpr size_t tile = 64;
  Real s[4];
  Real const *rows[4];
  for (size_t j0 = lo; j0 < hi; j0 += tile) {
    size_t const j1 = (std::min)(hi, j0 + tile);
    for (size_t i = j0; i < hi; i += 4) {
      size_t const nb = (std::min)(size_t(4), hi - i), skip = first > i ? first - i : 0;
      for (size_t b = 0; b < 4; ++b) {
        rows[b] = X + (i + (b < nb ? b : 0)) * ldx;
      }
      size_t const jend = skip < nb ? (std::min)(j1, i + nb) : j0;
      for (size_t j = j0; j < jend; ++j) {
        Real const *xj = X + j * ldx;
        if (j < i) {
          dot4(rows, xj, len, s);
          for (size_t b = skip; b < nb; ++b) {
            C[(i + b) * ldc + j] -= s[b];
          }
        }
        else {
          for (size_t b = (std::max)(j - i, skip); b < nb; ++b) {
            C[(i + b) * ldc + j] -= dot(rows[b], xj, len);
          }
        }
      }
    }
  }
}

// The Cholesky factorization A = L L^T in place: on input the lower triangle of the n x n matrix A (leading dimension
// n), on output that of L. The upper triangle is not referenced. Returns false if A is not numerically positive
// definite. Right-looking by blocks of 64 columns; within a block, left-looking four rows at a time. If rows i < first
// are those of L already, and the others those of A, the result is that of the whole factorization. On failure
// *failed, if given, is the column j of the pivot that is not positive: A[j][j] then holds that pivot, and the rows
// before it and row j up to the diagonal are those of L.
BOOST_MATH_TEST_EXPORT template <class Real>
bool cholesky(Real *A, size_t n, size_t first = 0, size_t *failed = nullptr) {
  using std::sqrt;
  constexpr size_t block = 64;
  Real s[4];
  Real *rows[4];
  Real const *parts[4];
  for (size_t k0 = 0; k0 < n; k0 += block) {
    size_t const k1 = (std::min)(n, k0 + block);
    // Columns k0 to k1 of L, from the diagonal block down:
    for (size_t i = k0; i < n; i += 4) {
      size_t const nb = (std::min)(size_t(4), n - i), skip = first > i ? first - i : 0;
      for (size_t b = 0; b < 4; ++b) {
        rows[b] = A + (i + (b < nb ? b : 0)) * n;
        parts[b] = rows[b] + k0;
      }
      size_t const jend = skip < nb ? (std::min)(k1, i + nb) : k0;
      for (size_t j = k0; j < jend; ++j) {
        Real const *lj = A + j * n;
        if (j < i) {
          dot4(parts, lj + k0, j - k0, s);
          for (size_t b = skip; b < nb; ++b) {
            rows[b][j] = (rows[b][j] - s[b]) / lj[j];
          }
          continue;
        }
        // Column j meets these rows on the diagonal: row j first, then the rows below it.
        for (size_t b = (std::max)(j - i, skip); b < nb; ++b) {
          Real const d = rows[b][j] - dot(parts[b], lj + k0, j - k0);
          if (b == j - i) {
            if (!(d > 0)) {
              rows[b][j] = d;
              if (failed) {
                *failed = j;
              }
              return false;
            }
            rows[b][j] = sqrt(d);
          }
          else {
            rows[b][j] = d / lj[j];
          }
        }
      }
    }
    // The trailing submatrix loses the contribution of these columns:
    subtract_gram(A, n, A + k0, n, k1 - k0, k1, n, first);
  }
  return true;
}

// Solves L X = B in place for the n x m matrix B (leading dimension m), where L is the lower triangle of an n x n
// matrix (leading dimension n).
BOOST_MATH_TEST_EXPORT template <class Real> void forward_solve(Real const *L, size_t n, Real *B, size_t m) {
  for (size_t i = 0; i < n; ++i) {
    Real *bi = B + i * m;
    for (size_t k = 0; k < i; ++k) {
      Real const lik = L[i * n + k];
      Real const *bk = B + k * m;
      for (size_t c = 0; c < m; ++c) {
        bi[c] -= lik * bk[c];
      }
    }
    Real const inv = 1 / L[i * n + i];
    for (size_t c = 0; c < m; ++c) {
      bi[c] *= inv;
    }
  }
}

// Solves L x = b in place for a vector b.
BOOST_MATH_TEST_EXPORT template <class Real> void forward_solve(Real const *L, size_t n, Real *b) {
  for (size_t i = 0; i < n; ++i) {
    b[i] = (b[i] - dot(L + i * n, b, i)) / L[i * n + i];
  }
}

// Solves L^T x = b in place, by rows of L: once x_i is known, remove its contribution L_ik x_i from b_k for k < i.
BOOST_MATH_TEST_EXPORT template <class Real> void backward_solve(Real const *L, size_t n, Real *b) {
  for (size_t i = n; i-- > 0;) {
    Real const *li = L + i * n;
    Real const xi = b[i] / li[i];
    b[i] = xi;
    for (size_t k = 0; k < i; ++k) {
      b[k] -= li[k] * xi;
    }
  }
}

// Solves L L^T x = b in place.
BOOST_MATH_TEST_EXPORT template <class Real> void cholesky_solve(Real const *L, size_t n, Real *b) {
  forward_solve(L, n, b);
  backward_solve(L, n, b);
}

// The lower triangle of (L L^T)^{-1} = Z^T Z, where Z = L^{-1}. Z is n x n scratch. Four rows of Z are taken at a
// time, for which each row of Z has zeros in the three places after its diagonal.
template <class Real> void cholesky_inverse(Real const *L, size_t n, Real *Z, Real *Kinv) {
  for (size_t i = 0; i < n; ++i) {
    Real *zi = Z + i * n;
    std::fill(zi, zi + (std::min)(n, i + 4), Real(0));
    size_t k = 0;
    for (; k + 3 < i; k += 4) {
      Real const *l = L + i * n + k, *z0 = Z + k * n, *z1 = z0 + n, *z2 = z1 + n, *z3 = z2 + n;
      for (size_t j = 0; j <= k + 3; ++j) {
        zi[j] = zi[j] - l[0] * z0[j] - l[1] * z1[j] - l[2] * z2[j] - l[3] * z3[j];
      }
    }
    for (; k < i; ++k) {
      Real const lik = L[i * n + k];
      Real const *zk = Z + k * n;
      for (size_t j = 0; j <= k; ++j) {
        zi[j] -= lik * zk[j];
      }
    }
    Real const inv = 1 / L[i * n + i];
    for (size_t j = 0; j < i; ++j) {
      zi[j] *= inv;
    }
    zi[i] = inv;
  }
  for (size_t a = 0; a < n; ++a) {
    Real *row = Kinv + a * n;
    std::fill(row, row + a + 1, Real(0));
    size_t k = a;
    for (; k + 3 < n; k += 4) {
      Real const *z0 = Z + k * n, *z1 = z0 + n, *z2 = z1 + n, *z3 = z2 + n;
      Real const c0 = z0[a], c1 = z1[a], c2 = z2[a], c3 = z3[a];
      for (size_t b = 0; b <= a; ++b) {
        row[b] = row[b] + c0 * z0[b] + c1 * z1[b] + c2 * z2[b] + c3 * z3[b];
      }
    }
    for (; k < n; ++k) {
      Real const *zk = Z + k * n;
      Real const zka = zk[a];
      for (size_t b = 0; b <= a; ++b) {
        row[b] += zka * zk[b];
      }
    }
  }
}

} // namespace detail
} // namespace tools
BOOST_MATH_NAMESPACE_END

#endif // BOOST_MATH_TOOLS_DETAIL_LINEAR_ALGEBRA_CHOLESKY_HPP
