// Copyright Nick Thompson 2026.
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0. (See accompanying file
// LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include "math_unit_test.hpp"
#include <boost/math/tools/detail/linear_algebra/cholesky.hpp>
#include <boost/math/tools/detail/linear_algebra/gauss_jordan.hpp>
#include <boost/math/tools/detail/linear_algebra/lu.hpp>
#include <boost/math/tools/detail/linear_algebra/qr.hpp>
#include <boost/math/tools/detail/linear_algebra/symmetric_eigen.hpp>
#include <cmath>
#include <cstddef>
#include <limits>
#include <random>
#include <vector>

namespace la = boost::math::tools::detail;

// Entries uniform in [-1, 1], plus shift on the diagonal; n x n, row-major.
template <class Real>
std::vector<Real> random_matrix(std::size_t n, Real shift, unsigned seed)
{
   std::mt19937_64 gen(seed);
   std::uniform_real_distribution<double> u(-1, 1);
   std::vector<Real> A(n * n);
   for (Real& a : A)
      a = static_cast<Real>(u(gen));
   for (std::size_t i = 0; i < n; ++i)
      A[i * n + i] += shift;
   return A;
}

// A B, or A B^T, for n x n matrices.
template <class Real>
std::vector<Real> mul(std::vector<Real> const& A, std::vector<Real> const& B, std::size_t n, bool transpose_b = false)
{
   std::vector<Real> C(n * n, Real(0));
   for (std::size_t i = 0; i < n; ++i)
      for (std::size_t j = 0; j < n; ++j)
         for (std::size_t k = 0; k < n; ++k)
            C[i * n + j] += A[i * n + k] * (transpose_b ? B[j * n + k] : B[k * n + j]);
   return C;
}

// max |A - B|, where only the lower triangles are compared if lower is set, and B = I if it is empty.
template <class Real>
Real distance(std::vector<Real> const& A, std::vector<Real> const& B, std::size_t n, bool lower = false)
{
   using std::abs;
   Real d = 0;
   for (std::size_t i = 0; i < n; ++i)
      for (std::size_t j = 0; j < (lower ? i + 1 : n); ++j)
         d = (std::max)(d, Real(abs(A[i * n + j] - (B.empty() ? Real(i == j) : B[i * n + j]))));
   return d;
}

template <class Real>
void test_lu()
{
   std::size_t const n = 12;
   Real const eps = std::numeric_limits<Real>::epsilon();
   std::vector<Real> const A = random_matrix<Real>(n, Real(n), 1);
   std::vector<Real> x(n), b(n, Real(0)), lu = A, colmax(n);
   std::vector<std::size_t> piv(n);
   for (std::size_t i = 0; i < n; ++i)
   {
      x[i] = Real(i) - 3;
      for (std::size_t j = 0; j < n; ++j)
         b[i] += A[i * n + j] * (Real(j) - 3);
   }
   CHECK_TRUE(la::lu_factor(lu, piv, colmax, n));
   la::lu_solve(lu, piv, n, b);
   for (std::size_t i = 0; i < n; ++i)
      CHECK_LE(std::abs(b[i] - x[i]), 64 * eps);

   // A zero column is refused whatever eps is; a pivot of 2 eps only by the default:
   lu = A;
   for (std::size_t i = 0; i < n; ++i)
      lu[i * n + 5] = 0;
   std::vector<Real> singular = lu;
   CHECK_TRUE(!la::lu_factor(lu, piv, colmax, n));
   CHECK_TRUE(!la::lu_factor(singular, piv, colmax, n, Real(0)));
   std::vector<Real> near{1, 1, 1, 1 + 2 * eps};
   lu = near;
   CHECK_TRUE(!la::lu_factor(lu, piv, colmax, 2));
   CHECK_TRUE(la::lu_factor(near, piv, colmax, 2, Real(0)));
   std::vector<Real> c{2, 2 + 2 * eps};
   la::lu_solve(near, piv, 2, c);
   CHECK_ULP_CLOSE(Real(1), c[0], 0);
   CHECK_ULP_CLOSE(Real(1), c[1], 0);
}

template <class Real>
void test_gauss_jordan()
{
   std::size_t const n = 12;
   Real const eps = std::numeric_limits<Real>::epsilon();
   std::vector<Real> const A = random_matrix<Real>(n, Real(n), 2);
   std::vector<Real> a = A, X(n * n);
   CHECK_TRUE(la::gauss_jordan_inverse(a.data(), X.data(), n));
   CHECK_LE(distance(mul(A, X, n), std::vector<Real>(), n), 64 * eps);
   a = A;
   for (std::size_t i = 0; i < n; ++i)
      a[i * n + 5] = 0;
   CHECK_TRUE(!la::gauss_jordan_inverse(a.data(), X.data(), n));
}

// n = 150 has more than two blocks of the factorization.
template <class Real>
void test_cholesky(std::size_t n)
{
   Real const eps = std::numeric_limits<Real>::epsilon();
   std::vector<Real> const B = random_matrix<Real>(n, Real(0), 3);
   std::vector<Real> A = mul(B, B, n, true);
   for (std::size_t i = 0; i < n; ++i)
      A[i * n + i] += Real(n);
   std::vector<Real> L = A;
   CHECK_TRUE(la::cholesky(L.data(), n));
   for (std::size_t i = 0; i < n; ++i)
      for (std::size_t j = i + 1; j < n; ++j)
         L[i * n + j] = 0;
   CHECK_LE(distance(mul(L, L, n, true), A, n), 16 * Real(n) * eps);

   // Rows that are those of L already are kept:
   std::vector<Real> partial = A;
   std::copy(L.begin(), L.begin() + static_cast<std::ptrdiff_t>((n / 2) * n), partial.begin());
   CHECK_TRUE(la::cholesky(partial.data(), n, n / 2));
   CHECK_LE(distance(partial, L, n, true), Real(0));

   std::vector<Real> x(n), b(n, Real(0));
   for (std::size_t i = 0; i < n; ++i)
   {
      x[i] = Real(i % 7) - 3;
      for (std::size_t j = 0; j < n; ++j)
         b[i] += A[i * n + j] * (Real(j % 7) - 3);
   }
   std::vector<Real> y = b;
   la::cholesky_solve(L.data(), n, y.data());
   for (std::size_t i = 0; i < n; ++i)
      CHECK_LE(std::abs(y[i] - x[i]), 64 * eps);

   // L^{-1} A = L^T, by the solve with several right-hand sides:
   std::vector<Real> Lt = A;
   la::forward_solve(L.data(), n, Lt.data(), n);
   for (std::size_t i = 0; i < n; ++i)
      for (std::size_t j = i; j < n; ++j)
         CHECK_LE(std::abs(Lt[i * n + j] - L[j * n + i]), 512 * eps);

   std::vector<Real> Z(n * n), Ainv(n * n);
   la::cholesky_inverse(L.data(), n, Z.data(), Ainv.data());
   for (std::size_t i = 0; i < n; ++i)
      for (std::size_t j = i + 1; j < n; ++j)
         Ainv[i * n + j] = Ainv[j * n + i];
   CHECK_LE(distance(mul(A, Ainv, n), std::vector<Real>(), n), 64 * eps);

   // Not positive definite from column n - 2 on: the pivot is reported, and the rows before it are those of L.
   std::vector<Real> indefinite = A;
   indefinite[(n - 2) * n + n - 2] -= 2 * Real(n);
   std::size_t failed = 0;
   CHECK_TRUE(!la::cholesky(indefinite.data(), n, 0, &failed));
   CHECK_EQUAL(failed, n - 2);
   CHECK_LE(indefinite[failed * n + failed], Real(0));
   CHECK_ULP_CLOSE(L[(n - 3) * n + n - 3], indefinite[(n - 3) * n + n - 3], 0);
   CHECK_ULP_CLOSE(L[(n - 2) * n + n - 3], indefinite[(n - 2) * n + n - 3], 0);
}

template <class Real>
void test_symmetric_eigen()
{
   std::size_t const n = 12;
   Real const eps = std::numeric_limits<Real>::epsilon();
   std::vector<Real> const B = random_matrix<Real>(n, Real(0), 4);
   std::vector<Real> const A = mul(B, B, n, true);
   std::vector<Real> V, d;
   CHECK_TRUE(la::symmetric_eigendecomposition(A, n, V, d));
   std::vector<Real> VD = V;
   for (std::size_t i = 0; i < n; ++i)
      for (std::size_t j = 0; j < n; ++j)
         VD[i * n + j] *= d[j];
   CHECK_LE(distance(mul(VD, V, n, true), A, n), 64 * Real(n) * eps);
   CHECK_LE(distance(mul(V, V, n, true), std::vector<Real>(), n), 64 * eps);
}

template <class Real>
void test_householder_qr(bool pivot)
{
   std::size_t const n = 12;
   Real const eps = std::numeric_limits<Real>::epsilon();
   std::vector<Real> const A = random_matrix<Real>(n, Real(0), 5);
   std::vector<Real> R = A, Q(n * n);
   la::householder_qr(R.data(), n, n, n, pivot, Q.data());
   CHECK_LE(distance(mul(Q, Q, n, true), std::vector<Real>(), n), 64 * eps);
   for (std::size_t i = 0; i < n; ++i)
      for (std::size_t j = 0; j < i; ++j)
         R[i * n + j] = 0;
   if (!pivot)
      CHECK_LE(distance(mul(Q, R, n), A, n), 64 * eps);
   else
      for (std::size_t i = 0; i + 1 < n; ++i)
         CHECK_LE(std::abs(R[(i + 1) * n + i + 1]), std::abs(R[i * n + i]) * (1 + 64 * eps));
   // Without Q, the same R:
   std::vector<Real> R2 = A;
   la::householder_qr(R2.data(), n, n, n, pivot, static_cast<Real*>(nullptr));
   for (std::size_t i = 0; i < n; ++i)
      for (std::size_t j = i; j < n; ++j)
         CHECK_ULP_CLOSE(R[i * n + j], R2[i * n + j], 0);
}

template <class Real>
void test_all()
{
   test_lu<Real>();
   test_gauss_jordan<Real>();
   test_cholesky<Real>(12);
   test_cholesky<Real>(150);
   test_symmetric_eigen<Real>();
   test_householder_qr<Real>(false);
   test_householder_qr<Real>(true);
}

int main()
{
   test_all<float>();
   test_all<double>();
   return boost::math::test::report_errors();
}
