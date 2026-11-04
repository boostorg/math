//  Copyright Nick Thompson 2026.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_FILTERS_DETAIL_WAVELET_BOUNDARY_HPP
#define BOOST_MATH_FILTERS_DETAIL_WAVELET_BOUNDARY_HPP

#ifdef _MSC_VER
#pragma once
#endif

#include <boost/math/tools/config.hpp>

#ifndef BOOST_MATH_BUILD_MODULE
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <vector>
#endif

// Boundary rows of the orthogonal wavelet transform on an interval, one set per level.
//
// At every level the interior rows are the rows of the scaling filter h (length 2p) and of its
// wavelet filter. Near each edge, 2p rows span the edge space E: the vectors supported in a window
// [0, w) orthogonal to every interior row that meets the window. Within E the rows are rotated so that
//  * the p high-pass rows annihilate the coefficient vectors that polynomials of degree < p have at
//    this level (at the finest level, the samples P(k)), and
//  * the p low-pass rows span the rest of E,
// each set chosen closest (orthogonal Procrustes) to the Cohen-Daubechies-Vial edge rows. The low-pass
// rows give the polynomials' coefficient vectors at the next level, so the construction repeats.
// Every level is orthogonal, and polynomials of degree < p have no detail coefficients at any level.
//
// Every step that decides orthogonality is an orthogonal transformation or a projector with a spectral
// gap, so the rows are orthonormal to working precision for all p and levels:
//  * E is the range of the window block M = I - C^T C of the projector onto E (C: the interior rows
//    meeting the window), from Householder QR with column pivoting of M; no rank tolerances.
//  * The class subspace and its complement in E come from one complete Householder QR.
//  * The Procrustes rotations are polar factors of square matrices by scaled Newton and Newton-Schulz,
//    orthogonal even when the matrices are ill conditioned.
//  * The class is kept in the eigenpolynomials of the interior map, so it needs no inverse or rescaling.
//  * The CDV rows come from square-root Smith doubling (no Gram matrix) and a staircase QR of the
//    complement of the low rows in E with the known CDV support ends.
//  * Last, the rows are projected onto E and given one Newton-Schulz step.
//
// Any length: the left edge is the same for every n. The right edge, read in reflected coordinates with
// the reversed filter, sees the interior rows shifted by one sample when the level's length is odd. An odd
// level keeps p scaling rows there and gets p - 1 wavelet rows: its edge space is built the same way from
// the shifted interior rows (dimension 2p - 1), and its targets are the CDV rows shifted by one sample, the
// shortest CDV wavelet dropped. The class frame follows: tau <- tau/2 after an even level, (tau - 1)/2
// after an odd one. So the right-edge rows depend on the parities of all the levels before them.
//
// A. Cohen, I. Daubechies, P. Vial, "Wavelets on the interval and fast wavelet transforms",
// Appl. Comput. Harmon. Anal. 1 (1993) 54-81.
//
// boundary_generator is the only entry point; everything else here is its implementation.

BOOST_MATH_NAMESPACE_BEGIN namespace filters { namespace detail {

// The rows of one edge at one level, as for the left edge of [0, inf): row i of lowpass (highpass) is
// the scaling (wavelet) coefficient at edge position i, i = 0..p-1, over the samples 0..width-1.
template <class Real>
struct wavelet_boundary_level
{
   std::size_t width = 0, nhigh = 0; // nhigh: p, or p - 1 at an odd level of the right edge
   std::vector<Real> lowpass, highpass; // p x width and nhigh x width, row-major
};

namespace wavelet_boundary {

template <class Real>
using matrix = std::vector<std::vector<Real>>;

// The tolerance unit: epsilon, but no finer than the filter tables, which have about 75 significant digits.
template <class Real>
Real tolerance_unit() { return (std::max)(Real(std::numeric_limits<Real>::epsilon()), Real(1e-74L)); }

template <class Real>
matrix<Real> zeros(std::size_t r, std::size_t c) { return matrix<Real>(r, std::vector<Real>(c, Real(0))); }

// a b, or a b^T
template <class Real>
matrix<Real> mul(const matrix<Real>& a, const matrix<Real>& b, bool transpose_b = false)
{
   const std::size_t n = a.size(), m = transpose_b ? b.size() : b[0].size(), inner = transpose_b ? b[0].size() : b.size();
   matrix<Real> c = zeros<Real>(n, m);
   for (std::size_t i = 0; i < n; ++i)
      for (std::size_t j = 0; j < m; ++j)
      {
         Real s = 0;
         for (std::size_t l = 0; l < inner; ++l) s += a[i][l] * (transpose_b ? b[j][l] : b[l][j]);
         c[i][j] = s;
      }
   return c;
}

template <class Real>
matrix<Real> transpose(const matrix<Real>& a)
{
   matrix<Real> t = zeros<Real>(a[0].size(), a.size());
   for (std::size_t i = 0; i < a.size(); ++i)
      for (std::size_t j = 0; j < a[0].size(); ++j) t[j][i] = a[i][j];
   return t;
}

// Gauss-Jordan with partial pivoting; false if a is singular.
template <class Real>
bool inverse(matrix<Real> a, matrix<Real>& x)
{
   using std::abs;
   const std::size_t n = a.size();
   x = zeros<Real>(n, n);
   for (std::size_t i = 0; i < n; ++i) x[i][i] = 1;
   for (std::size_t c = 0; c < n; ++c)
   {
      std::size_t piv = c;
      for (std::size_t r = c + 1; r < n; ++r)
         if (abs(a[r][c]) > abs(a[piv][c])) piv = r;
      if (a[piv][c] == 0) return false;
      std::swap(a[c], a[piv]);
      std::swap(x[c], x[piv]);
      const Real d = a[c][c];
      for (std::size_t j = 0; j < n; ++j) { a[c][j] /= d; x[c][j] /= d; }
      for (std::size_t r = 0; r < n; ++r)
      {
         if (r == c || a[r][c] == 0) continue;
         const Real f = a[r][c];
         for (std::size_t j = 0; j < n; ++j) { a[r][j] -= f * a[c][j]; x[r][j] -= f * x[c][j]; }
      }
   }
   return true;
}

// The first r columns of Q from the Householder QR of A (m x n), optionally with column pivoting.
template <class Real>
matrix<Real> householder_q(matrix<Real> A, std::size_t r, bool pivot)
{
   using std::sqrt;
   const std::size_t m = A.size(), n = A[0].size(), steps = (std::min)((std::min)(m, n), r);
   matrix<Real> vs;
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
            for (std::size_t i = k; i < m; ++i) s += A[i][j] * A[i][j];
            if (s > bn) { bn = s; best = j; }
         }
         if (best != k)
            for (std::size_t i = 0; i < m; ++i) std::swap(A[i][k], A[i][best]);
      }
      std::vector<Real> v(m, Real(0));
      Real nx = 0;
      for (std::size_t i = k; i < m; ++i) { v[i] = A[i][k]; nx += v[i] * v[i]; }
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
            for (std::size_t i = k; i < m; ++i) s += v[i] * A[i][j];
            s *= b;
            for (std::size_t i = k; i < m; ++i) A[i][j] -= s * v[i];
         }
      }
      vs.push_back(v);
      beta.push_back(b);
   }
   matrix<Real> Q = zeros<Real>(m, r);
   for (std::size_t j = 0; j < r; ++j) Q[j][j] = 1;
   for (std::size_t k = vs.size(); k-- > 0;)
      for (std::size_t j = 0; j < r; ++j)
      {
         Real s = 0;
         for (std::size_t i = k; i < m; ++i) s += vs[k][i] * Q[i][j];
         s *= beta[k];
         for (std::size_t i = k; i < m; ++i) Q[i][j] -= s * vs[k][i];
      }
   return Q;
}

// The orthogonal polar factor of a square nonsingular matrix: scaled Newton, then Newton-Schulz
// X <- X (3 I - X^T X) / 2 until X^T X = I to working precision. False if X is singular.
template <class Real>
bool polar(matrix<Real>& X)
{
   using std::abs; using std::sqrt;
   const std::size_t n = X.size();
   const Real eps = tolerance_unit<Real>();
   auto fro = [](const matrix<Real>& A) -> Real { Real s = 0; for (const std::vector<Real>& r : A) for (const Real& v : r) s += v * v; return sqrt(s); };
   for (int it = 0; it < 60; ++it)
   {
      matrix<Real> Xi;
      if (!inverse(X, Xi)) return false;
      const matrix<Real> Y = transpose(Xi);
      const Real g = sqrt(fro(Y) / fro(X));
      Real diff = 0;
      for (std::size_t i = 0; i < n; ++i)
         for (std::size_t j = 0; j < n; ++j)
         {
            const Real v = (g * X[i][j] + Y[i][j] / g) / 2;
            diff = (std::max)(diff, Real(abs(v - X[i][j])));
            X[i][j] = v;
         }
      if (diff < sqrt(eps)) break; // quadratic convergence: Newton-Schulz finishes
   }
   for (int it = 0; it < 4; ++it)
   {
      const matrix<Real> XtX = mul(transpose(X), X);
      Real dev = 0;
      matrix<Real> S = zeros<Real>(n, n);
      for (std::size_t i = 0; i < n; ++i)
         for (std::size_t j = 0; j < n; ++j)
         {
            S[i][j] = (i == j ? Real(3) : Real(0)) - XtX[i][j];
            dev = (std::max)(dev, Real(abs(XtX[i][j] - (i == j ? Real(1) : Real(0)))));
         }
      if (dev < 4 * eps && it > 0) break;
      X = mul(X, S);
      for (std::vector<Real>& r : X) for (Real& v : r) v /= 2;
   }
   return true;
}

// The scaling filter in CDV indexing: phi supported on [-p+1, p], h_i for i in [-p+1, p]; and the wavelet filter.
template <class Real>
struct cdv_filter
{
   const std::vector<Real>& h;
   long N;
   Real lo(long i) const { const long m = i + N - 1; return (m < 0 || m >= 2 * N) ? Real(0) : h[std::size_t(m)]; }
   Real hi(long i) const { return (i % 2 == 0) ? lo(1 - i) : Real(-lo(1 - i)); }
};

// An orthonormal basis B (2K - delta rows) of the edge space E, its width w and the projector M onto it.
// The interior rows k >= K cover [2k - p + 1 - delta, 2k + p - delta]. The window is the first where
// M = I - C^T C is a projector of rank 2K - delta. False if none is found.
template <class Real>
bool edge_space(const cdv_filter<Real>& f, std::size_t K, int delta, matrix<Real>& B, std::size_t& w_out, matrix<Real>& M)
{
   using std::abs;
   const long N = f.N;
   const Real eps = tolerance_unit<Real>();
   const std::size_t r = 2 * K - std::size_t(delta);
   for (std::size_t w = r; w <= 4 * K + 2 * std::size_t(N) + 2; ++w)
   {
      M = zeros<Real>(w, w);
      for (std::size_t i = 0; i < w; ++i) M[i][i] = 1;
      for (long k = long(K); 2 * k - N + 1 - delta < long(w); ++k)
      {
         const long lo = (std::max)(0L, 2 * k - N + 1 - delta), hi = (std::min)(long(w), 2 * k + N + 1 - delta);
         for (int pass = 0; pass < 2; ++pass)
            for (long a = lo; a < hi; ++a)
            {
               const Real ra = pass == 0 ? f.lo(a - 2 * k + delta) : f.hi(a - 2 * k + delta);
               for (long b = lo; b < hi; ++b) M[std::size_t(a)][std::size_t(b)] -= ra * (pass == 0 ? f.lo(b - 2 * k + delta) : f.hi(b - 2 * k + delta));
            }
      }
      Real tr = 0, dev = 0;
      for (std::size_t i = 0; i < w; ++i) tr += M[i][i];
      const matrix<Real> M2 = mul(M, M);
      for (std::size_t i = 0; i < w; ++i)
         for (std::size_t j = 0; j < w; ++j) dev = (std::max)(dev, Real(abs(M2[i][j] - M[i][j])));
      if (abs(tr - Real(r)) > Real(1) / 2 || dev > Real(1000) * Real(w) * eps) continue;
      const matrix<Real> Q = householder_q(M, r, true);              // w x r
      B = transpose(householder_q(mul(M, Q), r, false));             // one subspace step on the projector
      w_out = w;
      return true;
   }
   return false;
}

// Staircase basis: orthonormal rows r_0..r_{m-1} spanning the rows of Y, r_a ending at column ends[a]
// (increasing) with a positive last entry, by Householder reflections from the last row up.
template <class Real>
matrix<Real> staircase(matrix<Real> Y, const std::vector<std::size_t>& ends)
{
   using std::sqrt;
   const std::size_t m = Y.size(), w = Y[0].size();
   for (std::size_t a = m; a-- > 0;)
   {
      const std::size_t j = ends[a];
      for (std::size_t i = 0; i <= a; ++i)
         for (std::size_t c = j + 1; c < w; ++c) Y[i][c] = 0;
      std::vector<Real> v(a + 1);
      Real nx = 0;
      for (std::size_t i = 0; i <= a; ++i) { v[i] = Y[i][j]; nx += v[i] * v[i]; }
      nx = sqrt(nx);
      v[a] += (v[a] >= 0 ? Real(1) : Real(-1)) * nx;
      Real uu = 0;
      for (const Real& x : v) uu += x * x;
      if (uu > 0)
         for (std::size_t c = 0; c <= j; ++c)
         {
            Real d = 0;
            for (std::size_t i = 0; i <= a; ++i) d += v[i] * Y[i][c];
            d = 2 * d / uu;
            for (std::size_t i = 0; i <= a; ++i) Y[i][c] -= d * v[i];
         }
      for (std::size_t i = 0; i < a; ++i) Y[i][j] = 0;
      if (Y[a][j] < 0)
         for (std::size_t c = 0; c <= j; ++c) Y[a][c] = -Y[a][c];
   }
   return Y;
}

// The CDV rows of the left edge (p x w each), the Procrustes targets. The Gram matrix
// G = A G A^T + D D^T of the edge functions is never formed: a factor F (G = F F^T) comes from
// square-root Smith doubling, F <- R^T from the QR of [F, A^k F]^T; the reversed-order Cholesky
// factor of G is R^T from one more QR, applied by forward substitution. The edge wavelets are the
// complement of the low rows in E, in staircase form with the CDV support ends p + 2i.
template <class Real>
bool cdv_left_edge(const cdv_filter<Real>& f, const matrix<Real>& BE, std::size_t w, matrix<Real>& lo, matrix<Real>& wav)
{
   using std::abs; using std::sqrt;
   const long N = f.N;
   const std::size_t n = std::size_t(N);
   // P_r(k) = binomial(N-1-k, r) as a polynomial in k
   auto P = [&](std::size_t r, long k) {
      const Real x = Real(N - 1 - k);
      Real v = 1;
      for (std::size_t j = 0; j < r; ++j) v = v * (x - Real(j)) / Real(j + 1);
      return v;
   };
   auto floordiv = [](long a, long b) { return a >= 0 ? a / b : -((-a + b - 1) / b); };
   matrix<Real> V = zeros<Real>(n, n), Vinv;
   for (std::size_t k = 0; k < n; ++k)
      for (std::size_t s = 0; s < n; ++s) V[k][s] = P(s, long(k));
   if (!inverse(V, Vinv)) return false;
   // refinement E[P_r] = sum_s A[r][s] E_1[P_s] + sum_{l=N}^{3N-2} D[r][l-N] phi_{1,l}
   matrix<Real> A = zeros<Real>(n, n), D = zeros<Real>(n, 2 * n - 1);
   for (std::size_t r = 0; r < n; ++r)
   {
      std::vector<Real> c(n, Real(0));
      for (long l = 0; l < N; ++l)
         for (long k = -floordiv(N - l, 2) - 1; k <= (l + N) / 2 + 1; ++k) c[std::size_t(l)] += P(r, k) * f.lo(l - 2 * k);
      for (std::size_t s = 0; s < n; ++s)
         for (std::size_t k = 0; k < n; ++k) A[r][s] += Vinv[s][k] * c[k];
      for (long l = N; l <= 3 * N - 2; ++l)
         for (long k = (std::max)(0L, -floordiv(N - l, 2)); k <= N - 1; ++k) D[r][std::size_t(l - N)] += P(r, k) * f.lo(l - 2 * k);
   }
   // R^T from the QR of M^T (M: n x c): an n x n lower triangular Rt with M M^T = Rt Rt^T
   auto factor = [&](const matrix<Real>& M) {
      matrix<Real> Mt = transpose(M);
      const std::size_t c = Mt.size();
      for (std::size_t k = 0; k < n; ++k)
      {
         Real nx = 0;
         for (std::size_t i = k; i < c; ++i) nx += Mt[i][k] * Mt[i][k];
         nx = sqrt(nx);
         if (nx == 0) continue;
         std::vector<Real> v(c, Real(0));
         for (std::size_t i = k; i < c; ++i) v[i] = Mt[i][k];
         v[k] -= Mt[k][k] >= 0 ? -nx : nx;
         Real vv = 0;
         for (std::size_t i = k; i < c; ++i) vv += v[i] * v[i];
         for (std::size_t j = k; j < n; ++j)
         {
            Real d = 0;
            for (std::size_t i = k; i < c; ++i) d += v[i] * Mt[i][j];
            d = 2 * d / vv;
            for (std::size_t i = k; i < c; ++i) Mt[i][j] -= d * v[i];
         }
      }
      matrix<Real> Rt = zeros<Real>(n, n);
      for (std::size_t i = 0; i < n; ++i)
         for (std::size_t j = i; j < n; ++j) Rt[j][i] = Mt[i][j];
      return Rt;
   };
   matrix<Real> F = factor(D), Ak = A;
   const Real eps = tolerance_unit<Real>();
   for (int it = 0; it < 64; ++it)
   {
      const matrix<Real> AF = mul(Ak, F);
      Real m = 0, mf = 0;
      for (std::size_t i = 0; i < n; ++i)
         for (std::size_t j = 0; j < n; ++j) { m = (std::max)(m, Real(abs(AF[i][j]))); mf = (std::max)(mf, Real(abs(F[i][j]))); }
      if (m == 0 || m < eps * eps * mf) break;
      matrix<Real> Fn = zeros<Real>(n, 2 * n);
      for (std::size_t i = 0; i < n; ++i)
         for (std::size_t j = 0; j < n; ++j) { Fn[i][j] = F[i][j]; Fn[i][n + j] = AF[i][j]; }
      F = factor(Fn);
      Ak = mul(Ak, Ak);
   }
   // reversed order: G_rev = (J F)(J F)^T = L L^T, L lower triangular with positive diagonal; T^{-1} = J L J
   matrix<Real> JF = zeros<Real>(n, n), Tinv = zeros<Real>(n, n);
   for (std::size_t i = 0; i < n; ++i) JF[i] = F[n - 1 - i];
   matrix<Real> L = factor(JF);
   for (std::size_t j = 0; j < n; ++j)
   {
      if (!(L[j][j] != 0)) return false;
      if (L[j][j] < 0)
         for (std::size_t i = j; i < n; ++i) L[i][j] = -L[i][j];
   }
   for (std::size_t i = 0; i < n; ++i)
      for (std::size_t j = 0; j < n; ++j) Tinv[n - 1 - i][n - 1 - j] = L[i][j];
   const matrix<Real> ATi = mul(A, Tinv);
   // Z = [A T^{-1} | D]; the low rows in r order are J L^{-1} J Z: forward substitution on J Z
   const std::size_t lw = 3 * n - 1;
   matrix<Real> X = zeros<Real>(n, lw);
   for (std::size_t i = 0; i < n; ++i)
   {
      const std::size_t r = n - 1 - i;
      for (std::size_t j = 0; j < n; ++j) X[i][j] = ATi[r][j];
      for (std::size_t j = 0; j < 2 * n - 1; ++j) X[i][n + j] = D[r][j];
   }
   for (std::size_t i = 0; i < n; ++i)
   {
      for (std::size_t k = 0; k < i; ++k)
         for (std::size_t j = 0; j < lw; ++j) X[i][j] -= L[i][k] * X[k][j];
      for (std::size_t j = 0; j < lw; ++j) X[i][j] /= L[i][i];
   }
   // X[i] is the low row at position i, with its edge columns in r order
   lo = zeros<Real>(n, w);
   for (std::size_t i = 0; i < n; ++i)
   {
      for (std::size_t j = 0; j < n; ++j) lo[i][j] = X[i][n - 1 - j];
      for (std::size_t j = n; j < (std::min)(lw, w); ++j) lo[i][j] = X[i][j];
   }
   // edge wavelets: the complement of the low rows inside E, staircase with ends N + 2i
   const matrix<Real> Q = householder_q(transpose(mul(lo, BE, true)), 2 * n, false);
   matrix<Real> Hb = zeros<Real>(n, w);
   for (std::size_t i = 0; i < n; ++i)
      for (std::size_t c = 0; c < w; ++c)
      {
         Real v = 0;
         for (std::size_t k = 0; k < 2 * n; ++k) v += Q[k][n + i] * BE[k][c];
         Hb[i][c] = v;
      }
   std::vector<std::size_t> ends(n);
   for (std::size_t i = 0; i < n; ++i) ends[i] = n + 2 * i;
   wav = staircase(Hb, ends);
   return true;
}

} // namespace wavelet_boundary

// The boundary rows of one edge (the left edge of [0, inf) for the filter h of length 2p, p >= 2; the right
// edge is the left edge of the reversed filter, with the samples reversed). The setup, which does not depend
// on n, is done once; rows(odd) then gives the rows of the levels with the parities odd[0] (finest), odd[1], ...
// Odd levels need with_odd. valid() is false if the setup failed.
template <class Real>
class boundary_generator
{
public:
   using matrix = wavelet_boundary::matrix<Real>;

   boundary_generator(const std::vector<Real>& h, bool with_odd) : h_(h), N_(long(h.size() / 2)), K_(h.size() / 2)
   {
      namespace wb = wavelet_boundary;
      using std::abs;
      const wb::cdv_filter<Real> f{h_, N_};
      matrix Tl, Th;
      if (!wb::edge_space(f, K_, 0, geom_[0].B, geom_[0].w, geom_[0].PE) || !wb::cdv_left_edge(f, geom_[0].B, geom_[0].w, Tl, Th)) return;
      set_targets(geom_[0], Tl, Th, 0);
      if (with_odd)
      {
         if (!wb::edge_space(f, K_, 1, geom_[1].B, geom_[1].w, geom_[1].PE)) return;
         // the CDV rows shifted by one sample, the shortest wavelet dropped
         set_targets(geom_[1], Tl, Th, 1);
         geom_[1].valid = true;
      }
      // The class in the eigenpolynomials q_d of the interior map (S q)(k) = sum_i h_i q(2k + i),
      // S q_d = sqrt(2) 2^d q_d. In u = (k - c0) / s, S is triangular in the monomials u^t with well
      // separated diagonal, so q_d comes from a back substitution. The interior entries of the class are
      // q_d(c + tau) and the edge entries propagate as A <- (low X)[., d] / lambda_d.
      const std::size_t D = K_, w = geom_[0].w;
      c0_ = sc_ = Real(w - 1) / 2;
      matrix S = wb::zeros<Real>(D, D);
      {
         std::vector<Real> mom(D, Real(0)); // mom[j] = sum_i h_i ((c0 + i) / s)^j
         for (long i = -N_ + 1; i <= N_; ++i)
         {
            Real t = 1;
            const Real z = (c0_ + Real(i)) / sc_;
            for (std::size_t j = 0; j < D; ++j) { mom[j] += f.lo(i) * t; t *= z; }
         }
         for (std::size_t t = 0; t < D; ++t)
         {
            Real b = 1, pw = 1; // binomial(t, m), 2^m
            for (std::size_t m = 0; m <= t; ++m)
            {
               S[m][t] = b * pw * mom[t - m];
               b = b * Real(t - m) / Real(m + 1);
               pw *= 2;
            }
         }
      }
      beta_ = wb::zeros<Real>(D, D);
      lambda_.resize(D);
      for (std::size_t d = 0; d < D; ++d)
      {
         lambda_[d] = S[d][d];
         beta_[d][d] = 1;
         for (std::size_t m = d; m-- > 0;)
         {
            Real r = 0;
            for (std::size_t t = m + 1; t <= d; ++t) r += S[m][t] * beta_[d][t];
            beta_[d][m] = r / (lambda_[d] - S[m][m]);
         }
      }
      nu_.assign(D, Real(1)); // q_d normalized on the even window
      for (std::size_t d = 0; d < D; ++d)
      {
         Real mx = 0;
         for (std::size_t c = 0; c < w; ++c) mx = (std::max)(mx, Real(abs(q(d, Real(c)))));
         nu_[d] = mx;
      }
      geom_[0].valid = true;
   }

   bool valid() const { return geom_[0].valid; }
   bool valid_odd() const { return geom_[1].valid; }
   std::size_t width(bool odd) const { return geom_[odd ? 1 : 0].w; }

   // The class after the levels generated so far: its edge entries and frame shift.
   struct state
   {
      matrix A;
      Real tau = 0;
      bool ok = true;
   };
   state start() const
   {
      state st;
      st.A = wavelet_boundary::zeros<Real>(K_, K_);
      for (std::size_t i = 0; i < K_; ++i)
         for (std::size_t d = 0; d < K_; ++d) st.A[i][d] = qn(d, Real(i));
      st.ok = valid();
      return st;
   }

   // Appends to out the rows of the levels with these parities (finest first), continuing from st;
   // false (and st.ok false) if a step fails.
   bool extend(state& st, const std::vector<bool>& odd, std::vector<wavelet_boundary_level<Real>>& out) const
   {
      namespace wb = wavelet_boundary;
      const std::size_t K = K_, D = K_;
      matrix& A = st.A;
      Real& tau = st.tau;
      for (bool o : odd)
      {
         if (!st.ok) return false;
         const geometry& g = geom_[o ? 1 : 0];
         if (!g.valid) return st.ok = false;
         const std::size_t w = g.w, dim = g.B.size(), nh = dim - D;
         matrix X = wb::zeros<Real>(w, D);
         for (std::size_t c = 0; c < w; ++c)
            for (std::size_t d = 0; d < D; ++d) X[c][d] = c < K ? A[c][d] : qn(d, Real(c) + tau);
         const matrix Q = wb::householder_q(wb::mul(g.B, X), dim, false); // [O | O_perp]
         matrix O = wb::zeros<Real>(dim, D), Op = wb::zeros<Real>(dim, nh);
         for (std::size_t i = 0; i < dim; ++i)
         {
            for (std::size_t j = 0; j < D; ++j) O[i][j] = Q[i][j];
            for (std::size_t j = D; j < dim; ++j) Op[i][j - D] = Q[i][j];
         }
         matrix PL = wb::mul(g.TlE, O), PH = wb::mul(g.ThE, Op);
         if (!wb::polar(PL) || !wb::polar(PH)) return st.ok = false;
         // (X O^T) B with O^T formed explicitly: the same sums as the reference generator, rounded alike
         matrix low = wb::mul(wb::mul(PL, wb::transpose(O)), g.B), high = wb::mul(wb::mul(PH, wb::transpose(Op)), g.B);
         // R <- R P_E removes the O(eps) leakage onto the interior rows; then R <- (3I - R R^T) R / 2
         matrix R(low);
         R.insert(R.end(), high.begin(), high.end());
         R = wb::mul(R, g.PE);
         const matrix RRt = wb::mul(R, R, true);
         matrix C = wb::zeros<Real>(dim, dim);
         for (std::size_t i = 0; i < dim; ++i)
            for (std::size_t j = 0; j < dim; ++j) C[i][j] = ((i == j ? Real(1) : Real(0)) - RRt[i][j]) / 2;
         const matrix CR = wb::mul(C, R);
         for (std::size_t i = 0; i < dim; ++i)
            for (std::size_t c = 0; c < w; ++c) (i < K ? low[i][c] : high[i - K][c]) = R[i][c] + CR[i][c];
         wavelet_boundary_level<Real> rows;
         rows.width = w;
         rows.nhigh = nh;
         for (const std::vector<Real>& r : low) rows.lowpass.insert(rows.lowpass.end(), r.begin(), r.end());
         for (const std::vector<Real>& r : high) rows.highpass.insert(rows.highpass.end(), r.begin(), r.end());
         out.push_back(std::move(rows));
         A = wb::mul(low, X);
         for (std::size_t i = 0; i < K; ++i)
            for (std::size_t d = 0; d < D; ++d) A[i][d] /= lambda_[d];
         tau = o ? Real((tau - 1) / 2) : Real(tau / 2);
      }
      return st.ok;
   }

   // The rows for the levels with these parities, finest first; empty if a step fails.
   std::vector<wavelet_boundary_level<Real>> rows(const std::vector<bool>& odd) const
   {
      std::vector<wavelet_boundary_level<Real>> out;
      state st = start();
      if (!extend(st, odd, out)) out.clear();
      return out;
   }

private:
   struct geometry
   {
      bool valid = false;
      std::size_t w = 0;
      matrix B, PE, TlE, ThE; // E (rows), its projector, the targets in E coordinates
   };

   // The targets in E coordinates: the rows first.. of T, shifted left by `shift` samples, in g's window.
   static matrix window(const matrix& T, std::size_t first, std::size_t shift, std::size_t w)
   {
      matrix o;
      for (std::size_t i = first; i < T.size(); ++i)
      {
         std::vector<Real> r(w, Real(0));
         for (std::size_t c = 0; c < w && c + shift < T[i].size(); ++c) r[c] = T[i][c + shift];
         o.push_back(r);
      }
      return o;
   }
   void set_targets(geometry& g, const matrix& tl, const matrix& th, std::size_t shift)
   {
      g.TlE = wavelet_boundary::mul(window(tl, 0, shift, g.w), g.B, true);
      g.ThE = wavelet_boundary::mul(window(th, shift, shift, g.w), g.B, true);
   }
   Real q(std::size_t d, const Real& x) const
   {
      const Real u = (x - c0_) / sc_;
      Real v = 0;
      for (std::size_t m = K_; m-- > 0;) v = v * u + beta_[d][m];
      return v;
   }
   Real qn(std::size_t d, const Real& x) const { return q(d, x) / nu_[d]; }

   std::vector<Real> h_;
   long N_;
   std::size_t K_;
   geometry geom_[2]; // even, odd
   matrix beta_;
   std::vector<Real> lambda_, nu_;
   Real c0_ = 0, sc_ = 1;
};

}} BOOST_MATH_NAMESPACE_END

#endif
