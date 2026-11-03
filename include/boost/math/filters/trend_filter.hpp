//  Copyright Nick Thompson 2026.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_FILTERS_TREND_FILTER_HPP
#define BOOST_MATH_FILTERS_TREND_FILTER_HPP

#ifdef _MSC_VER
#pragma once
#endif

#include <boost/math/tools/config.hpp>

#ifndef BOOST_MATH_BUILD_MODULE
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iterator>
#include <limits>
#include <type_traits>
#include <utility>
#include <vector>
#endif
#include <boost/math/policies/error_handling.hpp>
#include <boost/math/policies/policy.hpp>

// Trend filtering: x = argmin 1/2 ||y - x||^2 + lambda ||D^(k+1) x||_1, D^(k+1) the (k+1)-th difference.
// k = 0 is total-variation denoising (Condat's direct algorithm, exact). k = 1, 2 use the primal-dual
// interior-point method of Kim, Koh, Boyd & Gorinevsky on the dual box-constrained QP, with banded Cholesky.

BOOST_MATH_NAMESPACE_BEGIN namespace filters { namespace detail {

// Condat, "A direct algorithm for 1D total variation denoising", IEEE SPL 20 (2013).
template <class Real>
void trend_filter_tv(const Real* y, Real* x, std::size_t n, const Real& lambda)
{
   std::size_t k = 0, k0 = 0, kplus = 0, kminus = 0;
   Real umin = lambda, umax = -lambda, vmin = y[0] - lambda, vmax = y[0] + lambda;
   const Real twolambda = 2 * lambda, minlambda = -lambda;
   for (;;)
   {
      while (k == n - 1)
      {
         if (umin < 0)
         {
            do { x[k0++] = vmin; } while (k0 <= kminus);
            umax = (vmin = y[kminus = k = k0]) + (umin = lambda) - vmax;
         }
         else if (umax > 0)
         {
            do { x[k0++] = vmax; } while (k0 <= kplus);
            umin = (vmax = y[kplus = k = k0]) + (umax = minlambda) - vmin;
         }
         else
         {
            vmin += umin / Real(k - k0 + 1);
            do { x[k0++] = vmin; } while (k0 <= k);
            return;
         }
      }
      if ((umin += y[k + 1] - vmin) < minlambda)
      {
         do { x[k0++] = vmin; } while (k0 <= kminus);
         vmax = (vmin = y[kplus = kminus = k = k0]) + twolambda;
         umin = lambda; umax = minlambda;
      }
      else if ((umax += y[k + 1] - vmax) > lambda)
      {
         do { x[k0++] = vmax; } while (k0 <= kplus);
         vmin = (vmax = y[kplus = kminus = k = k0]) - twolambda;
         umin = lambda; umax = minlambda;
      }
      else
      {
         ++k;
         if (umin >= lambda) { vmin += (umin - lambda) / Real((kminus = k) - k0 + 1); umin = lambda; }
         if (umax <= minlambda) { vmax += (umax + lambda) / Real((kplus = k) - k0 + 1); umax = minlambda; }
      }
   }
}

// The R-th difference: (D x)_i = sum_j d_j x_{i+j}, d_j = (-1)^(R-j) binomial(R, j); D is (n-R) x n.
template <class Real, int R>
struct difference
{
   Real d[R + 1];
   difference()
   {
      long b = 1;
      for (int j = 0; j <= R; ++j) { d[j] = Real((R - j) % 2 ? -b : b); b = b * (R - j) / (j + 1); }
   }
   void apply(const Real* x, Real* out, std::size_t n) const
   {
      for (std::size_t i = 0; i + R < n; ++i)
      {
         Real s = 0;
         for (int j = 0; j <= R; ++j) s += d[j] * x[i + j];
         out[i] = s;
      }
   }
   // out = D^T z, z of length n - R
   void apply_transpose(const Real* z, Real* out, std::size_t n) const
   {
      const std::size_t m = n - R;
      for (std::size_t i = 0; i < n; ++i)
      {
         Real s = 0;
         const std::size_t j0 = i >= m ? i - m + 1 : 0, j1 = (std::min)(i, std::size_t(R));
         for (std::size_t j = j0; j <= j1; ++j) s += d[j] * z[i - j];
         out[i] = s;
      }
   }
};

// Symmetric positive definite band matrix of half-bandwidth B; row i holds A(i, i-B..i). Cholesky in place.
template <class Real, int B>
struct band_matrix
{
   std::size_t m = 0;
   std::vector<Real> a, inv_diag;
   Real& at(std::size_t i, int j) { return a[i * (B + 1) + j]; }
   const Real& at(std::size_t i, int j) const { return a[i * (B + 1) + j]; }
   bool cholesky()
   {
      using std::sqrt;
      inv_diag.resize(m);
      for (std::size_t i = 0; i < m; ++i)
      {
         const int j0 = i >= std::size_t(B) ? 0 : int(B - i);
         for (int jj = j0; jj <= B; ++jj)
         {
            const std::size_t j = i - B + jj;
            Real s = at(i, jj);
            for (int kk = j0; kk < jj; ++kk)
            {
               const std::size_t c = i - B + kk;
               if (c + B >= j) s -= at(i, kk) * at(j, int(c + B - j));
            }
            if (jj == B)
            {
               if (!(s > 0)) return false;
               at(i, B) = sqrt(s);
               inv_diag[i] = 1 / at(i, B);
            }
            else
               at(i, jj) = s * inv_diag[j];
         }
      }
      return true;
   }
   void solve(Real* x) const
   {
      for (std::size_t i = 0; i < m; ++i)
      {
         Real s = x[i];
         for (int jj = i >= std::size_t(B) ? 0 : int(B - i); jj < B; ++jj) s -= at(i, jj) * x[i - B + jj];
         x[i] = s * inv_diag[i];
      }
      for (std::size_t i = m; i-- > 0;)
      {
         Real s = x[i];
         for (int q = 1; q <= B && i + q < m; ++q) s -= at(i + q, B - q) * x[i + q];
         x[i] = s * inv_diag[i];
      }
   }
};

// Primal-dual interior point (Kim, Koh, Boyd & Gorinevsky, SIAM Review 51 (2009)) on the dual
// min 1/2 ||D^T z||^2 - y^T D^T z subject to |z_i| <= lambda; then x = y - D^T z. R = k + 1.
// Stops when the duality gap is <= tol * scale. At large lambda the gap cannot be evaluated that
// accurately (D x = D y - D D^T z cancels terms of size lambda), so the method also stops when the
// surrogate gap sum mu .* (-f) is <= tol * scale and a step of at least 1/2 moved x by less than
// tol * scale in squared norm (Boyd & Vandenberghe, Convex Optimization, 11.7.3).
// With z of size lambda, x = y - D^T z is only known to about eps lambda, so no representable z has a
// gap below about eps lambda 4^R sum |z_i|. When that floor exceeds tol * scale (large lambda and n) and
// the gap has stopped falling (no new minimum in 10 iterations), the method stops once the gap and the
// surrogate gap are both within it.
// Returns false, with x filled by the error handler's result, if it fails.
template <class Real, int R, class Policy>
bool trend_filter_pdip(const std::vector<Real>& y, const Real& lambda, const Real& scale, std::vector<Real>& x, const Policy& pol)
{
   using std::abs; using std::sqrt;
   static const char* function = "boost::math::filters::trend_filter<%1%>";
   const std::size_t n = y.size(), m = n - R;
   const difference<Real, R> D;
   const Real eps = std::numeric_limits<Real>::epsilon();
   const Real tol = sqrt(eps); // relative duality gap: |x - x*|^2 <= 2 tol scale
   const Real alpha = Real(1) / 100, beta = Real(1) / 2, mu_factor = 2;
   auto fail = [&](const char* message) {
      x.assign(n, policies::raise_evaluation_error(function, message, lambda, pol));
      return false;
   };

   std::vector<Real> Dy(m), z(m, Real(0)), mu1(m, Real(1)), mu2(m, Real(1)), f1(m, -lambda), f2(m, -lambda), DTz(n), DDTz(m),
       dz(m), dmu1(m), dmu2(m), diag(m), nz(m), nmu1(m), nmu2(m), wn(n), rb(m);
   D.apply(y.data(), Dy.data(), n);

   // D D^T is banded Toeplitz: c[s] = sum_j d_j d_{j+s}. S = D D^T + diag(mu1/(-f1) + mu2/(-f2)).
   Real c[R + 1];
   for (int s = 0; s <= R; ++s)
   {
      c[s] = 0;
      for (int j = 0; j + s <= R; ++j) c[s] += D.d[j + s] * D.d[j];
   }
   band_matrix<Real, R> S;
   S.m = m;
   S.a.assign(m * (R + 1), Real(0));
   auto fill = [&](const Real& shift) {
      for (std::size_t i = 0; i < m; ++i)
      {
         for (int s = 0; s <= R; ++s) S.at(i, R - s) = std::size_t(s) <= i ? c[s] : Real(0);
         S.at(i, R) += diag[i] + shift;
      }
   };

   // norm of the residual of the perturbed KKT conditions at (zz, m1, m2) for barrier parameter inv_t
   auto residual = [&](const std::vector<Real>& zz, const std::vector<Real>& m1, const std::vector<Real>& m2, const Real& inv_t) {
      D.apply_transpose(zz.data(), wn.data(), n);
      D.apply(wn.data(), rb.data(), n);
      Real s = 0;
      for (std::size_t i = 0; i < m; ++i)
      {
         const Real rd = rb[i] - Dy[i] + m1[i] - m2[i];
         const Real c1 = m1[i] * (lambda - zz[i]) - inv_t, c2 = m2[i] * (lambda + zz[i]) - inv_t;
         s += rd * rd + c1 * c1 + c2 * c2;
      }
      return sqrt(s);
   };

   Real t = Real(1) / 10000000000, step = (std::numeric_limits<Real>::max)(), dx2 = (std::numeric_limits<Real>::max)();
   std::vector<Real> gaps;
   for (int it = 0; it < 200; ++it)
   {
      D.apply_transpose(z.data(), DTz.data(), n);
      D.apply(DTz.data(), DDTz.data(), n);
      Real p = 0, dual = 0, dtz2 = 0, eta = 0;
      for (std::size_t i = 0; i < m; ++i)
      {
         p += abs(Dy[i] - DDTz[i]);
         dual += Dy[i] * z[i];
         eta += mu1[i] * (lambda - z[i]) + mu2[i] * (lambda + z[i]);
      }
      for (std::size_t i = 0; i < n; ++i) dtz2 += DTz[i] * DTz[i];
      // primal objective at x = y - D^T z minus the dual objective
      const Real gap = dtz2 + lambda * p - dual;
      Real zsum = 0;
      for (std::size_t i = 0; i < m; ++i) zsum += abs(z[i]);
      const Real gap_floor = eps * lambda * Real(1 << (2 * R)) * zsum;
      gaps.push_back(gap);
      const bool stalled = it >= 10 && *std::min_element(gaps.end() - 10, gaps.end()) >= *std::min_element(gaps.begin(), gaps.end() - 10);
      if (gap <= tol * scale || (eta <= tol * scale && dx2 <= tol * scale && step >= beta && step <= 1)
          || (stalled && gap_floor > tol * scale && (std::max)(gap, eta) <= tol * scale + gap_floor))
      {
         x.resize(n);
         for (std::size_t i = 0; i < n; ++i) x[i] = y[i] - DTz[i];
         return true;
      }
      if (step >= Real(1) / 5) t = (std::max)(2 * Real(m) * mu_factor / (std::max)(gap, (std::numeric_limits<Real>::min)()), t * 6 / 5);
      const Real inv_t = 1 / t;
      for (std::size_t i = 0; i < m; ++i)
      {
         diag[i] = -(mu1[i] / f1[i] + mu2[i] / f2[i]);
         dz[i] = Dy[i] - DDTz[i] + inv_t / f1[i] - inv_t / f2[i];
      }
      fill(Real(0));
      if (!S.cholesky())
      {
         // definiteness lost to rounding only: retry with a relative shift of the diagonal
         Real mx = 0;
         for (std::size_t i = 0; i < m; ++i) mx = (std::max)(mx, c[0] + diag[i]);
         fill(4096 * eps * mx);
         if (!S.cholesky())
            return fail("Newton system not positive definite at lambda = %1%.");
      }
      S.solve(dz.data());
      Real s = 1;
      for (std::size_t i = 0; i < m; ++i)
      {
         dmu1[i] = -(mu1[i] + (inv_t + dz[i] * mu1[i]) / f1[i]);
         dmu2[i] = -(mu2[i] + (inv_t - dz[i] * mu2[i]) / f2[i]);
         if (dmu1[i] < 0) s = (std::min)(s, -mu1[i] / dmu1[i] * 99 / 100);
         if (dmu2[i] < 0) s = (std::min)(s, -mu2[i] / dmu2[i] * 99 / 100);
      }
      const Real r0 = residual(z, mu1, mu2, inv_t);
      bool accepted = false;
      for (int ls = 0; ls < 60 && !accepted; ++ls, s *= beta)
      {
         bool feasible = true;
         for (std::size_t i = 0; i < m; ++i)
         {
            nz[i] = z[i] + s * dz[i];
            nmu1[i] = mu1[i] + s * dmu1[i];
            nmu2[i] = mu2[i] + s * dmu2[i];
            if (!(abs(nz[i]) < lambda)) feasible = false;
         }
         accepted = feasible && residual(nz, nmu1, nmu2, inv_t) <= (1 - alpha * s) * r0;
         if (accepted) step = s;
      }
      if (!accepted)
         return fail("Line search failed at lambda = %1%.");
      for (std::size_t i = 0; i < m; ++i) rb[i] = nz[i] - z[i];
      D.apply_transpose(rb.data(), wn.data(), n);
      dx2 = 0;
      for (std::size_t i = 0; i < n; ++i) dx2 += wn[i] * wn[i];
      z.swap(nz); mu1.swap(nmu1); mu2.swap(nmu2);
      for (std::size_t i = 0; i < m; ++i) { f1[i] = z[i] - lambda; f2[i] = -z[i] - lambda; }
   }
   return fail("No convergence in 200 interior-point iterations at lambda = %1%.");
}

// The least-squares fit of y by polynomials of degree <= k (modified Gram-Schmidt, twice, on
// monomials in t in [-1, 1]). D^(k+1) annihilates it, so the trend filter of y is this fit plus
// the trend filter of y minus it.
template <class Real>
std::vector<Real> polynomial_fit(const std::vector<Real>& y, int k)
{
   using std::sqrt;
   const std::size_t n = y.size();
   std::vector<Real> p(n, Real(0)), q(n);
   std::vector<std::vector<Real>> basis;
   for (int j = 0; j <= k; ++j)
   {
      for (std::size_t i = 0; i < n; ++i)
      {
         const Real t = n > 1 ? Real(2 * i) / Real(n - 1) - 1 : Real(0);
         q[i] = j == 0 ? Real(1) : j == 1 ? t : t * t;
      }
      for (int pass = 0; pass < 2; ++pass)
         for (const std::vector<Real>& b : basis)
         {
            Real s = 0;
            for (std::size_t i = 0; i < n; ++i) s += b[i] * q[i];
            for (std::size_t i = 0; i < n; ++i) q[i] -= s * b[i];
         }
      Real norm = 0;
      for (const Real& v : q) norm += v * v;
      norm = sqrt(norm);
      for (Real& v : q) v /= norm;
      basis.push_back(q);
      Real s = 0;
      for (std::size_t i = 0; i < n; ++i) s += q[i] * y[i];
      for (std::size_t i = 0; i < n; ++i) p[i] += s * q[i];
   }
   return p;
}

// lambda_max = ||z||_inf for D^T z = r, r = y minus its least-squares polynomial of degree k: r lies in the
// range of D^T, so forward substitution solves it exactly (well conditioned, unlike D D^T z = D r). For
// lambda >= lambda_max the trend filter is that polynomial.
template <class Real>
Real trend_filter_lambda_max(const std::vector<Real>& r, int k)
{
   using std::abs;
   const int R = k + 1;
   Real d[4];
   long b = 1;
   for (int j = 0; j <= R; ++j) { d[j] = Real((R - j) % 2 ? -b : b); b = b * (R - j) / (j + 1); }
   std::vector<Real> z(r.size() - std::size_t(R));
   Real mx = 0;
   for (std::size_t i = 0; i < z.size(); ++i)
   {
      Real s = r[i];
      for (int j = 1; j <= R && std::size_t(j) <= i; ++j) s -= d[j] * z[i - std::size_t(j)];
      z[i] = s / d[0];
      mx = (std::max)(mx, Real(abs(z[i])));
   }
   return mx;
}

// The fit at fixed k in {0, 1, 2} and lambda >= 0, for y with at least k + 2 points. The degree-k
// polynomial trend is removed first, so the stopping tolerance is relative to what is left.
// ok: set to false if the interior-point solver failed (the result is then the error handler's).
template <class Real, class Policy = policies::policy<>>
std::vector<Real> trend_filter_fit(const std::vector<Real>& y, int k, const Real& lambda, const Policy& pol = Policy(), bool* ok_out = nullptr)
{
   if (lambda == 0)
      return y;
   const std::vector<Real> p = polynomial_fit(y, k);
   std::vector<Real> r(y.size()), x(y.size());
   Real scale = 0;
   for (std::size_t i = 0; i < y.size(); ++i) { r[i] = y[i] - p[i]; scale += r[i] * r[i]; }
   if (scale == 0)
      return y; // y is a polynomial of degree <= k
   if (k > 0 && lambda >= trend_filter_lambda_max(r, k))
      return p;
   bool ok = true;
   if (k == 0)
      trend_filter_tv(r.data(), x.data(), r.size(), lambda);
   else if (k == 1)
      ok = trend_filter_pdip<Real, 2>(r, lambda, scale, x, pol);
   else
      ok = trend_filter_pdip<Real, 3>(r, lambda, scale, x, pol);
   if (ok)
      for (std::size_t i = 0; i < y.size(); ++i) x[i] += p[i];
   if (ok_out) *ok_out = ok;
   return x;
}

// Two-fold CV (Nason 1996): the even samples y_0, y_2, ... (ceil(n/2) of them) and the odd ones,
// y_1, y_3, ...; each odd sample is predicted by the mean of the even fit at its neighbours and
// vice versa, or by the one neighbour at an end.
template <class Real>
Real cross_validation_score(const std::vector<Real>& xe, const std::vector<Real>& xo, const std::vector<Real>& ye, const std::vector<Real>& yo)
{
   const std::size_t he = ye.size(), ho = yo.size();
   Real s = 0;
   for (std::size_t i = 0; i < ho; ++i)
   {
      const Real po = i + 1 < he ? (xe[i] + xe[i + 1]) / 2 : xe[i];
      s += (po - yo[i]) * (po - yo[i]);
   }
   for (std::size_t i = 0; i < he; ++i)
   {
      const Real pe = i == 0 ? xo[0] : i < ho ? (xo[i - 1] + xo[i]) / 2 : xo[i - 1];
      s += (pe - ye[i]) * (pe - ye[i]);
   }
   return s;
}

// Robust noise level: the median absolute deviation (about the median) of the second differences,
// divided by sqrt(6) Phi^{-1}(3/4). Zero means the second differences are mostly constant: no noise.
template <class Real>
Real trend_filter_sigma(const std::vector<Real>& y)
{
   using std::abs; using std::sqrt;
   std::vector<Real> d(y.size() - 2);
   for (std::size_t i = 0; i + 2 < y.size(); ++i) d[i] = y[i] - 2 * y[i + 1] + y[i + 2];
   auto mid = d.begin() + d.size() / 2;
   std::nth_element(d.begin(), mid, d.end());
   const Real median = *mid;
   for (Real& v : d) v = abs(v - median);
   std::nth_element(d.begin(), mid, d.end());
   return *mid / (sqrt(Real(6)) * Real(0.6744897501960817));
}

// Even/odd cross-validation over k = 0, 1, 2 and lambda = sigma 10^L (sigma > 0). Per k: 6 values
// of L equally spaced over [lo_k, hi_k], one step of at most one decade beyond an end if the best is
// there, then two bisections around the best. The half samples use lambda / 2^(k+1). Above the lambda at which
// both halves' fits are polynomials every L gives the same fit, so L is capped there. best_score: the best CV score.
template <class Real, class Policy = policies::policy<>>
std::pair<int, Real> trend_filter_select(const std::vector<Real>& y, const Real& sigma, Real* best_score = nullptr, const Policy& pol = Policy())
{
   using std::pow; using std::log10;
   const std::size_t he = (y.size() + 1) / 2, ho = y.size() / 2;
   std::vector<Real> ye(he), yo(ho);
   for (std::size_t i = 0; i < y.size(); ++i) (i % 2 ? yo[i / 2] : ye[i / 2]) = y[i];
   std::pair<int, Real> best(0, Real(0));
   Real best_cv = (std::numeric_limits<Real>::max)();
   for (int k = 0; k <= 2 && ho >= std::size_t(k + 2); ++k)
   {
      const Real lo = Real(-2 + k) / 2, hi = Real(4 + k), coarse = (hi - lo) / 5;
      // the half-sample lambda above which both halves' fits are polynomials (k = 0: none needed)
      Real lmax = (std::numeric_limits<Real>::max)();
      if (k > 0)
      {
         lmax = 0;
         for (const std::vector<Real>* h : {&ye, &yo})
         {
            const std::vector<Real> p = polynomial_fit(*h, k);
            std::vector<Real> r(h->size());
            for (std::size_t i = 0; i < r.size(); ++i) r[i] = (*h)[i] - p[i];
            lmax = (std::max)(lmax, trend_filter_lambda_max(r, k));
         }
      }
      const Real Lcap = k > 0 && lmax > 0 ? Real(log10(lmax * Real(2 << k) / sigma)) : (std::numeric_limits<Real>::max)();
      auto score = [&](const Real& L) {
         const Real lh = (std::min)(Real(sigma * pow(Real(10), L) / Real(2 << k)), lmax);
         return cross_validation_score(trend_filter_fit(ye, k, lh, pol), trend_filter_fit(yo, k, lh, pol), ye, yo);
      };
      int bg = 0;
      Real bL = hi, bc = (std::numeric_limits<Real>::max)();
      for (int g = 0; g < 6; ++g)
      {
         const Real L = hi - g * coarse, c = score(L);
         if (c < bc) { bc = c; bL = L; bg = g; }
      }
      if ((bg == 0 && hi < Lcap) || bg == 5)
      {
         const Real L = bg == 0 ? hi + (std::min)(coarse, Real(1)) : lo - (std::min)(coarse, Real(1)), c = score(L);
         if (c < bc) { bc = c; bL = L; }
      }
      Real st = coarse;
      for (int r = 0; r < 2; ++r)
      {
         st /= 2;
         const Real c0 = bL;
         for (int sign = -1; sign <= 1; sign += 2)
         {
            const Real L = c0 + sign * st, c = score(L);
            if (c < bc) { bc = c; bL = L; }
         }
      }
      if (bc < best_cv)
      {
         best_cv = bc;
         best = std::make_pair(k, sigma * pow(Real(10), bL));
      }
   }
   if (best_score) *best_score = best_cv;
   return best;
}

// float is computed in double: the interior-point method needs more than 24 bits at large lambda
template <class Real>
using trend_filter_compute_type = typename std::conditional<std::is_same<Real, float>::value, double, Real>::type;

template <class RandomAccessContainer, class T>
RandomAccessContainer trend_filter_output(const RandomAccessContainer& v, const std::vector<T>& x)
{
   using Real = typename RandomAccessContainer::value_type;
   RandomAccessContainer out = v;
   std::transform(x.begin(), x.end(), std::begin(out), [](const T& a) { return static_cast<Real>(a); });
   return out;
}

template <class RandomAccessContainer>
RandomAccessContainer trend_filter_error(const RandomAccessContainer& v, const char* message, const typename RandomAccessContainer::value_type& value)
{
   RandomAccessContainer out = v;
   std::fill(std::begin(out), std::end(out), policies::raise_domain_error("boost::math::filters::trend_filter<%1%>", message, value, policies::policy<>()));
   return out;
}

} // namespace detail

// The trend filter of order k in {0, 1, 2} at lambda >= 0:
//    argmin_x 1/2 sum (v_i - x_i)^2 + lambda sum |(D^(k+1) x)_i|.
template <class RandomAccessContainer>
RandomAccessContainer trend_filter(const RandomAccessContainer& v, int k, const typename RandomAccessContainer::value_type& lambda)
{
   using Real = typename RandomAccessContainer::value_type;
   using T = detail::trend_filter_compute_type<Real>;
   const std::size_t n = std::size_t(std::distance(std::begin(v), std::end(v)));
   if (k < 0 || k > 2)
      return detail::trend_filter_error(v, "The order k must be 0, 1 or 2, but got %1%.", Real(k));
   if (!(lambda >= 0) || !(lambda <= (std::numeric_limits<Real>::max)()))
      return detail::trend_filter_error(v, "lambda must be finite and >= 0, but got %1%.", lambda);
   if (n < std::size_t(k + 2))
      return detail::trend_filter_error(v, "Order k needs at least k + 2 samples, but got %1%.", Real(n));
   const std::vector<T> y(std::begin(v), std::end(v));
   return detail::trend_filter_output(v, detail::trend_filter_fit(y, k, T(lambda)));
}

// The trend filter with k in {0, 1, 2} and lambda chosen by even/odd cross-validation; needs at
// least 4 samples. Assumes independent noise. Data without noise (sigma estimate 0) is returned unchanged.
template <class RandomAccessContainer>
RandomAccessContainer trend_filter(const RandomAccessContainer& v)
{
   using Real = typename RandomAccessContainer::value_type;
   using T = detail::trend_filter_compute_type<Real>;
   const std::vector<T> y(std::begin(v), std::end(v));
   if (y.size() < 4)
      return detail::trend_filter_error(v, "Cross-validation needs at least 4 samples, but got %1%.", Real(y.size()));
   const T sigma = detail::trend_filter_sigma(y);
   if (sigma == 0)
      return v;
   const std::pair<int, T> kl = detail::trend_filter_select(y, sigma);
   return detail::trend_filter_output(v, detail::trend_filter_fit(y, kl.first, kl.second));
}

} BOOST_MATH_NAMESPACE_END

#endif
