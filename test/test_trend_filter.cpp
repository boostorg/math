// Copyright Nick Thompson, 2026
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

#include "math_unit_test.hpp"
#include <array>
#include <cmath>
#include <limits>
#include <random>
#include <stdexcept>
#include <type_traits>
#include <vector>
#include <boost/math/filters/trend_filter.hpp>
#if !defined(BOOST_MATH_STANDALONE)
#include <boost/multiprecision/cpp_bin_float.hpp>
#endif

using boost::math::filters::trend_filter;

// sqrt(eps) of the precision the solver computes in (float is computed in double), and at least
// a few ulps of Real: the scale of the test tolerances
template <class Real>
Real tol_scale()
{
   using std::sqrt;
   using T = typename std::conditional<std::is_same<Real, float>::value, double, Real>::type;
   const Real s = Real(sqrt(std::numeric_limits<T>::epsilon()));
   return (std::max)(s, 16 * std::numeric_limits<Real>::epsilon());
}

// (k+1)-th difference coefficients, d_j = (-1)^(k+1-j) binomial(k+1, j)
template <class Real>
std::vector<Real> diff(int k)
{
   std::vector<Real> d(k + 2);
   long b = 1;
   for (int j = 0; j <= k + 1; ++j) { d[j] = Real((k + 1 - j) % 2 ? -b : b); b = b * (k + 1 - j) / (j + 1); }
   return d;
}

template <class Real>
std::vector<Real> apply_D(const std::vector<Real>& x, int k)
{
   const auto d = diff<Real>(k);
   std::vector<Real> out(x.size() - k - 1);
   for (std::size_t i = 0; i < out.size(); ++i)
      for (int j = 0; j <= k + 1; ++j) out[i] += d[j] * x[i + j];
   return out;
}

// Least-squares fit by polynomials of degree <= k (modified Gram-Schmidt, twice).
template <class Real>
std::vector<Real> polynomial_fit(const std::vector<Real>& y, int k)
{
   using std::sqrt;
   const std::size_t n = y.size();
   std::vector<std::vector<Real>> Q;
   for (int j = 0; j <= k; ++j)
   {
      std::vector<Real> q(n);
      for (std::size_t i = 0; i < n; ++i) { Real t = Real(2 * i) / Real(n - 1) - 1; q[i] = 1; for (int p = 0; p < j; ++p) q[i] *= t; }
      for (int pass = 0; pass < 2; ++pass)
         for (const auto& p : Q)
         {
            Real s = 0;
            for (std::size_t i = 0; i < n; ++i) s += p[i] * q[i];
            for (std::size_t i = 0; i < n; ++i) q[i] -= s * p[i];
         }
      Real nn = 0;
      for (const Real& v : q) nn += v * v;
      nn = sqrt(nn);
      for (Real& v : q) v /= nn;
      Q.push_back(q);
   }
   std::vector<Real> p(n);
   for (const auto& q : Q)
   {
      Real s = 0;
      for (std::size_t i = 0; i < n; ++i) s += q[i] * y[i];
      for (std::size_t i = 0; i < n; ++i) p[i] += s * q[i];
   }
   return p;
}

// (KKT checks are not run for float: rounding x to float, amplified by the forward substitution below, which grows
// like n^(k+1), swamps the check; float is computed in double, which is checked.)
// Dual certificate: z with D^T z = y - x by forward substitution; returns lambda_max for x = the
// polynomial fit, or checks the KKT conditions |z_i| <= lambda and z_i = lambda sign((D x)_i) at the knots.
template <class Real>
std::vector<Real> dual_from_residual(const std::vector<Real>& y, const std::vector<Real>& x, int k)
{
   const auto d = diff<Real>(k);
   const std::size_t n = y.size(), m = n - k - 1;
   std::vector<Real> z(m);
   for (std::size_t i = 0; i < m; ++i)
   {
      Real s = y[i] - x[i];
      for (std::size_t j = 1; j <= (std::min)(i, std::size_t(k + 1)); ++j) s -= d[j] * z[i - j];
      z[i] = s / d[0];
   }
   return z;
}

template <class Real>
std::vector<Real> test_signal(std::size_t n, unsigned seed, Real noise)
{
   using std::sin;
   std::mt19937_64 gen(seed);
   std::uniform_real_distribution<double> u(-1, 1);
   std::vector<Real> y(n);
   for (std::size_t i = 0; i < n; ++i)
   {
      const Real t = Real(i) / Real(n);
      y[i] = (t < Real(0.3) ? Real(1) : t < Real(0.6) ? Real(-2) : Real(0.5)) + 3 * t * t + sin(9 * t) + noise * Real(u(gen));
   }
   return y;
}

template <class Real>
void test_exact_cases()
{
   using std::abs;
   const auto y = test_signal<Real>(64, 1, Real(0.5));
   for (int k = 0; k <= 2; ++k)
   {
      // lambda = 0: the identity, exactly
      const auto x0 = trend_filter(y, k, Real(0));
      for (std::size_t i = 0; i < y.size(); ++i) CHECK_EQUAL(x0[i], y[i]);
      // lambda >= lambda_max = |z|_inf for the polynomial fit: the fit is the polynomial fit
      const auto p = polynomial_fit(y, k);
      Real lmax = 0;
      for (const Real& v : dual_from_residual(y, p, k)) lmax = (std::max)(lmax, abs(v));
      const auto x = trend_filter(y, k, 2 * lmax);
      Real err = 0, scale = 0;
      for (std::size_t i = 0; i < y.size(); ++i) { err = (std::max)(err, abs(x[i] - p[i])); scale = (std::max)(scale, abs(p[i])); }
      CHECK_LE(err, 1000 * tol_scale<Real>() * scale);
   }
}

// n = k + 2: D is one row d, and x = y - d z with z = clamp(d.y / d.d, -lambda, lambda).
template <class Real>
void test_closed_form()
{
   using std::abs;
   const Real tol = 100 * tol_scale<Real>();
   const std::vector<std::vector<Real>> data{{Real(0), Real(1)}, {Real(0), Real(3), Real(-1)}, {Real(1), Real(-2), Real(4), Real(0.5)}};
   for (int k = 0; k <= 2; ++k)
   {
      const std::vector<Real>& y = data[k];
      const auto d = diff<Real>(k);
      Real dy = 0, dd = 0, ymax = 0;
      for (int j = 0; j <= k + 1; ++j) { dy += d[j] * y[j]; dd += d[j] * d[j]; ymax = (std::max)(ymax, abs(y[j])); }
      // lambda below |d.y / d.d| saturates z; above it, the fit is the polynomial through the data's projection
      for (Real lambda : {abs(dy / dd) / 4, abs(dy / dd) * 4})
      {
         const Real z = (std::max)(-lambda, (std::min)(lambda, dy / dd));
         const auto x = trend_filter(y, k, lambda);
         for (int j = 0; j <= k + 1; ++j) CHECK_LE(abs(x[j] - (y[j] - d[j] * z)), tol * ymax);
      }
   }
}

// A polynomial trend of degree <= k far larger than the noise must not hide the noise from the solver.
template <class Real>
void test_polynomial_trend()
{
   using std::abs;
   const std::vector<Real> y0{0, 1, 0, 1};
   std::vector<Real> y = y0;
   for (std::size_t i = 0; i < y.size(); ++i) y[i] += Real(1000000) * Real(i);
   const auto x0 = trend_filter(y0, 1, Real(1)), x = trend_filter(y, 1, Real(1));
   Real moved = 0;
   for (std::size_t i = 0; i < y.size(); ++i)
   {
      CHECK_LE(abs((x[i] - Real(1000000) * Real(i)) - x0[i]), Real(1e-6));
      moved = (std::max)(moved, abs(x0[i] - y0[i]));
   }
   CHECK_GE(moved, Real(0.1)); // it is filtered, not returned as is
}

// Odd n: the even half has the extra sample, the last sample counts, and reversing the data
// reverses the fit (the even samples of the reversed data are the reversed even samples).
template <class Real>
void test_odd_length()
{
   using std::abs;
   for (std::size_t n : {std::size_t(5), std::size_t(7), std::size_t(301)})
   {
      const auto y = test_signal<Real>(n, 6, Real(0.4));
      const std::vector<Real> r(y.rbegin(), y.rend());
      const auto x = trend_filter(y), xr = trend_filter(r);
      Real e = 0, m = 0;
      for (std::size_t i = 0; i < n; ++i) { e = (std::max)(e, abs(x[i] - xr[n - 1 - i])); m = (std::max)(m, abs(x[i])); }
      CHECK_LE(e, 1000 * tol_scale<Real>() * m);
   }
}

// No noise (piecewise linear: the second differences vanish but at the kinks): returned unchanged.
template <class Real>
void test_noise_free()
{
   std::vector<Real> y(100);
   for (std::size_t i = 0; i < y.size(); ++i) y[i] = i < 40 ? Real(i) : Real(80) - Real(i);
   const auto x = trend_filter(y);
   for (std::size_t i = 0; i < y.size(); ++i) CHECK_EQUAL(x[i], y[i]);
}

template <class Real>
void test_kkt()
{
   using std::abs;
   const Real tol = 100 * tol_scale<Real>();
   const auto y = test_signal<Real>(200, 3, Real(0.4));
   for (int k = 0; k <= 2; ++k)
      for (Real lambda : {Real(0.3), Real(3), Real(30)})
      {
         const auto x = trend_filter(y, k, lambda);
         const auto z = dual_from_residual(y, x, k);
         const auto Dx = apply_D(x, k);
         Real dmax = 0;
         for (const Real& v : Dx) dmax = (std::max)(dmax, abs(v));
         Real zmax = 0, knot_defect = 0;
         for (std::size_t i = 0; i < z.size(); ++i)
         {
            zmax = (std::max)(zmax, abs(z[i]));
            if (abs(Dx[i]) > 1000 * tol * dmax) knot_defect = (std::max)(knot_defect, abs(z[i] - (Dx[i] > 0 ? lambda : -lambda)));
         }
         // the box |z| <= lambda, and z = lambda sign(D x) where D x is clearly nonzero
         CHECK_LE(zmax, lambda * (1 + 100 * tol));
         CHECK_LE(knot_defect, lambda * 1000 * tol);
         // y - x is orthogonal to the polynomials of degree <= k (the null space of D)
         const auto p = polynomial_fit([&] { std::vector<Real> r(y.size()); for (std::size_t i = 0; i < r.size(); ++i) r[i] = y[i] - x[i]; return r; }(), k);
         Real pmax = 0, ymax = 0;
         for (std::size_t i = 0; i < y.size(); ++i) { pmax = (std::max)(pmax, abs(p[i])); ymax = (std::max)(ymax, abs(y[i])); }
         CHECK_LE(pmax, 100 * tol * ymax);
      }
}

template <class Real>
void test_equivariance()
{
   using std::abs;
   const Real tol = 1000 * tol_scale<Real>();
   const auto y = test_signal<Real>(150, 4, Real(0.4));
   for (int k = 0; k <= 2; ++k)
   {
      // adding a polynomial of degree <= k adds it to the fit; scaling y and lambda scales the fit
      std::vector<Real> yp = y, ys = y;
      for (std::size_t i = 0; i < y.size(); ++i)
      {
         const Real t = Real(i) / Real(y.size());
         Real q = 2;
         if (k >= 1) q += 3 * t;
         if (k >= 2) q -= 5 * t * t;
         yp[i] += q;
         ys[i] *= 7;
      }
      const Real lambda = Real(2);
      const auto x = trend_filter(y, k, lambda), xp = trend_filter(yp, k, lambda), xs = trend_filter(ys, k, 7 * lambda);
      Real ep = 0, es = 0, xm = 0;
      for (std::size_t i = 0; i < y.size(); ++i)
      {
         ep = (std::max)(ep, abs((xp[i] - (yp[i] - y[i])) - x[i]));
         es = (std::max)(es, abs(xs[i] / 7 - x[i]));
         xm = (std::max)(xm, abs(x[i]));
      }
      CHECK_LE(ep, tol * xm);
      CHECK_LE(es, tol * xm);
   }
   // the cross-validated filter is scale equivariant too
   const auto y2 = test_signal<Real>(256, 5, Real(0.4));
   std::vector<Real> y2s = y2;
   for (Real& v : y2s) v *= 3;
   const auto c = trend_filter(y2), cs = trend_filter(y2s);
   Real e = 0, m = 0;
   for (std::size_t i = 0; i < y2.size(); ++i) { e = (std::max)(e, abs(cs[i] / 3 - c[i])); m = (std::max)(m, abs(c[i])); }
   CHECK_LE(e, tol * m);
}

// Noisy Blocks and Doppler (Donoho & Johnstone), SNR 7: the CV fit must cut the noise MSE by a wide margin.
template <class Real>
void test_cross_validation()
{
   using std::sin; using std::sqrt;
   const std::size_t n = 2048;
   const Real pi = Real(3.141592653589793238462643383279502884L);
   const std::array<Real, 11> tj{Real(0.1), Real(0.13), Real(0.15), Real(0.23), Real(0.25), Real(0.4), Real(0.44), Real(0.65), Real(0.76), Real(0.78), Real(0.81)};
   const std::array<Real, 11> hj{4, -5, 3, -4, 5, Real(-4.2), Real(2.1), Real(4.3), Real(-3.1), Real(2.1), Real(-4.2)};
   for (int s = 0; s < 2; ++s)
   {
      std::vector<Real> f(n);
      for (std::size_t i = 0; i < n; ++i)
      {
         const Real t = Real(i) / Real(n);
         if (s == 0) { for (std::size_t j = 0; j < tj.size(); ++j) if (t >= tj[j]) f[i] += hj[j]; }
         else f[i] = sqrt(t * (1 - t)) * sin(2 * pi * Real(1.05) / (t + Real(0.05)));
      }
      Real mean = 0, var = 0;
      for (const Real& v : f) mean += v / Real(n);
      for (const Real& v : f) var += (v - mean) * (v - mean) / Real(n);
      const Real sigma = sqrt(var) / 7;
      std::mt19937_64 gen(17 + s);
      std::uniform_real_distribution<double> u(0, 1);
      std::vector<Real> y(n);
      for (std::size_t i = 0; i < n; ++i)
      {  // Box-Muller, so the noise does not depend on the standard library's normal_distribution
         const double a = u(gen), b = u(gen);
         y[i] = f[i] + sigma * Real(std::sqrt(-2 * std::log(1 - a)) * std::cos(2 * 3.141592653589793 * b));
      }
      const auto x = trend_filter(y);
      Real mse = 0;
      for (std::size_t i = 0; i < n; ++i) mse += (x[i] - f[i]) * (x[i] - f[i]) / Real(n);
      // noise MSE is sigma^2; the CV fit reaches well under a quarter of it on both signals
      CHECK_LE(mse, sigma * sigma / 4);
   }
}

template <class Real>
void test_errors()
{
   const std::vector<Real> y{1, 2, 3, 5};
   CHECK_THROW(trend_filter(y, 3, Real(1)), std::domain_error);
   CHECK_THROW(trend_filter(y, -1, Real(1)), std::domain_error);
   CHECK_THROW(trend_filter(y, 1, Real(-1)), std::domain_error);
   CHECK_THROW(trend_filter(y, 1, std::numeric_limits<Real>::quiet_NaN()), std::domain_error);
   CHECK_THROW(trend_filter(std::vector<Real>{1, 2, 3}, 2, Real(1)), std::domain_error);
   CHECK_THROW(trend_filter(std::vector<Real>{1, 2, 3}), std::domain_error);
   // k + 2 samples is enough
   CHECK_NO_THROW(trend_filter(y, 2, Real(1)));
}

int main()
{
   test_errors<double>();
   test_exact_cases<float>();
   test_exact_cases<double>();
   test_closed_form<float>();
   test_closed_form<double>();
   test_polynomial_trend<double>();
   test_odd_length<double>();
   test_noise_free<double>();
   test_kkt<double>();
   test_equivariance<float>();
   test_equivariance<double>();
   test_cross_validation<float>();
   test_cross_validation<double>();
#if !defined(BOOST_MATH_STANDALONE)
   using mp = boost::multiprecision::cpp_bin_float_50;
   test_exact_cases<mp>();
   test_closed_form<mp>();
   test_kkt<mp>();
   test_equivariance<mp>();
#endif
   return boost::math::test::report_errors();
}
