//  (C) Copyright Nick Thompson 2026.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_SPECIAL_SPHEROIDAL_HPP
#define BOOST_MATH_SPECIAL_SPHEROIDAL_HPP

// Prolate and oblate spheroidal wave functions: the eigenvalue λ_n^m(±c²),
// the angular function Ps_n^m(η, ±c²) in the Meixner–Schäfke normalization
// of DLMF 30.4, and the radial functions of the first and second kind
// normalized as in DLMF 30.11 (Flammer's R_mn^(1), R_mn^(2)).

#include <boost/math/tools/config.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/math/tools/promotion.hpp>
#include <boost/math/tools/precision.hpp>
#include <boost/math/policies/error_handling.hpp>
#include <boost/math/special_functions/math_fwd.hpp>
#include <boost/math/special_functions/fpclassify.hpp>
#ifndef BOOST_MATH_BUILD_MODULE
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <initializer_list>
#include <limits>
#include <type_traits>
#include <vector>
#endif

BOOST_MATH_NAMESPACE_BEGIN

namespace detail { namespace spheroidal {

// Recurrence for the Legendre-expansion coefficients d_k (A&S 21.7.3):
//   α_k d_{k+2} + (β_k - λ) d_k + γ_k d_{k-2} = 0,
// with s = c² for prolate and s = -c² for oblate functions.
template <class T>
T alpha(unsigned m, unsigned k, const T& s)
{
   const T mk = T(2 * m + k);
   return (mk + 2) * (mk + 1) * s / ((T(2 * m + 2 * k) + 3) * (T(2 * m + 2 * k) + 5));
}

template <class T>
T beta(unsigned m, unsigned k, const T& s)
{
   const T l = T(m + k);
   const T mm = T(m);
   return l * (l + 1) + (2 * l * (l + 1) - 2 * mm * mm - 1) * s / ((2 * l - 1) * (2 * l + 3));
}

template <class T>
T gamma(unsigned m, unsigned k, const T& s)
{
   const T kk = T(k);
   const T l = T(m + k);
   return kk * (kk - 1) * s / ((2 * l - 3) * (2 * l - 1));
}

// Expansion coefficients d_k, k = p, p + 2, ..., of Ps_n^m(η, s) in Ferrers functions P_{m+k}^m(η).
template <class T>
struct coefficients
{
   unsigned n;
   unsigned m;
   unsigned p;   // parity of n - m
   T s;
   T lambda;
   T lambda_lo = 0;   // λ - lambda, from a double-length residual (eigenvalue_low_part)
   std::vector<T> d;  // d[j] = d_{p+2j}
};

// Number of eigenvalues of the truncated, symmetrized tridiagonal block that are below x.
template <class T>
std::size_t sturm_count(unsigned m, unsigned p, std::size_t size, const T& s, const T& x)
{
   BOOST_MATH_STD_USING
   std::size_t count = 0;
   T q = beta(m, p, s) - x;
   for (std::size_t j = 0;; ++j)
   {
      if (q < 0)
      {
         ++count;
      }
      if (j + 1 == size)
      {
         break;
      }
      if (q == 0)
      {
         q = tools::epsilon<T>() * (abs(x) + 1);
      }
      const unsigned k = p + 2 * static_cast<unsigned>(j);
      q = beta(m, k + 2, s) - x - alpha(m, k, s) * gamma(m, k + 2, s) / q;
   }
   return count;
}

// Double-length arithmetic for the eigenvalue's low part. λ ≈ -c² is a difference of terms of size c², so
// rounding it costs an absolute error of about ε c²; where the solutions grow exponentially, which is where
// the integrations below run, that becomes a relative error of several ε in the result. As Arb does with its
// guard bits, the one step that amplifies the rounding is redone at higher precision: Bouwkamp's residual
// at the rounded λ, in pairs (hi, lo) with |lo| <= ulp(hi)/2, gives λ - λ_hi by one Newton step.
template <class T>
struct twofold
{
   T hi;
   T lo;
   twofold(const T& h = T(0), const T& l = T(0)) : hi(h), lo(l) {}
   static twofold sum(const T& a, const T& b)
   {
      const T s = a + b;
      const T bb = s - a;
      return { s, (a - (s - bb)) + (b - bb) };
   }
   friend twofold operator+(const twofold& x, const twofold& y)
   {
      const twofold s = sum(x.hi, y.hi);
      return sum(s.hi, T(s.lo + x.lo + y.lo));
   }
   friend twofold operator-(const twofold& x) { return { -x.hi, -x.lo }; }
   friend twofold operator-(const twofold& x, const twofold& y) { return x + -y; }
   friend twofold operator*(const twofold& x, const twofold& y)
   {
      BOOST_MATH_STD_USING
      const T p = fma(x.hi, y.hi, T(0)); // fma(a, b, 0): the rounded product, which can't be fused into the sum below
      return sum(p, T(fma(x.hi, y.hi, -p) + x.hi * y.lo + x.lo * y.hi));
   }
   friend twofold operator*(const twofold& x, const T& y)
   {
      BOOST_MATH_STD_USING
      const T p = fma(x.hi, y, T(0));
      return sum(p, T(fma(x.hi, y, -p) + x.lo * y));
   }
   friend twofold operator/(const twofold& x, const twofold& y)
   {
      const T q1 = x.hi / y.hi;
      const twofold r = x - y * twofold(q1);
      return sum(q1, T(r.hi / y.hi));
   }
};

template <class T>
twofold<T> twofold_sum(const T& a, const T& b)
{
   return twofold<T>::sum(a, b);
}

// The leading part of a number: itself, or the high half of a double-length one.
template <class T>
T leading(const T& x)
{
   return x;
}

template <class T>
T leading(const twofold<T>& x)
{
   return x.hi;
}

template <class T>
twofold<T> exact(unsigned long long v)
{
   const T hi = T(v);
   return { hi, T(static_cast<long long>(v - static_cast<unsigned long long>(hi))) };
}

// α, β and γ in double length: integer ratios times s.
template <class T>
void recurrence_twofold(unsigned m, unsigned k, const T& s, twofold<T>* a, twofold<T>* b, twofold<T>* g)
{
   const twofold<T> S = { s, T(0) };
   const unsigned long long mk = 2ull * m + k;
   const unsigned long long l = static_cast<unsigned long long>(m) + k;
   const unsigned long long d = 2 * l;
   *a = exact<T>((mk + 2) * (mk + 1)) * S / exact<T>((d + 3) * (d + 5));
   // 2l(l+1) - 2m² - 1 >= 2m(m+1) - 2m² - 1 >= -1; the denominator is (2l-1)(2l+3), negative only for l = 0.
   const long long num = static_cast<long long>(2 * l * (l + 1)) - 2 * static_cast<long long>(m) * m - 1;
   const long long den = (static_cast<long long>(d) - 1) * (static_cast<long long>(d) + 3);
   const twofold<T> ratio = (num < 0 ? -exact<T>(static_cast<unsigned long long>(-num)) : exact<T>(static_cast<unsigned long long>(num)))
                            / (den < 0 ? -exact<T>(static_cast<unsigned long long>(-den)) : exact<T>(static_cast<unsigned long long>(den)));
   *b = exact<T>(l * (l + 1)) + ratio * S;
   *g = k < 2 ? twofold<T>{ T(0), T(0) } : exact<T>(static_cast<unsigned long long>(k) * (k - 1)) * S / exact<T>((d - 3) * (d - 1));
}

// Bouwkamp's function: the residual of the recurrence at k = n - m, with the ratios d_{k±2}/d_k
// taken from the continued fractions that are stable away from k = n - m.
template <class T>
T bouwkamp(unsigned m, unsigned r, unsigned p, unsigned top, const T& s, const T& x, T* derivative)
{
   T rho = 0;
   T drho = 0;
   for (unsigned k = top; k > r; k -= 2)
   {
      const T a = alpha(m, k, s);
      const T denom = beta(m, k, s) - x + a * rho;
      const T ddenom = -1 + a * drho;
      const T g = gamma(m, k, s);
      rho = -g / denom;
      drho = g * ddenom / (denom * denom);
   }
   T sigma = 0;
   T dsigma = 0;
   for (unsigned k = p; k < r; k += 2)
   {
      const T g = gamma(m, k, s);
      const T denom = beta(m, k, s) - x + g * sigma;
      const T ddenom = -1 + g * dsigma;
      const T a = alpha(m, k, s);
      sigma = -a / denom;
      dsigma = a * ddenom / (denom * denom);
   }
   *derivative = -1 + alpha(m, r, s) * drho + gamma(m, r, s) * dsigma;
   return beta(m, r, s) - x + alpha(m, r, s) * rho + gamma(m, r, s) * sigma;
}

// λ - x for the eigenvalue x, rounded to T: Bouwkamp's residual at x in double length, divided by its slope.
template <class T>
T eigenvalue_low_part(unsigned m, unsigned r, unsigned p, unsigned top, const T& s, const T& x)
{
   BOOST_MATH_STD_USING
   const twofold<T> X = { x, T(0) };
   twofold<T> a;
   twofold<T> b;
   twofold<T> g;
   twofold<T> rho = { T(0), T(0) };
   for (unsigned k = top; k > r; k -= 2)
   {
      recurrence_twofold(m, k, s, &a, &b, &g);
      rho = -g / (b + -X + a * rho);
   }
   twofold<T> sigma = { T(0), T(0) };
   for (unsigned k = p; k < r; k += 2)
   {
      recurrence_twofold(m, k, s, &a, &b, &g);
      sigma = -a / (b + -X + g * sigma);
   }
   recurrence_twofold(m, r, s, &a, &b, &g);
   const twofold<T> f = b + -X + a * rho + g * sigma;
   T df;
   bouwkamp(m, r, p, top, s, x, &df);
   const T lo = -f.hi / df;
   // A converged eigenvalue is within a few ε of the scale |λ| + |s| of the recurrence; a larger step is no correction.
   return (BOOST_MATH_NAMESPACE::isfinite)(lo) && abs(lo) <= 64 * tools::epsilon<T>() * (abs(x) + abs(s)) ? lo : T(0);
}

template <class T>
unsigned initial_top(unsigned m, unsigned n, const T& s)
{
   BOOST_MATH_STD_USING
   const T c = sqrt(abs(s));
   const unsigned extra = static_cast<unsigned>(ceil(c)) + 2 * static_cast<unsigned>(tools::digits<T>() / 3) + 20;
   // Same parity as n - m.
   return (n - m) + 2 * extra;
}

template <class T, class Policy>
T eigenvalue(unsigned m, unsigned n, const T& s, unsigned top, const Policy& pol, const char* function)
{
   BOOST_MATH_STD_USING
   const unsigned r = n - m;
   const unsigned p = r & 1u;
   const std::size_t size = (top - p) / 2 + 1;
   const std::size_t index = (r - p) / 2;
   // Gershgorin bounds on the symmetrized tridiagonal matrix.
   T lo = 0;
   T hi = 0;
   for (std::size_t j = 0; j < size; ++j)
   {
      const unsigned k = p + 2 * static_cast<unsigned>(j);
      T radius = 0;
      if (j + 1 < size)
      {
         radius += sqrt(abs(alpha(m, k, s) * gamma(m, k + 2, s)));
      }
      if (j > 0)
      {
         radius += sqrt(abs(alpha(m, k - 2, s) * gamma(m, k, s)));
      }
      const T b = beta(m, k, s);
      if (j == 0 || b - radius < lo)
      {
         lo = b - radius;
      }
      if (j == 0 || b + radius > hi)
      {
         hi = b + radius;
      }
   }
   // Bisect to roughly single precision, then polish with Newton's method on Bouwkamp's function.
   const T coarse = sqrt(tools::epsilon<T>());
   for (int i = 0; i < 2 * tools::digits<T>() && hi - lo > coarse * (abs(lo) + abs(hi) + 1); ++i)
   {
      const T mid = (lo + hi) / 2;
      if (sturm_count(m, p, size, s, mid) > index)
      {
         hi = mid;
      }
      else
      {
         lo = mid;
      }
   }
   T x = (lo + hi) / 2;
   T last_step = tools::max_value<T>();
   for (unsigned i = 0; i < 64; ++i)
   {
      T df;
      const T f = bouwkamp(m, r, p, top, s, x, &df);
      T next = x - f / df;
      if (!(next >= lo && next <= hi))
      {
         // Newton left the bracket; fall back to a bisection step.
         if (sturm_count(m, p, size, s, x) > index)
         {
            hi = x;
         }
         else
         {
            lo = x;
         }
         next = (lo + hi) / 2;
      }
      const T step = abs(next - x);
      x = next;
      // Quadratic convergence stops once the steps are rounding noise.
      if (step <= 2 * tools::epsilon<T>() * (abs(x) + 1) || (i > 2 && step >= last_step / 2))
      {
         return x;
      }
      last_step = step;
   }
   return policies::raise_evaluation_error<T>(function, "Eigenvalue iteration did not converge; last value %1%", x, pol);
}

// One Taylor step of length h from a for (w, w') along
//   q² w'' + 2x q w' - [(λ - c2 x²) q + σ m²] w = 0,  q = x² - σ,
// leaving in u the scaled coefficients u_k = w^(k)(a) h^k / k!, so that w(a + t h) = Σ u_k t^k on [0, 1].
// N is T, or twofold<T> with λ and c2 passed in double length.
template <class T, class N>
void taylor_step(const N& lambda, unsigned m, int sigma, const N& c2, const T& a, const T& h, N* R, N* Rd, std::vector<N>& u, const N& lambda_lo = N(0))
{
   BOOST_MATH_STD_USING
   const N A = a;
   const N H = h;
   const N one = T(1);
   // Polynomial coefficients in t = ξ - a.
   const N q0 = sigma > 0 ? (A - one) * (A + one) : A * A + one;
   const N q[3] = { q0, N(2 * a), one };
   const N P2[5] = { q0 * q0, N(4 * a) * q0, N(4 * a) * A + N(T(2)) * q0, N(4 * a), one };
   const N P1[4] = { N(2 * a) * q0, N(T(2)) * q0 + N(4 * a) * A, N(6 * a), N(T(2)) };
   const N Lc[3] = { N(lambda - c2 * A * A) + lambda_lo, -(N(2 * a) * c2), -c2 };
   N P0[5] = { T(0), T(0), T(0), T(0), T(0) };
   for (int i = 0; i < 3; ++i)
   {
      for (int k = 0; k < 3; ++k)
      {
         P0[i + k] = P0[i + k] - Lc[i] * q[k];
      }
   }
   P0[0] = P0[0] - N(T(sigma) * T(m) * T(m));
   // u_k neither overflows nor underflows near the singular point; the recurrence's coefficients scale by powers of h.
   N S2[5];
   N S1[4];
   N S0[5];
   N hp = one;
   for (int i = 0; i < 5; ++i)
   {
      S2[i] = P2[i] * hp;
      if (i < 4)
      {
         S1[i] = P1[i] * hp * H;
      }
      S0[i] = P0[i] * hp * H * H;
      hp = hp * H;
   }
   const T eps = tools::epsilon<T>();
   const T tol = std::is_same<N, T>::value ? T(eps / 4) : T(eps * eps / 4);
   u.assign({ *R, H * *Rd });
   N value = u[0] + u[1];
   N dvalue = u[1];
   for (unsigned k = 0; k < 100000; ++k)
   {
      N acc = T(0);
      for (unsigned i = 0; i <= 4; ++i)
      {
         if (i > 0 && i <= k + 2)
         {
            acc = acc + S2[i] * u[k - i + 2] * T(T(k - i + 2) * T(k - i + 1));
         }
         if (i < 4 && i <= k + 1)
         {
            acc = acc + S1[i] * u[k - i + 1] * T(k - i + 1);
         }
         if (i <= k)
         {
            acc = acc + S0[i] * u[k - i];
         }
      }
      const N next = -acc / (S2[0] * T(T(k + 2) * T(k + 1)));
      u.push_back(next);
      value = value + next;
      dvalue = dvalue + next * T(k + 2);
      // Both sums must have converged: near a singular point the terms fall only by about half per order, and the
      // derivative's terms carry an extra factor k, so a test on the value alone leaves it tens of ε short.
      const T v = tol * abs(leading(value));
      const T dv = tol * abs(leading(dvalue));
      if (k > 8 && abs(leading(next)) <= v && abs(leading(u[k + 1])) <= v && T(k + 2) * abs(leading(next)) <= dv && T(k + 1) * abs(leading(u[k + 1])) <= dv)
      {
         break;
      }
   }
   *R = value;
   *Rd = dvalue / H;
}

// The next point from a toward b: inside half the radius of convergence and, so that the series does not
// cancel, within about one local wavelength; representable, so that the series is summed exactly there.
template <class T>
T taylor_next(int sigma, const T& wavelength, const T& a, const T& b)
{
   BOOST_MATH_STD_USING
   const T radius = sigma > 0 ? T(abs(abs(a) - 1)) : T(sqrt(a * a + 1));
   const T h = (std::min)(T(radius / 2), wavelength);
   return abs(b - a) <= h ? b : T(b < a ? a - h : a + h);
}

// Continue (w, w') from a to b by Taylor steps. With σ = 1 and c2 = ±c² this is the prolate radial equation
// (x = ξ > 1) and also the angular one (x = η, |η| < 1); with σ = -1 and c2 = c² it is the oblate radial equation.
template <class T, class N>
void taylor_continue(const N& lambda, unsigned m, int sigma, const N& c2, T a, const T& b, N* R, N* Rd, const N& lambda_lo = N(0))
{
   BOOST_MATH_STD_USING
   // In double length the series may grow by e^4 before it cancels.
   const T wavelength = (std::is_same<N, T>::value ? T(1) : T(4)) / sqrt((std::max)(T(abs(leading(c2))), T(abs(leading(lambda)))) + 1);
   std::vector<N> u;
   while (a != b)
   {
      const T next = taylor_next(sigma, wavelength, a, b);
      taylor_step(lambda, m, sigma, c2, a, T(next - a), R, Rd, u, lambda_lo);
      a = next;
   }
}

// Σ d_k P_{m+k}^m(η) and its derivative for Ferrers functions, written as A(η) Σ d_k g_{m+k}(η) with
// A = (-1)^m (2m-1)!! (1-η²)^{m/2} and g_l = P_l^m / A a polynomial obeying the Legendre recurrence.
// abs_sum and abs_dsum receive Σ |d_k P_{m+k}^m(η)| and the same for the derivative: the scales against
// which the sums may cancel.
template <class T, class Policy>
T legendre_series(unsigned m, unsigned p, const std::vector<T>& d, const T& eta, T* derivative, T* abs_sum, T* abs_dsum, const Policy& pol, const char* function)
{
   BOOST_MATH_STD_USING
   const T w = (1 - eta) * (1 + eta);
   const T root = sqrt(w);
   T sum = 0;
   T dsum = 0;
   T asum = 0;
   T adsum = 0;
   T g_prev = 0;
   T g = 1;
   T dg_prev = 0;
   T dg = 0;
   std::size_t j = 0;
   for (unsigned l = m; j < d.size(); ++l)
   {
      if (l - m == p + 2 * j)
      {
         sum += d[j] * g;
         dsum += d[j] * dg;
         asum += abs(d[j] * g);
         adsum += abs(d[j] * dg);
         ++j;
      }
      const T L = T(l);
      const T g_next = ((2 * L + 1) * eta * g - (L + m) * g_prev) / (L - m + 1);
      const T dg_next = ((2 * L + 1) * (g + eta * dg) - (L + m) * dg_prev) / (L - m + 1);
      g_prev = g;
      g = g_next;
      dg_prev = dg;
      dg = dg_next;
   }
   // A and A' = -m η (-1)^m (2m-1)!! (1-η²)^{m/2-1}
   T A = 1;
   T Ad = 0;
   if (m > 0)
   {
      T B = 1;  // (2m-1)!! (1-η²)^{(m-2)/2}, built without spurious overflow
      for (unsigned i = 1; i <= m; ++i)
      {
         A *= (2 * T(i) - 1) * root;
         B *= (2 * T(i) - 1);
         if (i + 2 <= m)
         {
            B *= root;
         }
      }
      if (m == 1 && derivative)
      {
         B = (w == 0) ? policies::raise_overflow_error<T>(function, nullptr, pol) : T(B / root);
      }
      if (m & 1u)
      {
         A = -A;
         B = -B;
      }
      Ad = -T(m) * eta * B;
   }
   if (derivative)
   {
      *derivative = A * dsum + (Ad == 0 ? T(0) : T(Ad * sum));
   }
   if (abs_sum)
   {
      *abs_sum = abs(A) * asum;
   }
   if (abs_dsum)
   {
      *abs_dsum = abs(A) * adsum + abs(Ad) * asum;
   }
   return A * sum;
}

// Ps(η)/Ps(0) (n - m even) or Ps(η)/Ps'(0) (n - m odd), and its derivative, from the angular equation
//   (1-η²)² w'' - 2η(1-η²) w' + [(λ - sη²)(1-η²) - m²] w = 0
// integrated outward from η = 0, where the solution is dominant for the functions whose Legendre
// series cancels there (oblate, large c).
template <class T>
T angular_from_origin(const T& lambda, unsigned m, unsigned p, const T& s, const T& eta, T* derivative, const T& lambda_lo = T(0))
{
   T w = p == 0 ? T(1) : T(0);
   T wd = p == 0 ? T(0) : T(1);
   taylor_continue(lambda, m, 1, s, T(0), eta, &w, &wd, lambda_lo);
   *derivative = wd;
   return w;
}

// Ratio [(2m+k)!/k!] / [(2m+r)!/r!] for k ≡ r (mod 2), k = p + 2j.
template <class T>
std::vector<T> factorial_ratios(unsigned m, unsigned r, unsigned p, std::size_t size)
{
   std::vector<T> F(size);
   const std::size_t i0 = (r - p) / 2;
   F[i0] = 1;
   for (std::size_t j = i0; j + 1 < size; ++j)
   {
      const unsigned k = p + 2 * static_cast<unsigned>(j);
      F[j + 1] = F[j] * (T(2 * m + k) + 1) * (T(2 * m + k) + 2) / ((T(k) + 1) * (T(k) + 2));
   }
   for (std::size_t j = i0; j > 0; --j)
   {
      const unsigned k = p + 2 * static_cast<unsigned>(j) - 2;
      F[j - 1] = F[j] * ((T(k) + 1) * (T(k) + 2)) / ((T(2 * m + k) + 1) * (T(2 * m + k) + 2));
   }
   return F;
}

template <class T, class Policy>
coefficients<T> compute_coefficients(unsigned n, unsigned m, const T& s, const Policy& pol, const char* function, bool normalize = true)
{
   BOOST_MATH_STD_USING
   coefficients<T> result;
   result.n = n;
   result.m = m;
   result.s = s;
   const unsigned r = n - m;
   const unsigned p = r & 1u;
   result.p = p;
   unsigned top = initial_top(m, n, s);
   for (;;)
   {
      result.lambda = eigenvalue(m, n, s, top, pol, function);
      const T& x = result.lambda;
      const std::size_t size = (top - p) / 2 + 1;
      std::vector<T>& d = result.d;
      // Build d from d[match] = 1 with ratios that are stable in their direction: the continued fraction
      // from the top for j > match and from the bottom for j < match. Matching at n - m and then again at
      // the largest |d_k| keeps both recurrences away from the side where the solution is subdominant.
      // For large c the denominators β_k - λ + ... are small differences of terms of size c², so once the
      // eigenvalue's low part is known the ratios are recomputed in double length (with λ = x + λ_lo), which
      // keeps the low-order d_k to full relative accuracy.
      bool precise = false;
      auto ratio_step = [&](unsigned k, const twofold<T>& previous, bool upward) {
         if (!precise)
         {
            return twofold<T>{ upward ? T(-gamma(m, k, s) / (beta(m, k, s) - x + alpha(m, k, s) * previous.hi))
                                      : T(-alpha(m, k, s) / (beta(m, k, s) - x + gamma(m, k, s) * previous.hi)), T(0) };
         }
         twofold<T> a;
         twofold<T> b;
         twofold<T> g;
         recurrence_twofold(m, k, s, &a, &b, &g);
         const twofold<T> denom = b + twofold<T>{ T(-x), T(-result.lambda_lo) } + (upward ? a : g) * previous;
         return -(upward ? g : a) / denom;
      };
      auto build = [&](std::size_t match) {
         d.assign(size, T(0));
         d[match] = 1;
         twofold<T> ratio = { T(0), T(0) };
         std::vector<T> rho(size, T(0));
         for (std::size_t j = size - 1; j > match; --j)
         {
            const unsigned k = p + 2 * static_cast<unsigned>(j);
            ratio = ratio_step(k, ratio, true);
            rho[j] = ratio.hi + ratio.lo;  // d_k / d_{k-2}
         }
         for (std::size_t j = match + 1; j < size; ++j)
         {
            d[j] = d[j - 1] * rho[j];
         }
         ratio = { T(0), T(0) };
         for (std::size_t j = 0; j < match; ++j)
         {
            const unsigned k = p + 2 * static_cast<unsigned>(j);
            ratio = ratio_step(k, ratio, false);
            rho[j] = ratio.hi + ratio.lo;  // d_k / d_{k+2}
         }
         for (std::size_t j = match; j > 0; --j)
         {
            d[j - 1] = d[j] * rho[j - 1];
         }
      };
      build((r - p) / 2);
      std::size_t peak = 0;
      for (std::size_t j = 1; j < size; ++j)
      {
         if (abs(d[j]) > abs(d[peak]))
         {
            peak = j;
         }
      }
      std::size_t match = (r - p) / 2;
      if (peak != match)
      {
         match = peak;
         build(match);
      }
      T dmax = 0;
      for (const T& v : d)
      {
         dmax = (std::max)(dmax, T(abs(v)));
      }
      // The tail must be negligible, otherwise the truncation was too aggressive.
      if (abs(d.back()) <= tools::epsilon<T>() * tools::epsilon<T>() * dmax || top > 100000)
      {
         if (top > 100000)
         {
            policies::raise_evaluation_error<T>(function, "Expansion coefficients failed to converge for c^2 = %1%", s, pol);
         }
         result.lambda_lo = eigenvalue_low_part(m, r, p, top, s, result.lambda);
         if (normalize)
         {
            precise = true;
            build(match);
         }
         break;
      }
      top = 2 * top + p;
   }
   if (!normalize)
   {
      return result;
   }
   std::vector<T>& d = result.d;
   const std::size_t size = d.size();
   // Meixner–Schäfke normalization (DLMF 30.4.4): ||Ps_n^m|| = ||P_n^m||.
   const std::vector<T> F = factorial_ratios<T>(m, r, p, size);
   T norm = 0;
   for (std::size_t j = 0; j < size; ++j)
   {
      const unsigned k = p + 2 * static_cast<unsigned>(j);
      norm += d[j] * d[j] * F[j] * (T(2 * n) + 1) / (T(2 * m + 2 * k) + 1);
   }
   // Sign (DLMF 30.4.1): Ps_n^m(0) (n - m even) or Ps_n^m'(0) (n - m odd) has the sign of P_n^m(0) or
   // P_n^m'(0). The Legendre series at η = 0 can cancel almost completely (oblate, large c); then the
   // sign comes from matching the series where it does not cancel to a solution integrated from η = 0.
   T scale = 1 / sqrt(norm);
   const T target = (((n + m - p) / 2) & 1u) ? T(-1) : T(1);  // sign of P_n^m(0) or P_n^m'(0)
   T derivative_at_zero;
   T abs_value;
   T abs_derivative;
   const T v0 = legendre_series(m, p, d, T(0), &derivative_at_zero, &abs_value, &abs_derivative, pol, function);
   const T value_at_zero = p == 0 ? v0 : derivative_at_zero;
   const T abs_at_zero = p == 0 ? abs_value : abs_derivative;
   T sgn;
   if (abs(value_at_zero) > 1024 * tools::epsilon<T>() * abs_at_zero)
   {
      sgn = value_at_zero * target > 0 ? T(1) : T(-1);
   }
   else
   {
      sgn = 0;
      for (T eta = T(0.5); eta < 1; eta = (1 + eta) / 2)
      {
         T S_d;
         T S_abs;
         const T S = legendre_series(m, p, d, eta, &S_d, &S_abs, static_cast<T*>(nullptr), pol, function);
         if (abs(S) < 1024 * tools::epsilon<T>() * S_abs)
         {
            continue;
         }
         T u_d;
         const T u = angular_from_origin(result.lambda, m, p, s, eta, &u_d);
         // S = A u with A the value (or slope) at 0, so sign(A) = sign(S u + S' u').
         const T A = S * u + S_d * u_d;
         sgn = A * target > 0 ? T(1) : T(-1);
         break;
      }
      if (sgn == 0)
      {
         policies::raise_evaluation_error<T>(function, "Could not determine the sign of the angular function for c^2 = %1%", s, pol);
         sgn = 1;
      }
   }
   scale *= sgn;
   for (T& v : d)
   {
      v *= scale;
   }
   return result;
}

// The regular solution (1-η²)^{m/2} g(η) of the angular equation near η = 1, with g(1) = 1 and g analytic:
// in t = η - 1 the equation (1-η²) g'' - 2(m+1)η g' + (λ - sη² - m(m+1)) g = 0 gives
//   2(k+1)(k+m+1) e_{k+1} = (λ - m(m+1) - s - k(k-1) - 2(m+1)k) e_k - 2s e_{k-1} - s e_{k-2}.
// Returns the solution and its derivative at η = 1 + t, t <= 0.
template <class T>
T angular_near_one(const T& lambda, unsigned m, const T& s, const T& t, T* derivative, const T& lambda_lo = T(0))
{
   BOOST_MATH_STD_USING
   const T mm = T(m);
   const T shift = T(lambda - mm * (mm + 1) - s) + lambda_lo;
   if (t == 0)
   {
      // Ps(1) for m = 0 (the only case with a nonzero value there), with g'(1) = e_1.
      *derivative = shift / (2 * (mm + 1));
      return 1;
   }
   T e_2 = 0;
   T e_1 = 0;
   T e = 1;
   T g = 1;
   T gd = 0;
   T tk = 1;  // t^k
   for (unsigned k = 0; k < 100000; ++k)
   {
      const T K = T(k);
      const T next = ((shift - K * (K - 1) - 2 * (mm + 1) * K) * e - 2 * s * e_1 - s * e_2) / (2 * (K + 1) * (K + mm + 1));
      gd += (K + 1) * next * tk;
      tk *= t;
      const T term = next * tk;
      g += term;
      e_2 = e_1;
      e_1 = e;
      e = next;
      if (k > 4 && abs(term) <= tools::epsilon<T>() * abs(g) / 4 && abs(e_1 * tk / t) <= tools::epsilon<T>() * abs(g))
      {
         break;
      }
   }
   const T eta = 1 + t;
   const T w = -t * (2 + t);  // 1 - η²
   const T A = pow(w, mm / 2);
   *derivative = A * gd - (m == 0 ? T(0) : T(mm * eta * g * pow(w, mm / 2 - 1)));
   return A * g;
}

// Ps_n^m(η, s) and its derivative. Where the Legendre series cancels, which happens near η = 0 for
// oblate and near η = ±1 for prolate functions with large c, the function is a multiple of a solution
// that is dominant in the direction of integration: started at η = 0 (by parity) or at η = ±1 (the
// regular solution), integrated to η and on to a point where the series is accurate, and scaled to match.
template <class T, class Policy>
T angular(const coefficients<T>& co, const T& eta, T* derivative, const Policy& pol, const char* function)
{
   BOOST_MATH_STD_USING
   if (abs(eta) == 1 && co.m > 0)
   {
      // Ps vanishes at the endpoints for m > 0, where the series is exact; its derivative is infinite for m = 1.
      return legendre_series(co.m, co.p, co.d, eta, derivative, static_cast<T*>(nullptr), static_cast<T*>(nullptr), pol, function);
   }
   // The series is accurate where it cancels by no more than the function's local envelope allows:
   // max(|Ps|, |Ps'|/k), k the local wavenumber of the angular equation. Near a zero of an oscillating Ps
   // that tolerates the cancellation its conditioning permits; where Ps grows exponentially it does not,
   // since there the cancellation would amplify the errors of the coefficients d_k.
   const T mm = T(co.m);
   auto accurate = [&](const T& x, const T& S, const T& Sd, const T& a, const T& ad) {
      const T w = (1 - x) * (1 + x);
      const T V = co.lambda - co.s * x * x - mm * mm / w;
      const T k = sqrt(abs(V) / w) + 1;
      const T Sdd = (2 * x * Sd - V * S) / w;
      // Where the equation oscillates (V > 0) cancellation near a zero is within the conditioning |x S'|;
      // where it grows exponentially the cancellation would amplify the errors of the d_k.
      const T cond = V > 0 ? T(abs(x * Sd)) : T(0);
      const T dcond = V > 0 ? T(abs(x * Sdd)) : T(0);
      return a <= 2 * (std::max)({ T(abs(S)), T(abs(Sd) / k), cond }) && ad <= 2 * (std::max)({ T(abs(Sd)), T(abs(Sdd) / k), dcond });
   };
   T d;
   T a;
   T ad;
   T value = legendre_series(co.m, co.p, co.d, eta, &d, &a, &ad, pol, function);
   if (accurate(eta, value, d, a, ad))
   {
      if (derivative)
      {
         *derivative = d;
      }
      return value;
   }
   // Work at |η| and restore the parity Ps(-η) = (-1)^{n-m} Ps(η) at the end.
   const T x = abs(eta);
   const bool flip = eta < 0;
   const T parity = co.p ? T(-1) : T(1);
   // Integrate in the direction in which the function grows: inward from η = 1 for prolate functions
   // (outward from 0 their rounding errors would excite the growing solution), outward from 0 for oblate ones.
   const bool from_end = co.s > 0;
   T u;
   T ud;
   T xs = x;
   T us;
   T usd;
   if (from_end)
   {
      // The power series at η = 1 up to 1 - h0, then integration.
      const T h0 = (std::min)(T(0.5), T(1 / (2 * sqrt((std::max)(co.s, T(abs(co.lambda)))))));
      us = angular_near_one(co.lambda, co.m, co.s, T(-h0), &usd, co.lambda_lo);
      xs = 1 - h0;
      if (x < xs)
      {
         taylor_continue(co.lambda, co.m, 1, co.s, xs, x, &us, &usd, co.lambda_lo);
         xs = x;
         u = us;
         ud = usd;
      }
      else
      {
         u = angular_near_one(co.lambda, co.m, co.s, T(x - 1), &ud, co.lambda_lo);
      }
   }
   else
   {
      us = u = angular_from_origin(co.lambda, co.m, co.p, co.s, x, &ud, co.lambda_lo);
      usd = ud;
   }
   // Move toward the region where the series is accurate and match there.
   bool matched = false;
   for (unsigned i = 0; i < 256 && !matched; ++i)
   {
      // Toward where the series converges without cancellation: η = 0 for prolate, η = 1 for oblate functions.
      const T next = co.s > 0 && xs > 0 ? T(xs * 3 / 4) : T(xs + (1 - xs) / 4);
      taylor_continue(co.lambda, co.m, 1, co.s, xs, next, &us, &usd, co.lambda_lo);
      xs = next;
      T Sd;
      T Sa;
      T Sad;
      const T S = legendre_series(co.m, co.p, co.d, xs, &Sd, &Sa, &Sad, pol, function);
      // The scale A needs S or S' to full relative accuracy, not just within its conditioning.
      const bool use_value = Sa <= 2 * abs(S);
      if (use_value || Sad <= 2 * abs(Sd))
      {
         const T A = use_value ? T(S / us) : T(Sd / usd);
         value = A * u;
         d = A * ud;
         matched = (BOOST_MATH_NAMESPACE::isfinite)(A);
      }
   }
   if (!matched)
   {
      return policies::raise_evaluation_error<T>(function, "No accurate matching point for the angular function at eta = %1%", eta, pol);
   }
   if (flip)
   {
      value *= parity;
      d *= -parity;
   }
   if (derivative)
   {
      *derivative = d;
   }
   return value;
}

// z = c sqrt(x² - τ) as z + *z_lo, so that the phase of oscillating functions of z is not lost to rounding:
// for large cξ an ulp of z is many ulps of cos(z). x² and the products use error-free transforms.
template <class T>
T phase(const T& c, const T& x, const T& tau, T* z_lo)
{
   BOOST_MATH_STD_USING
   T q;
   T q_lo;
   if (tau == 1 && x < 2)
   {
      // (x - 1)(x + 1): x - 1 is exact here, x + 1 and the product are not.
      const T a = x - 1;
      const T s = x + 1;
      const T s_lo = (x - (s - 1)) + (1 - (s - (s - 1)));
      q = fma(a, s, T(0));
      q_lo = fma(a, s, -q) + a * s_lo;
   }
   else
   {
      const T p = fma(x, x, T(0));
      const T p_lo = fma(x, x, -p);
      q = p - tau;
      const T b = q - p;
      q_lo = (p - (q - b)) + (-tau - b) + p_lo;
   }
   const T root = sqrt(q);
   const T root_lo = (fma(-root, root, q) + q_lo) / (2 * root);
   const T z = fma(c, root, T(0));
   *z_lo = fma(c, root, -z) + c * root_lo;
   return z;
}

// Corrects b_l(z), l = 0..L-1 (spherical Bessel or Neumann functions) to b_l(z + z_lo) to first order,
// with b_l' = (l/z) b_l - b_{l+1}.
template <class T>
void shift_phase(std::vector<T>& b, const T& z, const T& z_lo)
{
   if (z_lo == 0)
   {
      return;
   }
   for (std::size_t l = 0; l + 1 < b.size(); ++l)
   {
      b[l] += z_lo * (T(l) / z * b[l] - b[l + 1]);
   }
}

// Spherical Bessel functions j_l(z), l = 0..L, by Miller's backward recurrence normalized with j_0 or j_1.
template <class T>
std::vector<T> sph_bessel_j_all(unsigned L, const T& z)
{
   BOOST_MATH_STD_USING
   std::vector<T> j(L + 1, T(0));
   if (z * z < tools::epsilon<T>())
   {
      // j_l(z) = z^l/(2l+1)!! (1 + O(z²)); the recurrence below would overflow for tiny z.
      T term = 1;
      for (unsigned l = 0; l <= L && term != 0; ++l)
      {
         j[l] = term;
         term *= z / T(2 * l + 3);
      }
      return j;
   }
   const unsigned start = L + 20 + static_cast<unsigned>(z) + static_cast<unsigned>(tools::digits<T>() / 2);
   const T big = sqrt(tools::max_value<T>()) * tools::epsilon<T>();
   T f_next = 0;
   T f = tools::min_value<T>() / tools::epsilon<T>();
   for (unsigned l = start; l > 0; --l)
   {
      const T f_prev = (2 * T(l) + 1) / z * f - f_next;
      f_next = f;
      f = f_prev;
      if (l - 1 <= L)
      {
         j[l - 1] = f;
      }
      if (abs(f) > big)
      {
         f /= big;
         f_next /= big;
         for (unsigned i = l - 1; i <= L; ++i)
         {
            j[i] /= big;
         }
      }
   }
   const T j0 = sin(z) / z;
   const T j1 = (sin(z) / z - cos(z)) / z;
   const T scale = abs(j0) >= abs(j1) ? T(j0 / j[0]) : T(j1 / j[1]);
   for (T& v : j)
   {
      v *= scale;
   }
   return j;
}

// Oblate R^(1) from the expansion in spherical Bessel functions (DLMF 30.11.3 with γ = -ic, z = iξ):
//   R = (ξ² + 1)^{m/2} c^m Σ' (-1)^{(k-r)/2} v_k j_{m+k}(cξ)/(cξ)^m / Σ' v_k,  v_k = d_k (2m+k)!/k!.
// The oblate d_k all have one sign, so Σ' v_k does not cancel (the prolate one does; see prolate_radial_1).
// Writing ξ^{-m} j_{m+k}(cξ) = c^m j_{m+k}(z)/z^m keeps the small-ξ behavior free of cancellation.
template <class T>
T oblate_radial_1(const coefficients<T>& co, const T& c, const T& xi, T* derivative)
{
   BOOST_MATH_STD_USING
   const unsigned m = co.m;
   const unsigned p = co.p;
   const unsigned r = co.n - co.m;
   const std::size_t size = co.d.size();
   const std::vector<T> F = factorial_ratios<T>(m, r, p, size);
   const T z = fma(c, xi, T(0));
   std::vector<T> j = sph_bessel_j_all(m + p + 2 * static_cast<unsigned>(size - 1) + 2, z);
   shift_phase(j, z, T(fma(c, xi, -z)));
   T zm = 1;           // z^m
   T odd = 1;          // (2m+1)!!
   T cm = 1;           // c^m
   for (unsigned i = 0; i < m; ++i)
   {
      zm *= z;
      odd *= T(2 * i + 3);
      cm *= c;
   }
   T norm = 0;
   T sum = 0;
   T dsum = 0;
   for (std::size_t i = 0; i < size; ++i)
   {
      const unsigned k = p + 2 * static_cast<unsigned>(i);
      const unsigned l = m + k;
      const T v = co.d[i] * F[i];
      norm += v;
      const T sv = (((k - r) / 2) & 1u) ? T(-v) : v;
      if (z * z < tools::epsilon<T>())
      {
         // j_l(z)/z^m = z^k/(2l+1)!! (1 + O(z²)), k = l - m, and its z-derivative k z^{k-1}/(2l+1)!!.
         T zk = 1;
         T odd_l = odd;
         for (unsigned e = 0; e < k; ++e)
         {
            zk *= z;
            odd_l *= T(2 * (m + e) + 3);
         }
         sum += sv * zk / odd_l;
         dsum += k == 0 ? T(0) : T(sv * T(k) * (z == 0 ? T(k == 1 ? 1 : 0) : T(zk / z)) / odd_l);
      }
      else
      {
         sum += sv * j[l] / zm;
         // d/dz [j_l / z^m] = ((l - m)/z j_l - j_{l+1}) / z^m
         dsum += sv * (T(l - m) / z * j[l] - j[l + 1]) / zm;
      }
   }
   const T q = xi * xi + 1;
   const T f = pow(q, T(m) / 2) * cm;
   if (derivative)
   {
      *derivative = (T(m) * xi / q * f * sum + f * c * dsum) / norm;
   }
   return f * sum / norm;
}

// Prolate R^(1) from the expansion in j_{m+k}(c sqrt(ξ²-1)) normalized by Flammer's sum
// (Van Buren & Boisvert 2002); unlike the expansion in j_{m+k}(cξ) its denominator does not cancel:
//   R = (ξ/sqrt(ξ²-1))^p Σ' d_k a_k j_{m+k}(z) / Σ' (-1)^{(k-r)/2} d_k a_k,
//   a_k = (k+2m-1)!!/k!! (p = 0) or (k+2m)!!/(k-1)!! (p = 1).
template <class T>
T prolate_radial_1(const coefficients<T>& co, const T& c, const T& xi, T* derivative)
{
   BOOST_MATH_STD_USING
   const unsigned m = co.m;
   const unsigned p = co.p;
   const unsigned r = co.n - co.m;
   const std::size_t size = co.d.size();
   const T q = (xi - 1) * (xi + 1);
   const T root = sqrt(q);
   T z_lo;
   const T z = phase(c, xi, T(1), &z_lo);
   std::vector<T> j = sph_bessel_j_all(m + p + 2 * static_cast<unsigned>(size - 1) + 2, z);
   shift_phase(j, z, z_lo);
   std::vector<T> a(size);
   const std::size_t i0 = (r - p) / 2;
   a[i0] = 1;
   for (std::size_t i = i0 + 1; i < size; ++i)
   {
      const unsigned k = p + 2 * static_cast<unsigned>(i);
      a[i] = a[i - 1] * (p == 0 ? T(k + 2 * m - 1) / T(k) : T(k + 2 * m) / T(k - 1));
   }
   for (std::size_t i = i0; i > 0; --i)
   {
      const unsigned k = p + 2 * static_cast<unsigned>(i);
      a[i - 1] = a[i] / (p == 0 ? T(k + 2 * m - 1) / T(k) : T(k + 2 * m) / T(k - 1));
   }
   T den = 0;
   for (std::size_t i = 0; i < size; ++i)
   {
      const unsigned k = p + 2 * static_cast<unsigned>(i);
      const T w = co.d[i] * a[i];
      den += ((((k - r) / 2) & 1u) ? -w : w);
   }
   T num = 0;
   T dnum = 0;
   if (z == 0)
   {
      // ξ = 1. With j_l(z) ~ z^l/(2l+1)!!, only the lowest orders survive; the derivative is infinite for
      // m = 1, which the caller reports as an overflow.
      const T infinite = std::numeric_limits<T>::has_infinity ? std::numeric_limits<T>::infinity() : tools::max_value<T>();
      for (std::size_t i = 0; i < size; ++i)
      {
         const unsigned l = m + p + 2 * static_cast<unsigned>(i);
         const T w = co.d[i] * a[i];
         if (p == 0)
         {
            // R = N, R' = c² lim j_l'(z)/z: -1/3 (l = 0), 1/0 (l = 1), 2/15 (l = 2)
            num += l == 0 ? w : T(0);
            dnum += w * (l == 0 ? T(-1) / 3 : l == 1 ? infinite : l == 2 ? T(2) / 15 : T(0));
         }
         else
         {
            // R = c lim j_l/z, R' = c [lim j_l/z + c² lim (((l-1)/z) j_l - j_{l+1})/z²]
            num += l == 1 ? T(w / 3) : T(0);
            dnum += w * (l == 1 ? T(-1) / 15 : l == 2 ? infinite : l == 3 ? T(2) / 105 : T(0));
         }
      }
      if (derivative)
      {
         *derivative = p == 0 ? T(c * c * dnum / den) : T(c * (num + c * c * dnum) / den);
      }
      return p == 0 ? T(num / den) : T(c * num / den);
   }
   for (std::size_t i = 0; i < size; ++i)
   {
      const unsigned l = m + p + 2 * static_cast<unsigned>(i);
      const T w = co.d[i] * a[i];
      // p = 0: N and N'(z), with j_l' = (l/z) j_l - j_{l+1}.
      // p = 1: Ñ = N/z and M̃ = Σ w (((l-1)/z) j_l - j_{l+1}), so that dÑ/dz = M̃/z.
      num += w * j[l] / (p ? z : T(1));
      dnum += w * (T(l - p) / z * j[l] - j[l + 1]);
   }
   if (p == 0)
   {
      if (derivative)
      {
         // dz/dξ = cξ/sqrt(ξ²-1)
         *derivative = c * xi / root * dnum / den;
      }
      return num / den;
   }
   // R = (ξ/sqrt(ξ²-1)) N/Den = cξ Ñ/Den and R' = c (Ñ + ξ² M̃/q)/Den, with no cancellation near ξ = 1.
   if (derivative)
   {
      *derivative = c * (num + xi * xi * dnum / q) / den;
   }
   return c * xi * num / den;
}

// R^(2) from the expansion in spherical Neumann functions (DLMF 30.11.3):
//   R = ((ξ² - σ)/ξ²)^{m/2} Σ' (-1)^{(k-r)/2} v_k y_{m+k}(cξ) / Σ' v_k,  v_k = d_k (2m+k)!/k!.
// For prolate functions with large c the sums cancel (see prolate_radial_1), and near the singular
// point the series converges slowly; *condition reports the cancellation and false means no convergence.
template <class T>
bool neumann_radial_2(const coefficients<T>& co, const T& c, const T& xi, T* value, T* derivative, T* condition)
{
   BOOST_MATH_STD_USING
   const unsigned m = co.m;
   const unsigned p = co.p;
   const unsigned r = co.n - co.m;
   const T sigma = co.s > 0 ? T(1) : T(-1);
   const std::size_t size = co.d.size();
   const std::vector<T> F = factorial_ratios<T>(m, r, p, size);
   const T z = fma(c, xi, T(0));
   const unsigned L = m + p + 2 * static_cast<unsigned>(size - 1) + 2;
   // y_l(z) by forward recurrence, which is stable; stop before overflow.
   std::vector<T> y(L + 1, T(0));
   y[0] = -cos(z) / z;
   y[1] = (y[0] - sin(z)) / z;
   {
      // Shift the seeds to z + z_lo (y_0' = -y_1, y_1' = y_0 - 2 y_1 / z); the recurrence keeps the phase.
      const T z_lo = fma(c, xi, -z);
      const T y0 = y[0];
      y[0] -= z_lo * y[1];
      y[1] += z_lo * (y0 - 2 * y[1] / z);
   }
   unsigned top = 1;
   const T limit = tools::max_value<T>() * z / (8 * (T(L) + 2));
   while (top < L && abs(y[top]) < limit)
   {
      y[top + 1] = (2 * T(top) + 1) / z * y[top] - y[top - 1];
      ++top;
   }
   T norm = 0;
   T abs_norm = 0;
   T sum = 0;
   T dsum = 0;
   T abs_sum = 0;
   T last = 0;
   for (std::size_t i = 0; i < size; ++i)
   {
      const unsigned k = p + 2 * static_cast<unsigned>(i);
      const unsigned l = m + k;
      const T v = co.d[i] * F[i];
      norm += v;
      abs_norm += abs(v);
      if (l + 1 > top)
      {
         return false;
      }
      const T sv = (((k - r) / 2) & 1u) ? T(-v) : v;
      const T dy = T(l) / z * y[l] - y[l + 1];
      sum += sv * y[l];
      dsum += sv * dy;
      abs_sum += abs(sv * y[l]) + abs(sv * dy);
      last = abs(sv * y[l]) + abs(sv * dy);
   }
   if (!(last <= tools::epsilon<T>() * (abs(sum) + abs(dsum)) / 4))
   {
      return false;
   }
   const T q = sigma > 0 ? T((xi - 1) * (xi + 1)) : T(xi * xi + 1);
   const T f = pow(q / (xi * xi), T(m) / 2);
   *value = f * sum / norm;
   *derivative = (sigma * T(m) / (xi * q) * f * sum + f * c * dsum) / norm;
   *condition = (std::max)(T(abs_sum / sqrt(sum * sum + dsum * dsum)), T(abs_norm / abs(norm)));
   return true;
}

// A minimal complex type: std::complex is not specified for multiprecision types.
template <class T>
struct cplx
{
   T re;
   T im;
   cplx(const T& r = T(0), const T& i = T(0)) : re(r), im(i) {}
   cplx& operator+=(const cplx& o) { re += o.re; im += o.im; return *this; }
   cplx& operator-=(const cplx& o) { re -= o.re; im -= o.im; return *this; }
   friend cplx operator*(const cplx& a, const cplx& b) { return cplx(a.re * b.re - a.im * b.im, a.re * b.im + a.im * b.re); }
   friend cplx operator*(const cplx& a, const T& b) { return cplx(a.re * b, a.im * b); }
   friend cplx operator*(const T& b, const cplx& a) { return cplx(a.re * b, a.im * b); }
   friend cplx operator/(const cplx& a, const T& b) { return cplx(a.re / b, a.im / b); }
   friend cplx operator+(const cplx& a, const cplx& b) { return cplx(a.re + b.re, a.im + b.im); }
   friend cplx operator-(const cplx& a) { return cplx(-a.re, -a.im); }
   friend bool operator==(const cplx& a, const cplx& b) { return a.re == b.re && a.im == b.im; }
   cplx conj() const { return cplx(re, -im); }
   friend cplx operator/(const cplx& a, const cplx& b) { const T d = b.re * b.re + b.im * b.im; return cplx((a.re * b.re + a.im * b.im) / d, (a.im * b.re - a.re * b.im) / d); }
   T magnitude() const { BOOST_MATH_STD_USING return sqrt(re * re + im * im); }
};

// Asymptotic expansion of R^(3) = R^(1) + i R^(2) = e^{iz} Σ_K b_K ξ^{-K-1}, b_0 = e^{-i(n+1)π/2}/c,
// with the phase z = c sqrt(ξ² - τ), τ = σ - λ/c², matched to the WKB wavenumber so that the
// coefficients do not grow with c (τ ≈ 1 for low prolate modes, τ ≈ 0 for oblate ones).
// With p = ξ² - τ and q = ξ² - σ, v = Σ b_K ξ^{-K-1} satisfies
//   p q² v'' + (2ic ξ q² p^{1/2} + 2ξ q p) v' + [-icτ q² p^{-1/2} + 2ic ξ² q p^{1/2} + λτ q - σm² p] v = 0.
// Expanding p^{±1/2} in powers of ξ^{-2} turns this into a recurrence for b_K. Returns false if the
// series does not reach full accuracy without cancellation at this ξ.
template <class T>
bool radial_asymptotic(const coefficients<T>& co, const T& c, const T& xi, T* R2, T* R2d, T* amplification, std::vector<cplx<T>>* series = nullptr)
{
   BOOST_MATH_STD_USING
   typedef cplx<T> C;
   const T sigma = co.s > 0 ? T(1) : T(-1);
   const T tau = sigma - co.lambda / (c * c);
   if (xi * xi <= 2 * abs(tau))
   {
      return false;  // outside the region where the series in ξ^{-2} converges well
   }
   const T mm = T(co.m);
   const C ic(T(0), c);
   const int max_terms = 400;
   const int J = max_terms / 2 + 8;
   // Laurent coefficients A_d[i] of ξ^{6-i}, for derivative orders d = 0, 1, 2.
   std::vector<C> A[3];
   for (auto& a : A)
   {
      a.assign(2 * J + 12, C());
   }
   auto add = [&](int d, int power, const C& coef) {
      const int i = 6 - power;
      if (i >= 0 && i < static_cast<int>(A[d].size()))
      {
         A[d][static_cast<std::size_t>(i)] += coef;
      }
   };
   // Polynomials q² = ξ⁴ - 2σξ² + 1, q = ξ² - σ, p = ξ² - τ.
   const T q2[3] = { 1, -2 * sigma, 1 };  // coefficients of ξ⁴, ξ², ξ⁰
   const T q1[2] = { 1, -sigma };         // ξ², ξ⁰
   const T p1[2] = { 1, -tau };           // ξ², ξ⁰
   // A2 = p q²
   for (int a = 0; a < 2; ++a)
      for (int b = 0; b < 3; ++b)
         add(2, 6 - 2 * a - 2 * b, C(T(p1[a] * q2[b])));
   // A1 polynomial part 2ξ q p, and A0 polynomial part λτ q - σ m² p
   for (int a = 0; a < 2; ++a)
   {
      for (int b = 0; b < 2; ++b)
      {
         add(1, 5 - 2 * a - 2 * b, C(T(2 * q1[a] * p1[b])));
      }
      add(0, 2 - 2 * a, C(T(co.lambda * tau * q1[a] - sigma * mm * mm * p1[a])));
   }
   // p^{1/2} = ξ Σ s_j ξ^{-2j}, p^{-1/2} = ξ^{-1} Σ t_j ξ^{-2j}
   T sj = 1;
   T tj = 1;
   for (int j = 0; j < J; ++j)
   {
      const int pj = 2 * j;
      for (int b = 0; b < 3; ++b)
      {
         // 2ic ξ q² p^{1/2} v'  and  -icτ q² p^{-1/2} v
         add(1, 6 - 2 * b - pj, T(2 * q2[b] * sj) * ic);
         add(0, 3 - 2 * b - pj, T(-tau * q2[b] * tj) * ic);
      }
      for (int a = 0; a < 2; ++a)
      {
         // 2ic ξ² q p^{1/2} v
         add(0, 5 - 2 * a - pj, T(2 * q1[a] * sj) * ic);
      }
      sj *= (T(0.5) - T(j)) / T(j + 1) * (-tau);
      tj *= (T(-0.5) - T(j)) / T(j + 1) * (-tau);
   }
   // The leading 2ic ξ⁶ v' and 2ic ξ⁵ v combine into -2icK b_K.
   A[1][0] -= T(2) * ic;
   A[0][1] -= T(2) * ic;
   std::vector<C> b;
   const unsigned quarter = (co.n + 1) & 3u;  // e^{-i(n+1)π/2} = (-i)^{n+1}
   const C units[4] = { C(T(1)), C(T(0), T(-1)), C(T(-1)), C(T(0), T(1)) };
   b.push_back(units[quarter] / c);
   const T inv = 1 / xi;
   C w = b[0] * inv;
   C dw = -b[0] * inv * inv;
   T power = inv;   // ξ^{-K-1}
   T largest = b[0].magnitude() * inv;
   for (int K = 1; K < max_terms; ++K)
   {
      C sum;
      for (int d = 0; d < 3; ++d)
      {
         // ξ^{6-i} times the order-d derivative of b_jj ξ^{-jj-1} lands on ξ^{4-K} when jj = K + 1 - i - d.
         for (int i = 0; i < static_cast<int>(A[d].size()); ++i)
         {
            const int jj = K + 1 - i - d;
            if (jj < 0)
            {
               break;
            }
            const C& a = A[d][static_cast<std::size_t>(i)];
            if (jj >= K || a == C())
            {
               continue;
            }
            T fall = 1;
            for (int k = 1; k <= d; ++k)
            {
               fall *= -T(jj + k);
            }
            sum += a * b[static_cast<std::size_t>(jj)] * fall;
         }
      }
      // sum / (2icK) = -i sum / (2cK)
      b.push_back(C(sum.im, -sum.re) / T(2 * c * K));
      power *= inv;
      const C& bK = b.back();
      const T size = bK.magnitude() * power;
      largest = (std::max)(largest, size);
      w += bK * power;
      dw -= bK * T(K + 1) * power * inv;
      if (size <= tools::epsilon<T>() * w.magnitude() / 8 && b[b.size() - 2].magnitude() * power * xi <= tools::epsilon<T>() * w.magnitude())
      {
         if (largest > 16 * w.magnitude())
         {
            return false;  // too much cancellation
         }
         *amplification = largest / w.magnitude();
         if (series)
         {
            *series = b;
         }
         T z_lo;
         const T z = phase(c, xi, tau, &z_lo);
         const T root = z / c;
         // e^{i(z + z_lo)}
         const C e(cos(z) - z_lo * sin(z), sin(z) + z_lo * cos(z));
         const C R = e * w;
         const C Rd = e * (dw + ic * T(xi / root) * w);
         *R2 = R.im;
         *R2d = Rd.im;
         return true;
      }
      if (K > 2 * c * xi + 40)
      {
         return false;  // past the smallest term of the divergent series
      }
   }
   return false;
}

template <class T, class Policy>
T radial_2(const coefficients<T>& co, const T& c, const T& xi, bool oblate, T* derivative, const Policy& pol, const char* function)
{
   BOOST_MATH_STD_USING
   // R^(2): from the Neumann series or the asymptotic expansion, each at the nearest point where it converges,
   // whichever promises the smaller error (its cancellation plus the integration steps back to ξ); then
   // integrated inward.
   if (!oblate && c > 8 && xi - 1 < T(1) / 8)
   {
      // The integration below leaves a few ε times R1, from the start value and from each step, which its
      // rescale cannot remove; near ξ = 1, where R1 is several times R2, that is far outside the conditioning
      // in ξ - 1. So, as for the eigenvalue, the amplifying part runs in double length, from where the
      // asymptotic series is good to about ε of the amplitude; it adds no multiple of R2, so needs no rescale.
      T a = (std::min)(T(1 + c / 8), T(4));
      T A2;
      T A2d;
      T amplification;
      for (unsigned attempt = 0; !radial_asymptotic(co, c, a, &A2, &A2d, &amplification); ++attempt)
      {
         if (attempt == 200)
         {
            return policies::raise_evaluation_error<T>(function, "Asymptotic expansion failed to converge at xi = %1%", xi, pol);
         }
         a = 1 + 3 * (a - 1) / 2;
      }
      twofold<T> R = A2;
      twofold<T> Rd = A2d;
      taylor_continue(twofold<T>(co.lambda, co.lambda_lo), co.m, 1, twofold<T>(c) * twofold<T>(c), a, xi, &R, &Rd);
      if (derivative)
      {
         *derivative = Rd.hi + Rd.lo;
      }
      return R.hi + R.lo;
   }
   const T wavenumber = sqrt((std::max)(T(c * c), T(abs(co.lambda))) + 1);
   auto steps = [&](const T& from) {
      return abs(from - xi) * wavenumber + (oblate ? T(0) : T(log2((from - 1) / (xi - 1))));
   };
   T start = (std::max)(xi, T(2));
   T R2;
   T R2d;
   T condition = 0;
   bool neumann = false;
   for (; start <= 8 * (std::max)(xi, T(2)); start *= 2)
   {
      if (neumann_radial_2(co, c, start, &R2, &R2d, &condition) && condition < 64)
      {
         neumann = true;
         break;
      }
   }
   if (!neumann || condition > 4)
   {
      T a = oblate ? (std::max)(xi, T(1 / c)) : xi;
      T A2;
      T A2d;
      T amplification;
      bool found = false;
      for (unsigned attempt = 0; attempt <= 200; ++attempt)
      {
         if (radial_asymptotic(co, c, a, &A2, &A2d, &amplification))
         {
            found = true;
            break;
         }
         a = oblate ? T(a * 3 / 2) : T(1 + 3 * (a - 1) / 2);
      }
      if (!found && !neumann)
      {
         return policies::raise_evaluation_error<T>(function, "Asymptotic expansion failed to converge at xi = %1%", xi, pol);
      }
      if (found && (!neumann || amplification + steps(a) < condition + steps(start)))
      {
         start = a;
         R2 = A2;
         R2d = A2d;
      }
   }
   if (start != xi)
   {
      taylor_continue(co.lambda, co.m, oblate ? -1 : 1, T(c * c), start, xi, &R2, &R2d, co.lambda_lo);
      // The integration's rounding errors add α R1 + β R2. Since the Wronskian of R1 with itself vanishes,
      // the computed W(R1, R2) is exactly (1 + β) times 1/(c(ξ² ∓ 1)), so dividing by it removes β; α R1 is
      // what remains, and it shrinks wherever R2 dominates (toward ξ = 1 for prolate functions).
      T R1d;
      const T R1 = oblate ? oblate_radial_1(co, c, xi, &R1d) : prolate_radial_1(co, c, xi, &R1d);
      const T q = oblate ? T(xi * xi + 1) : T((xi - 1) * (xi + 1));
      const T scale = 1 / (c * q * (R1 * R2d - R1d * R2));
      if ((BOOST_MATH_NAMESPACE::isfinite)(scale) && abs(scale - 1) < T(1) / 1024)
      {
         R2 *= scale;
         R2d *= scale;
      }
   }
   if (derivative)
   {
      *derivative = R2d;
   }
   return R2;
}

// Oblate R^(2) near ξ = 0 as A u_e + B u_o, with u_e and u_o the even and odd solutions (u_e(0) = u_o'(0) = 1,
// u_e'(0) = u_o(0) = 0), so (ξ² + 1)(u_e u_o' - u_e' u_o) = 1. The Wronskian with R^(1), which is a multiple of
// u_e (n - m even) or of u_o (odd), gives B = 1/(c R1(0)) or A = -1/(c R1'(0)) exactly. The other coefficient,
// R2(0) or R2'(0), is exponentially small for large c and comes from matching at ξ0; it is known only to about
// an ulp of the scale there, so below that it is indistinguishable from 0 and is set to 0, which keeps R2 accurate
// relative to its own size near ξ = 0.
// Power series in u = 1/ξ, truncated to n terms.
template <class T>
std::vector<cplx<T>> series_multiply(const std::vector<cplx<T>>& a, const std::vector<cplx<T>>& b, std::size_t n)
{
   std::vector<cplx<T>> r(n);
   for (std::size_t i = 0; i < (std::min)(a.size(), n); ++i)
   {
      for (std::size_t j = 0; j < b.size() && i + j < n; ++j)
      {
         r[i + j] += a[i] * b[j];
      }
   }
   return r;
}

template <class T>
cplx<T> series_value(const std::vector<cplx<T>>& a, const T& u)
{
   cplx<T> r;
   for (std::size_t k = a.size(); k-- > 0;)
   {
      r = r * u + a[k];
   }
   return r;
}

// Gauss-Legendre nodes and weights on [-1, 1], by Newton's method on P_n.
template <class T>
void gauss_legendre(unsigned n, std::vector<T>& x, std::vector<T>& w)
{
   BOOST_MATH_STD_USING
   x.resize(n);
   w.resize(n);
   for (unsigned i = 0; i < n; ++i)
   {
      T t = cos(constants::pi<T>() * (T(i) + T(0.75)) / (T(n) + T(0.5)));
      T dp = 0;
      for (int it = 0; it < 100; ++it)
      {
         T p0 = 1;
         T p1 = t;
         for (unsigned k = 2; k <= n; ++k)
         {
            const T p2 = ((2 * T(k) - 1) * t * p1 - (T(k) - 1) * p0) / T(k);
            p0 = p1;
            p1 = p2;
         }
         dp = T(n) * (t * p1 - p0) / (t * t - 1);
         const T step = p1 / dp;
         t -= step;
         if (abs(step) <= tools::epsilon<T>() * 2)
         {
            break;
         }
      }
      x[i] = t;
      w[i] = 2 / ((1 - t * t) * dp * dp);
   }
}

// The exponentially small coefficient of the oblate R2 at ξ = 0 (R2(0) for even n - m, R2'(0) for odd), which
// matching cannot resolve. For large c the eigenvalues pair up: λ_e (n - m even) and λ_o = λ_{e+1}, split by an
// exponentially small Δ = λ_o - λ_e. With f = R2_n and g the first-kind function of the partner, the radial
// equation gives d/dξ[(ξ²+1) W(f, g)] = (λ_g - λ_f) f g, and W vanishes at infinity, where f and g have the
// same phase, so the coefficient is -Δ ∫_0^∞ f g dξ divided by g(0) or g'(0). Δ comes from the same identity for
// the angular functions on [0, 1], Δ = Ps_e(0) Ps_o'(0) / ∫_0^1 Ps_e Ps_o dη, whose values at 0
// are computed to full relative accuracy. ∫_0^X f g is integrated exactly along shared Taylor steps from ξ = 0
// (where f starts from the matched coefficient, whose error hardly affects the integral), and the tail from
// the asymptotic expansions of both functions.
template <class T, class Policy>
bool oblate_small_coefficient(const coefficients<T>& co, const T& c, const T& exact, T* result, const Policy& pol, const char* function)
{
   BOOST_MATH_STD_USING
   typedef cplx<T> C;
   const unsigned n = co.n;
   const unsigned m = co.m;
   const bool even = co.p == 0;
   const coefficients<T> partner = compute_coefficients(even ? n + 1 : n - 1, m, co.s, pol, function);
   const coefficients<T>& ce = even ? co : partner;
   const coefficients<T>& cpo = even ? partner : co;
   // Δ from the angular functions.
   T Se_d;
   const T Se = angular(ce, T(0), &Se_d, pol, function);
   T So_d;
   angular(cpo, T(0), &So_d, pol, function);
   std::vector<T> gx;
   std::vector<T> gw;
   gauss_legendre(20u, gx, gw);
   T J = 0;
   T lo = 0;
   const T width = 1 / (4 * c);
   while (lo < 1)
   {
      const T hi = (1 - lo) / 2 <= width ? T(1) : T(lo + (1 - lo) / 2);
      for (unsigned i = 0; i < gx.size(); ++i)
      {
         const T eta = (lo + hi) / 2 + (hi - lo) / 2 * gx[i];
         const T a = legendre_series(m, ce.p, ce.d, eta, static_cast<T*>(nullptr), static_cast<T*>(nullptr), static_cast<T*>(nullptr), pol, function);
         const T b = legendre_series(m, cpo.p, cpo.d, eta, static_cast<T*>(nullptr), static_cast<T*>(nullptr), static_cast<T*>(nullptr), pol, function);
         J += (hi - lo) / 2 * gw[i] * a * b;
      }
      lo = hi;
   }
   const T delta = Se * So_d / J;
   // g, the partner's first-kind function, and its value or slope at 0.
   T g0d;
   const T g0 = oblate_radial_1(partner, c, T(0), &g0d);
   // A point X beyond which both asymptotic series hold, with their coefficients.
   std::vector<C> bf;
   std::vector<C> bg;
   // The asymptotic series and the tail's integration by parts are good to about e^{-2cX}.
   T X = (std::max)(T(1), T((10 - T(0.6) * log(tools::epsilon<T>())) / c));
   T dummy;
   T amp;
   for (unsigned attempt = 0;; ++attempt)
   {
      if (radial_asymptotic(co, c, X, &dummy, &dummy, &amp, &bf) && radial_asymptotic(partner, c, X, &dummy, &dummy, &amp, &bg))
      {
         break;
      }
      if (attempt > 100)
      {
         return false;
      }
      X *= T(1.5);
   }
   // On [0, X], f = A u_e + B u_o with u_e, u_o the even and odd solutions at 0, so ∫_0^X f g = A K + B M with
   // K = ∫ u_e g and M = ∫ u_o g, integrated exactly along shared Taylor steps.
   T ue = 1;
   T ued = 0;
   T uo = 0;
   T uod = 1;
   T g = even ? T(0) : g0;
   T gd = even ? g0d : T(0);
   const T c2 = c * c;
   const T wavelength = 1 / sqrt((std::max)({ c2, T(abs(co.lambda)), T(abs(partner.lambda)) }) + 1);
   std::vector<T> ce_u;
   std::vector<T> co_u;
   std::vector<T> g_u;
   T K = 0;
   T M = 0;
   auto product = [](const std::vector<T>& x, const std::vector<T>& y) {
      T sum = 0;
      for (std::size_t i = 0; i < x.size(); ++i)
      {
         for (std::size_t j = 0; j < y.size(); ++j)
         {
            sum += x[i] * y[j] / T(i + j + 1);
         }
      }
      return sum;
   };
   for (T a = 0; a != X;)
   {
      const T next = taylor_next(-1, wavelength, a, X);
      const T h = next - a;
      taylor_step(co.lambda, m, -1, c2, a, h, &ue, &ued, ce_u, co.lambda_lo);
      taylor_step(co.lambda, m, -1, c2, a, h, &uo, &uod, co_u, co.lambda_lo);
      taylor_step(partner.lambda, m, -1, c2, a, h, &g, &gd, g_u, partner.lambda_lo);
      K += h * product(ce_u, g_u);
      M += h * product(co_u, g_u);
      a = next;
   }
   // The tail from f = Im(e^{i z_f} v_f), g = Re(e^{i z_g} v_g), z = c sqrt(ξ² - τ), v = Σ b_K u^{K+1}:
   // f g = Im(e^{i(z_f - z_g)} v_f conj(v_g))/2 + Im(e^{i(z_f + z_g)} v_f v_g)/2.
   const std::size_t L = (std::max)(bf.size(), bg.size()) + 3;
   std::vector<C> vf(L);
   std::vector<C> vg(L);
   std::vector<C> vgc(L);
   for (std::size_t k = 0; k < bf.size() && k + 1 < L; ++k)
   {
      vf[k + 1] = bf[k];
   }
   for (std::size_t k = 0; k < bg.size() && k + 1 < L; ++k)
   {
      vg[k + 1] = bg[k];
      vgc[k + 1] = bg[k].conj();
   }
   const T tf = -1 - co.lambda / c2;
   const T tg = -1 - partner.lambda / c2;
   // δ = z_f - z_g = c Σ_{j≥1} binomial(1/2, j) ((-τ_f)^j - (-τ_g)^j) u^{2j-1}, and ψ' = (z_f + z_g)'.
   std::vector<C> idelta(L);
   std::vector<C> dpsi(L);
   {
      T bin_half = 1;
      T bin_mhalf = 1;
      T pf = 1;
      T pg = 1;
      for (std::size_t j = 0; 2 * j < L; ++j)
      {
         if (j > 0 && 2 * j - 1 < L)
         {
            idelta[2 * j - 1] = C(T(0), c * bin_half * (pf - pg));
         }
         dpsi[2 * j] = C(c * bin_mhalf * (pf + pg));
         bin_half *= (T(0.5) - T(j)) / T(j + 1);
         bin_mhalf *= (T(-0.5) - T(j)) / T(j + 1);
         pf *= -tf;
         pg *= -tg;
      }
   }
   // e^{iδ} by the power-series exponential, E_k = (1/k) Σ_{j=1}^{k} j δ_j E_{k-j}.
   std::vector<C> E(L);
   E[0] = C(T(1));
   for (std::size_t k = 1; k < L; ++k)
   {
      C acc;
      for (std::size_t j = 1; j <= k; ++j)
      {
         acc += T(j) * idelta[j] * E[k - j];
      }
      E[k] = acc / T(k);
   }
   const T u = 1 / X;
   T tail = 0;
   {
      const std::vector<C> slow = series_multiply(E, series_multiply(vf, vgc, L), L);
      T up = u;  // X^{1-k}
      for (std::size_t k = 2; k < L; ++k)
      {
         up *= u;
         tail += slow[k].im * (up * X) / (2 * T(k - 1));
      }
   }
   {
      // ∫_X^∞ e^{iψ} q = -e^{iψ(X)} Σ_k T_k(X), T_0 = q/(iψ'), T_{k+1} = -T_k'/(iψ').
      std::vector<C> idpsi(L);
      for (std::size_t k = 0; k < L; ++k)
      {
         idpsi[k] = C(-dpsi[k].im, dpsi[k].re);  // i ψ'
      }
      std::vector<C> inv(L);
      inv[0] = C(T(1)) / idpsi[0];
      for (std::size_t k = 1; k < L; ++k)
      {
         C acc;
         for (std::size_t j = 1; j <= k; ++j)
         {
            acc += idpsi[j] * inv[k - j];
         }
         inv[k] = -(acc / idpsi[0]);
      }
      std::vector<C> Tk = series_multiply(series_multiply(vf, vg, L), inv, L);
      C total;
      T previous = tools::max_value<T>();
      for (unsigned k = 0; k < 60; ++k)
      {
         const C value = series_value(Tk, u);
         if (value.magnitude() > previous)
         {
            break;
         }
         total += value;
         previous = value.magnitude();
         if (previous <= tools::epsilon<T>() * total.magnitude())
         {
            break;
         }
         // -T_k', with d/dξ u^j = -j u^{j+1}
         std::vector<C> d(L);
         for (std::size_t j = 0; j + 1 < L; ++j)
         {
            d[j + 1] = T(j) * Tk[j];
         }
         Tk = series_multiply(d, inv, L);
      }
      T zf_lo;
      T zg_lo;
      const T psi = phase(c, X, tf, &zf_lo) + phase(c, X, tg, &zg_lo);
      const T psi_lo = zf_lo + zg_lo;
      const C e(cos(psi) - psi_lo * sin(psi), sin(psi) + psi_lo * cos(psi));
      tail += -(e * total).im / 2;
   }
   // -W(0) = (λ_g - λ_f) ∫_0^∞ f g, solved for the unknown coefficient:
   //   n - m even: A = -Δ (B M + tail) / (g'(0) + Δ K);  odd: B = -Δ (A K + tail) / (g(0) + Δ M).
   *result = even ? T(-delta * (exact * M + tail) / (g0d + delta * K)) : T(-delta * (exact * K + tail) / (g0 + delta * M));
   return (BOOST_MATH_NAMESPACE::isfinite)(*result);
}

template <class T, class Policy>
bool oblate_radial_2_near_zero(const coefficients<T>& co, const T& c, const T& xi, const T& xi0, T* value, T* derivative, const Policy& pol, const char* function)
{
   BOOST_MATH_STD_USING
   T R0d;
   const T R0 = radial_2(co, c, xi0, true, &R0d, pol, function);
   const T c2 = c * c;
   T ue = 1, ued = 0, uo = 0, uod = 1;
   taylor_continue(co.lambda, co.m, -1, c2, T(0), xi, &ue, &ued, co.lambda_lo);
   taylor_continue(co.lambda, co.m, -1, c2, T(0), xi, &uo, &uod, co.lambda_lo);
   T ve = ue, ved = ued, vo = uo, vod = uod;
   taylor_continue(co.lambda, co.m, -1, c2, xi, xi0, &ve, &ved, co.lambda_lo);
   taylor_continue(co.lambda, co.m, -1, c2, xi, xi0, &vo, &vod, co.lambda_lo);
   const T w = xi0 * xi0 + 1;
   T A = w * (R0 * vod - R0d * vo);
   T B = w * (R0d * ve - R0 * ved);
   // The size of the terms that make up each matched coefficient.
   const T scale_A = w * (abs(R0 * vod) + abs(R0d * vo));
   const T scale_B = w * (abs(R0d * ve) + abs(R0 * ved));
   // The coefficient fixed by the Wronskian with R1 is exact. The other is exponentially small for large c and
   // cancels in matching; when it is well below the size of its terms, it comes from the eigenvalue splitting
   // instead (oblate_small_coefficient). Otherwise the inward integration is accurate (to a few ε where the
   // splitting formula is not, at small c), and the caller uses it.
   const T rough = co.p == 0 ? A : B;
   if (!(abs(rough) < (co.p == 0 ? scale_A : scale_B) / 16))
   {
      return false;
   }
   T R1d;
   const T R1 = oblate_radial_1(co, c, T(0), &R1d);
   const T exact = co.p == 0 ? T(1 / (c * R1)) : T(-1 / (c * R1d));
   T small;
   if (!oblate_small_coefficient(co, c, exact, &small, pol, function))
   {
      return false;
   }
   if (co.p == 0)
   {
      A = small;
      B = exact;
   }
   else
   {
      A = exact;
      B = small;
   }
   *derivative = A * ued + B * uod;
   *value = A * ue + B * uo;
   return true;
}

template <class T, class Policy>
T radial(const coefficients<T>& co, const T& c, const T& xi, int kind, bool oblate, T* derivative, const Policy& pol, const char* function)
{
   BOOST_MATH_STD_USING
   if (kind == 1)
   {
      return oblate ? oblate_radial_1(co, c, xi, derivative) : prolate_radial_1(co, c, xi, derivative);
   }
   const T xi0 = (std::min)(T(1), T(4 / c));
   T value;
   T d;
   if (oblate && xi < xi0 && oblate_radial_2_near_zero(co, c, xi, xi0, &value, &d, pol, function))
   {
      if (derivative)
      {
         *derivative = d;
      }
      return value;
   }
   return radial_2(co, c, xi, oblate, derivative, pol, function);
}

enum class quantity { eigenvalue, angular, radial_1, radial_2 };

// DLMF's λ_n^m(γ²) (30.3); the recurrence uses Flammer's, which is λ + γ².
template <class T>
T dlmf_eigenvalue(const coefficients<T>& co)
{
   const twofold<T> t = twofold_sum(co.lambda, T(-co.s));
   return t.hi + T(t.lo + co.lambda_lo);
}

// The eigenvalue and coefficients of one (n, m, c), computed once, and the functions that use them.
template <class T, class Policy>
class wave
{
public:
   // Everything depends on c², except that R^(1) and R^(2) change by (-1)^n and (-1)^(n+1) with the sign of c.
   wave(unsigned n, unsigned m, const T& c, bool oblate, const char* function)
      : c_(c < 0 ? T(-c) : c), oblate_(oblate), flip_1_(c < 0 && n % 2 == 1), flip_2_(c < 0 && n % 2 == 0)
   {
      if (n < m)
      {
         invalid_ = policies::raise_domain_error<T>(function, "Degree n must be >= order m, but got n = %1%", T(n), Policy());
      }
      else if (!(BOOST_MATH_NAMESPACE::isfinite)(c))
      {
         invalid_ = policies::raise_domain_error<T>(function, "Spheroidal parameter c must be finite, but got c = %1%", c, Policy());
      }
      else
      {
         valid_ = true;
         co_ = compute_coefficients(n, m, oblate ? T(-c * c) : T(c * c), Policy(), function);
      }
   }

   T operator()(quantity what, const T& x, T* derivative, const char* function) const
   {
      BOOST_MATH_STD_USING
      const Policy pol;
      if (!valid_)
      {
         return invalid_;
      }
      if (what == quantity::eigenvalue)
      {
         return dlmf_eigenvalue(co_);
      }
      if (!(BOOST_MATH_NAMESPACE::isfinite)(x))
      {
         return policies::raise_domain_error<T>(function, "Argument must be finite, but got %1%", x, pol);
      }
      if (what == quantity::angular)
      {
         if (abs(x) > 1)
         {
            return policies::raise_domain_error<T>(function, "Angular argument must be in [-1, 1], but got eta = %1%", x, pol);
         }
         return angular(co_, x, derivative, pol, function);
      }
      if (c_ == 0)
      {
         return policies::raise_domain_error<T>(function, "Radial functions need c > 0, but got c = %1%", c_, pol);
      }
      if (oblate_ ? x < 0 : x < 1)
      {
         return policies::raise_domain_error<T>(function, oblate_ ? "Oblate radial argument must be >= 0, but got xi = %1%" : "Prolate radial argument must be >= 1, but got xi = %1%", x, pol);
      }
      if (what == quantity::radial_2 && !oblate_ && x == 1)
      {
         return policies::raise_overflow_error<T>(function, nullptr, pol);
      }
      T d;
      const T value = radial(co_, c_, x, what == quantity::radial_2 ? 2 : 1, oblate_, &d, pol, function);
      if (!(BOOST_MATH_NAMESPACE::isfinite)(value) || (derivative && !(BOOST_MATH_NAMESPACE::isfinite)(d)))
      {
         return policies::raise_overflow_error<T>(function, nullptr, pol);
      }
      const bool flip = what == quantity::radial_2 ? flip_2_ : flip_1_;
      if (derivative)
      {
         *derivative = flip ? T(-d) : d;
      }
      return flip ? T(-value) : value;
   }

private:
   coefficients<T> co_;
   T c_;
   bool oblate_;
   bool flip_1_;
   bool flip_2_;
   bool valid_ = false;
   T invalid_ = 0;
};

template <class T, class Policy>
T evaluate(quantity what, bool oblate, unsigned n, unsigned m, const T& c, const T& x, T* derivative, const Policy&, const char* function)
{
   if (what == quantity::eigenvalue && n >= m && (BOOST_MATH_NAMESPACE::isfinite)(c))
   {
      // Without the coefficients' normalization, which the eigenvalue does not need.
      return dlmf_eigenvalue(compute_coefficients(n, m, oblate ? T(-c * c) : T(c * c), Policy(), function, false));
   }
   return wave<T, Policy>(n, m, c, oblate, function)(what, x, derivative, function);
}

template <class T1, class T2, class Policy>
inline typename tools::promote_args<T1, T2>::type function_api(quantity what, bool oblate, unsigned n, unsigned m, T1 c, T2 x,
   typename tools::promote_args<T1, T2>::type* derivative, const Policy&, const char* function)
{
   typedef typename tools::promote_args<T1, T2>::type result_type;
   typedef typename policies::evaluation<result_type, Policy>::type value_type;
   typedef typename policies::normalise<Policy, policies::promote_float<false>, policies::promote_double<false>>::type forwarding_policy;
   value_type d = 0;
   const value_type value = evaluate(what, oblate, n, m, static_cast<value_type>(c), static_cast<value_type>(x), derivative ? &d : nullptr, forwarding_policy(), function);
   if (derivative)
   {
      *derivative = policies::checked_narrowing_cast<result_type, forwarding_policy>(d, function);
   }
   return policies::checked_narrowing_cast<result_type, forwarding_policy>(value, function);
}

// The public classes' arguments and results, promoted and narrowed as for the free functions.
template <class Real, class Policy, bool Oblate>
class wave_function
{
   typedef typename policies::evaluation<Real, Policy>::type value_type;
   typedef typename policies::normalise<Policy, policies::promote_float<false>, policies::promote_double<false>>::type forwarding_policy;

public:
   Real eigenvalue() const
   {
      return get(quantity::eigenvalue, Real(0), nullptr, Oblate ? "boost::math::oblate_spheroidal<%1%>::eigenvalue()" : "boost::math::prolate_spheroidal<%1%>::eigenvalue()");
   }
   Real angular(Real eta, Real* derivative = nullptr) const
   {
      return get(quantity::angular, eta, derivative, Oblate ? "boost::math::oblate_spheroidal<%1%>::angular(%1%)" : "boost::math::prolate_spheroidal<%1%>::angular(%1%)");
   }
   Real radial_1(Real xi, Real* derivative = nullptr) const
   {
      return get(quantity::radial_1, xi, derivative, Oblate ? "boost::math::oblate_spheroidal<%1%>::radial_1(%1%)" : "boost::math::prolate_spheroidal<%1%>::radial_1(%1%)");
   }
   Real radial_2(Real xi, Real* derivative = nullptr) const
   {
      return get(quantity::radial_2, xi, derivative, Oblate ? "boost::math::oblate_spheroidal<%1%>::radial_2(%1%)" : "boost::math::prolate_spheroidal<%1%>::radial_2(%1%)");
   }

protected:
   wave_function(unsigned n, unsigned m, Real c)
      : w_(n, m, static_cast<value_type>(c), Oblate, Oblate ? "boost::math::oblate_spheroidal<%1%>::oblate_spheroidal(unsigned, unsigned, %1%)" : "boost::math::prolate_spheroidal<%1%>::prolate_spheroidal(unsigned, unsigned, %1%)")
   {
   }

private:
   Real get(quantity what, Real x, Real* derivative, const char* function) const
   {
      value_type d = 0;
      const value_type value = w_(what, static_cast<value_type>(x), derivative ? &d : nullptr, function);
      if (derivative)
      {
         *derivative = policies::checked_narrowing_cast<Real, forwarding_policy>(d, function);
      }
      return policies::checked_narrowing_cast<Real, forwarding_policy>(value, function);
   }

   wave<value_type, forwarding_policy> w_;
};

}} // namespace detail::spheroidal

// The functions of one (n, m, c) at many arguments. The eigenvalue and expansion coefficients, most of the
// cost of each free function below, are computed once, by the constructor.
BOOST_MATH_EXPORT template <class Real, class Policy = policies::policy<>>
class prolate_spheroidal : public detail::spheroidal::wave_function<Real, Policy, false>
{
public:
   prolate_spheroidal(unsigned n, unsigned m, Real c) : detail::spheroidal::wave_function<Real, Policy, false>(n, m, c) {}
};

BOOST_MATH_EXPORT template <class Real, class Policy = policies::policy<>>
class oblate_spheroidal : public detail::spheroidal::wave_function<Real, Policy, true>
{
public:
   oblate_spheroidal(unsigned n, unsigned m, Real c) : detail::spheroidal::wave_function<Real, Policy, true>(n, m, c) {}
};

// Each public function, with and without a policy; the derivative with respect to the argument is
// returned through the optional pointer, at little extra cost.
#define BOOST_MATH_SPHEROIDAL_FUNCTION(name, what, oblate) \
   BOOST_MATH_EXPORT template <class T1, class T2, class Policy> \
   inline typename tools::promote_args<T1, T2>::type name(unsigned n, unsigned m, T1 c, T2 x, \
      typename tools::promote_args<T1, T2>::type* derivative, const Policy& pol) \
   { \
      return detail::spheroidal::function_api(detail::spheroidal::quantity::what, oblate, n, m, c, x, derivative, pol, \
         "boost::math::" #name "<%1%>(unsigned, unsigned, %1%, %1%)"); \
   } \
   BOOST_MATH_EXPORT template <class T1, class T2> \
   inline typename tools::promote_args<T1, T2>::type name(unsigned n, unsigned m, T1 c, T2 x, \
      typename tools::promote_args<T1, T2>::type* derivative = nullptr) \
   { \
      return name(n, m, c, x, derivative, policies::policy<>()); \
   }

// Angular functions Ps_n^m(η, ±c²) (DLMF 30.4).
BOOST_MATH_SPHEROIDAL_FUNCTION(prolate_spheroidal_angular, angular, false)
BOOST_MATH_SPHEROIDAL_FUNCTION(oblate_spheroidal_angular, angular, true)
// Radial functions of the first and second kind (DLMF 30.11; Flammer's R_mn^(1), R_mn^(2)).
BOOST_MATH_SPHEROIDAL_FUNCTION(prolate_spheroidal_radial_1, radial_1, false)
BOOST_MATH_SPHEROIDAL_FUNCTION(prolate_spheroidal_radial_2, radial_2, false)
BOOST_MATH_SPHEROIDAL_FUNCTION(oblate_spheroidal_radial_1, radial_1, true)
BOOST_MATH_SPHEROIDAL_FUNCTION(oblate_spheroidal_radial_2, radial_2, true)
#undef BOOST_MATH_SPHEROIDAL_FUNCTION

// Eigenvalues λ_n^m(c²) and λ_n^m(-c²) (DLMF 30.3).
#define BOOST_MATH_SPHEROIDAL_EIGENVALUE(name, oblate) \
   BOOST_MATH_EXPORT template <class T, class Policy> \
   inline typename tools::promote_args<T>::type name(unsigned n, unsigned m, T c, const Policy& pol) \
   { \
      return detail::spheroidal::function_api(detail::spheroidal::quantity::eigenvalue, oblate, n, m, c, T(0), \
         static_cast<typename tools::promote_args<T>::type*>(nullptr), pol, "boost::math::" #name "<%1%>(unsigned, unsigned, %1%)"); \
   } \
   BOOST_MATH_EXPORT template <class T> \
   inline typename tools::promote_args<T>::type name(unsigned n, unsigned m, T c) \
   { \
      return name(n, m, c, policies::policy<>()); \
   }

BOOST_MATH_SPHEROIDAL_EIGENVALUE(prolate_spheroidal_eigenvalue, false)
BOOST_MATH_SPHEROIDAL_EIGENVALUE(oblate_spheroidal_eigenvalue, true)
#undef BOOST_MATH_SPHEROIDAL_EIGENVALUE

BOOST_MATH_NAMESPACE_END

#endif
