//  (C) Copyright Nick Thompson 2026.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_SPECIAL_STRUVE_HPP
#define BOOST_MATH_SPECIAL_STRUVE_HPP

// Struve functions H_nu(x) and modified Struve functions L_nu(x) of real order and argument (DLMF chapter 11).

#include <boost/math/tools/config.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/math/tools/promotion.hpp>
#include <boost/math/tools/precision.hpp>
#include <boost/math/tools/compensated_arithmetic.hpp>
#include <boost/math/policies/error_handling.hpp>
#include <boost/math/special_functions/math_fwd.hpp>
#include <boost/math/special_functions/fpclassify.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/special_functions/log1p.hpp>
#include <boost/math/special_functions/digamma.hpp>
#include <boost/math/special_functions/bessel.hpp>
#include <boost/math/special_functions/bessel_prime.hpp>
#include <boost/math/quadrature/exp_sinh.hpp>
#ifndef BOOST_MATH_BUILD_MODULE
#include <cmath>
#include <complex>
#include <cstddef>
#include <limits>
#include <vector>
#endif

BOOST_MATH_NAMESPACE_BEGIN

namespace detail { namespace struve {

template <class T>
bool is_nonpositive_integer(const T& z)
{
   BOOST_MATH_STD_USING
   return z <= 0 && floor(z) == z;
}

// Gamma(b) as a double-word number: (1 + tgamma1pm1(b0 - 1)), b0 in [1, 2), times or divided by the exact
// factors b0 + i in double-word arithmetic. The result is good to well under an ulp; tgamma itself is several
// ulps off in places (and tens of epsilon for cpp_bin_float_50).
template <class T, class Policy>
double_word<T> gamma_parts(const T& b, const Policy& pol)
{
   BOOST_MATH_STD_USING
   const double_word_arithmetic<T> dw;
   if (b >= 1 && b < 2)
   {
      return dw.fast_two_sum(T(1), BOOST_MATH_NAMESPACE::tgamma1pm1(T(b - 1), pol));
   }
   if (b >= 2)
   {
      const T n = floor(b) - 1;
      const T b0 = b - n;   // exact
      double_word<T> g = gamma_parts(b0, pol);
      for (T i = 0; i < n; ++i)
      {
         g = dw.multiply(g, T(b0 + i));   // b0 + i is exact
      }
      return g;
   }
   // b < 1: Gamma(b) = Gamma(b + n) / prod_{i<n} (b + i).
   const T n = ceil(1 - b);
   double_word<T> p = { T(1), T(0) };
   for (T i = 0; i < n; ++i)
   {
      p = dw.multiply(p, T(b + i));   // b + i is exact for |b| >= 1
   }
   return dw.divide(gamma_parts(T(b + n), pol), p);
}

// (x/2)^(nu + ca) / (sqrt(pi) Gamma(nu + cb)), zero when nu + cb is a pole of Gamma. The sums nu + ca and nu + cb
// are carried to first order in their rounding errors: rounding nu + 1 alone costs |log(x/2)| ulps in the power,
// and rounding nu + 3/2 |digamma(nu + 3/2)| ulps in Gamma.
template <class T, class Policy>
T scale(const T& nu, const T& ca, const T& cb, const T& x2, const Policy& pol)
{
   BOOST_MATH_STD_USING
   const T a = nu + ca;
   const T b = nu + cb;
   const T da = two_sum_error(nu, ca, a);
   const T db = two_sum_error(nu, cb, b);
   if (is_nonpositive_integer(b) && db == 0)
   {
      return 0;
   }
   T correction = 0;
   if (da != 0)
   {
      correction += da * log(x2);
   }
   if (db != 0)
   {
      correction -= db * BOOST_MATH_NAMESPACE::digamma(b, pol);
   }
   if (fabs(b) < max_factorial<T>::value)
   {
      // times sqrt(pi) = rh + rl; sqrt(pi) is 1.772453850905516 - 7.666586499825799e-17 - 1.3058334907945429e-33 in doubles
      typedef long double ld;
      const T rh = constants::root_pi<T>();
      const T rl = std::is_floating_point<T>::value && tools::digits<T>() <= std::numeric_limits<ld>::digits
         ? static_cast<T>((static_cast<ld>(1.772453850905516) - static_cast<ld>(rh)) - static_cast<ld>(7.666586499825799e-17) - static_cast<ld>(1.3058334907945429e-33))
         : T(0);
      const double_word_arithmetic<T> dw;
      const double_word<T> d = dw.multiply(gamma_parts(b, pol), double_word<T>{ rh, rl });
      T p = pow(x2, a);
      if ((BOOST_MATH_NAMESPACE::isfinite)(p) && p != 0 && (BOOST_MATH_NAMESPACE::isnormal)(d.hi))
      {
         const T q = p / d.hi;
         return q + q * (correction - d.lo / d.hi);
      }
      const T h = pow(x2, a / 2);
      return (h / (d.hi + d.lo)) * h * (1 + correction);
   }
   int s = 1;
   T lg = BOOST_MATH_NAMESPACE::lgamma(b, &s, pol);
   return s * exp(a * log(x2) - lg) * (1 + correction) / constants::root_pi<T>();
}

// Power series DLMF 11.2.1 (modified = false) and 11.2.2 (modified = true):
//   sum_k (-+1)^k (x/2)^(2k+nu+1) / (Gamma(k+3/2) Gamma(k+nu+3/2)),
// summed term by term with the x-derivative alongside. *abs_sum receives sum |t_k|.
template <class T, class Policy>
T power_series(bool modified, const T& nu, const T& x, T* derivative, T* abs_sum, const Policy& pol)
{
   BOOST_MATH_STD_USING
   const char* function = modified ? "boost::math::struve_l<%1%>(%1%, %1%)" : "boost::math::struve_h<%1%>(%1%, %1%)";
   const T x2 = x / 2;
   const T q = modified ? x2 * x2 : -x2 * x2;
   long k = 0;
   if (is_nonpositive_integer(T(nu + T(1.5))))
   {
      k = BOOST_MATH_NAMESPACE::lltrunc(-(nu + T(1.5)), pol) + 1;
   }
   // Gamma(k + 3/2) = (sqrt(pi)/2) prod_{i=1}^k (i + 1/2)
   T t = 2 * scale(nu, static_cast<T>(2 * k + 1), static_cast<T>(k) + T(1.5), x2, pol);
   for (long i = 1; i <= k; ++i)
   {
      t /= static_cast<T>(i) + T(0.5);
   }
   if (!modified && (k & 1))
   {
      t = -t;
   }
   if (!(BOOST_MATH_NAMESPACE::isfinite)(t))
   {
      // The first term overflows, and those after it are far smaller.
      const T v = (t < 0 ? -1 : 1) * policies::raise_overflow_error<T>(function, nullptr, pol);
      if (derivative)
      {
         *derivative = 2 * static_cast<T>(k) + nu + 1 < 0 ? T(-v) : v;
      }
      if (abs_sum)
      {
         *abs_sum = fabs(v);
      }
      return v;
   }
   T sum = 0;
   T dsum = 0;
   T asum = 0;
   const T eps = tools::epsilon<T>();
   for (long n = 0;; ++n, ++k)
   {
      if (n == 100000)
      {
         return policies::raise_evaluation_error<T>(function, "Power series failed to converge, x = %1%", x, pol);
      }
      const T kk = static_cast<T>(k);
      sum += t;
      dsum += t * (2 * kk + nu + 1);
      asum += fabs(t);
      const T b = kk + nu + T(1.5);
      const T ratio = q / ((kk + T(1.5)) * b);
      if (b > 0 && fabs(ratio) < T(0.5) && fabs(t) <= eps * fabs(sum) / 4)
      {
         break;
      }
      t *= ratio;
   }
   if (derivative)
   {
      *derivative = dsum / x;
   }
   if (abs_sum)
   {
      *abs_sum = asum;
   }
   return sum;
}

// H_nu by the Neumann series DLMF 11.4.18,
//   H_nu(x) = sqrt(x/(2 pi)) sum_k (x/2)^k / (k! (k + 1/2)) J_{k+nu+1/2}(x),
// for mu = nu + 1/2 >= 0. The J_{mu+k} come from Miller's backward recurrence, normalized by
//   (x/2)^mu = sum_k (mu + 2k) Gamma(mu + k) / k! J_{mu+2k}(x)   (DLMF 10.23.15).
// While x is below the first zero of J_mu every term of both sums is positive.
template <class T, class Policy>
T neumann_series(const T& nu, const T& x, T* derivative, const Policy& pol)
{
   BOOST_MATH_STD_USING
   const T mu = nu + T(0.5);
   const T x2 = x / 2;
   const T eps = tools::epsilon<T>();
   // Forward test recurrence (Olver): p_m grows like Y_{mu+m}/J_{mu+m}; start where p_m > a_m / eps.
   T p_prev = 0;
   T p = 1;
   T a = 2;
   long m = 0;
   while (fabs(p) * eps < a || m < 4)
   {
      T p_next = 2 * (mu + static_cast<T>(m + 1)) / x * p - p_prev;
      p_prev = p;
      p = p_next;
      ++m;
      const T mm = static_cast<T>(m);
      a *= x2 / mm * (mm - T(0.5)) / (mm + T(0.5));
      if (m > 100000)
      {
         return policies::raise_evaluation_error<T>("boost::math::struve_h<%1%>(%1%, %1%)", "Neumann series failed to converge, x = %1%", x, pol);
      }
   }
   // The normalization sum needs J_{mu+2j} for as many j as the H sum needs J_{mu+k} for k.
   const std::size_t m0 = 2 * static_cast<std::size_t>(m) + 16;
   std::vector<T> f(m0 + 2);
   f[m0 + 1] = 0;
   f[m0] = tools::min_value<T>() / eps;
   for (std::size_t j = m0; j > 0; --j)
   {
      f[j - 1] = 2 * (mu + static_cast<T>(j)) / x * f[j] - f[j + 1];
      if (fabs(f[j - 1]) > tools::max_value<T>() * eps)
      {
         for (std::size_t i = j - 1; i <= m0 + 1; ++i)
         {
            f[i] *= eps;
         }
      }
   }
   // Both sums relative to f[0], each term from the one before: the coefficients alone overflow for large mu.
   T norm = 1;
   T t = 1;
   for (std::size_t j = 1; 2 * j <= m0 && f[2 * j] != 0; ++j)
   {
      const T jj = static_cast<T>(j);
      t *= (j == 1 ? mu + 2 : (mu + 2 * jj) * (mu + jj - 1) / ((mu + 2 * jj - 2) * jj)) * (f[2 * j] / f[2 * j - 2]);
      norm += t;
   }
   T sum = 0;
   T dsum = 0;
   t = 2;
   for (std::size_t k = 0; k <= m0 && f[k] != 0; ++k)
   {
      const T kk = static_cast<T>(k);
      if (k > 0)
      {
         t *= x2 / kk * (kk - T(0.5)) / (kk + T(0.5)) * (f[k] / f[k - 1]);
      }
      sum += t;
      dsum += t * ((2 * kk + mu + T(0.5)) / x - f[k + 1] / f[k]);
   }
   const T pre = scale(nu, T(1), T(1.5), x2, pol);
   if (derivative)
   {
      *derivative = pre * (dsum / norm);
   }
   return pre * (sum / norm);
}

// K_nu(x) = H_nu(x) - Y_nu(x) by the large-x expansion DLMF 11.6.1,
//   K_nu(x) ~ (1/pi) sum_k Gamma(k+1/2) (x/2)^(nu-2k-1) / Gamma(nu+1/2-k),
// and M_nu(x) = L_nu(x) - I_nu(x) by DLMF 11.6.2 (alternating signs). Returns false unless a term falls
// below eps/2 times |sum| + scale_floor before the terms start to grow.
template <class T, class Policy>
bool asymptotic_series(bool modified, const T& nu, const T& x, T& value, T& deriv, const T& scale_floor, const Policy& pol)
{
   BOOST_MATH_STD_USING
   const T x2 = x / 2;
   const T inv = 1 / (x2 * x2);
   // First nonzero term: Gamma(nu + 1/2 - k) has a pole for every k when nu + 1/2 is a nonpositive integer.
   const T b0 = nu + T(0.5);
   if (is_nonpositive_integer(b0))
   {
      value = 0;
      deriv = 0;
      return true;
   }
   T t = scale(nu, T(-1), T(0.5), x2, pol);
   if (modified)
   {
      t = -t;
   }
   T sum = 0;
   T dsum = 0;
   T last = tools::max_value<T>();
   for (long k = 0; k < 100000; ++k)
   {
      if (fabs(t) >= last)
      {
         return false;
      }
      const T kk = static_cast<T>(k);
      sum += t;
      dsum += t * (nu - 2 * kk - 1);
      if (t == 0 || fabs(t) <= tools::epsilon<T>() * (fabs(sum) + scale_floor) / 2)
      {
         value = sum;
         deriv = dsum / x;
         return true;
      }
      last = fabs(t);
      t *= (kk + T(0.5)) * (nu - T(0.5) - kk) * inv;
      if (modified)
      {
         t = -t;
      }
   }
   return false;
}

// K_nu(x) = H_nu(x) - Y_nu(x) by DLMF 11.5.2 with s = x t:
//   K_nu(x) = 2 (x/2)^nu / (sqrt(pi) Gamma(nu + 1/2) x) int_0^inf e^-s (1 + s^2/x^2)^(nu - 1/2) ds,
// a positive integrand, by exp-sinh quadrature with compensated summation. DLMF states it for nu > -1/2, but both sides are
// entire in nu for x > 0 (the integral converges uniformly), so it holds for every real nu. The derivative is
//   K'_nu(x) = 2 (x/2)^nu / (sqrt(pi) Gamma(nu + 1/2) x^2) int_0^inf e^-s (nu - s) (1 + s^2/x^2)^(nu - 1/2) ds.
// The two are integrated together when the derivative is wanted, and both must converge.
template <class T, class Policy>
T k_integral(const T& nu, const T& x, T* derivative, const Policy& pol)
{
   BOOST_MATH_STD_USING
   // 9 refinements suffice for up to 64 bits, and each further one doubles the digits.
   static const quadrature::exp_sinh<T, Policy> integrator([] {
      std::size_t levels = 9;
      for (int d = tools::digits<T>(); d > 64; d /= 2)
      {
         ++levels;
      }
      return levels;
   }(), true);
   auto f = [&](const T& s) {
      const T r = s / x;
      return exp((nu - T(0.5)) * BOOST_MATH_NAMESPACE::log1p(r * r, pol) - s);
   };
   // Alone, the value needs a tighter tolerance to be as accurate as with the derivative's terms in the error estimate.
   const T tol = (derivative ? 4 : 2) * tools::epsilon<T>();
   T error;
   T l1;
   T r;
   const T pre = 2 * scale(nu, T(0), T(0.5), T(x / 2), pol) / x;
   if (derivative)
   {
      // The integral and the derivative's integral as the real and imaginary parts:
      auto g = [&](const T& s) {
         const T v = f(s);
         return std::complex<T>(v, (nu - s) * v);
      };
      auto norm = [](const std::complex<T>& z) { return T(fabs(z.real()) + fabs(z.imag())); };
      const std::complex<T> z = integrator.integrate(g, std::complex<T>(0), norm, tol, &error, &l1);
      r = z.real();
      *derivative = pre * z.imag() / x;
   }
   else
   {
      r = integrator.integrate(f, tol, &error, &l1);
   }
   if (error > tol * l1)
   {
      return policies::raise_evaluation_error<T>("boost::math::struve_h<%1%>(%1%, %1%)", "Integral failed to converge, x = %1%", x, pol);
   }
   return pre * r;
}

// Where the Neumann series gives way to K + Y: below the first zero of J_mu (j_{mu,1} > mu + 1.85 mu^(1/3)), where
// every term of the Neumann series is positive, and late enough that |Y| is not much larger than |H|.
template <class T>
T first_bessel_zero_bound(const T& mu)
{
   BOOST_MATH_STD_USING
   return mu < 1 ? T(2.4) : T(mu + cbrt(mu));
}

// The power series is used where its terms do not cancel by more than this factor relative to
// max(|f|, |x f'|), the scale of the rounding error that the condition number allows.
template <class T>
bool series_is_accurate(const T& value, const T& x_derivative, const T& abs_sum)
{
   BOOST_MATH_STD_USING
   return abs_sum <= 4 * (std::max)(fabs(value), fabs(x_derivative));
}

// x = 0: H_nu(0) = L_nu(0) = 0 for nu > -1, 2/pi for nu = -1, 0 at nu = -3/2, -5/2, ... (where the functions are
// +-J_{n+1/2} and I_{n+1/2}), and infinite otherwise, with the sign of 1/Gamma(nu + 3/2). The derivative is
// 2/pi at nu = 0, zero for nu > 0, nu = -1 and nu = -3/2, -5/2, ..., and infinite otherwise.
template <class T, class Policy>
T value_at_zero(const T& nu, T* derivative, const char* function, const Policy& pol)
{
   BOOST_MATH_STD_USING
   const bool pole = is_nonpositive_integer(T(nu + T(1.5)));
   int sign = 1;
   if (nu < -1 && !pole)
   {
      BOOST_MATH_NAMESPACE::lgamma(T(nu + T(1.5)), &sign, pol);
   }
   if (derivative)
   {
      if (nu > 0 || nu == -1 || pole)
      {
         *derivative = 0;
      }
      else if (nu == 0)
      {
         *derivative = constants::two_div_pi<T>();
      }
      else
      {
         // d/dx of c x^(nu+1), c of the sign of 1/Gamma(nu + 3/2)
         *derivative = (nu > -1 ? sign : -sign) * policies::raise_overflow_error<T>(function, nullptr, pol);
      }
   }
   if (nu > -1 || pole)
   {
      return 0;
   }
   if (nu == -1)
   {
      return constants::two_div_pi<T>();
   }
   return sign * policies::raise_overflow_error<T>(function, nullptr, pol);
}

// x > 0
template <class T, class Policy>
T struve_h_imp(const T& nu, const T& x, T* derivative, const Policy& pol)
{
   BOOST_MATH_STD_USING
   if (nu == T(0.5))
   {
      // H_{1/2}(x) = sqrt(2/(pi x)) (1 - cos x) (DLMF 11.4.5) has double zeros at x = 2 pi k, where H = K + Y
      // would cancel; 1 - cos x = 2 sin^2(x/2).
      const T s = sqrt(2 / (constants::pi<T>() * x));
      const T sh = sin(x / 2);
      const T v = 2 * s * sh * sh;
      if (derivative)
      {
         *derivative = 2 * s * sh * cos(x / 2) - v / (2 * x);
      }
      return v;
   }
   const T mu = nu + T(0.5);
   // Power series, where it does not cancel.
   if (x < 20 || x * x < 25 * (nu + T(1.5)))
   {
      T d, abs_sum;
      const T v = power_series(false, nu, x, &d, &abs_sum, pol);
      if (series_is_accurate(v, T(x * d), abs_sum))
      {
         if (derivative)
         {
            *derivative = d;
         }
         return v;
      }
   }
   // Below the first zero of J_{nu+1/2}: the Neumann series in J_{nu+1/2+k}, all of whose terms are positive.
   if (mu >= 0 && x < first_bessel_zero_bound(mu))
   {
      return neumann_series(nu, x, derivative, pol);
   }
   // H = K + Y.
   T k, kd;
   if (!asymptotic_series(false, nu, x, k, kd, T(0), pol))
   {
      k = k_integral(nu, x, derivative ? &kd : nullptr, pol);
   }
   if (derivative)
   {
      *derivative = kd + BOOST_MATH_NAMESPACE::cyl_neumann_prime(nu, x, pol);
   }
   return k + BOOST_MATH_NAMESPACE::cyl_neumann(nu, x, pol);
}

// x > 0
template <class T, class Policy>
T struve_l_imp(const T& nu, const T& x, T* derivative, const Policy& pol)
{
   BOOST_MATH_STD_USING
   // Large x: L = I_nu + M_nu, M_nu from its asymptotic expansion; M is far smaller than I there. The expansion
   // determines M only up to exponentially small terms (it is also that of L - I_{-nu}, and I_{-nu} - I_nu is a
   // multiple of K_nu), so it needs K_nu / I_nu ~ pi e^(-2x) below eps/16.
   // For nu < 0 those terms also include multiples of K_|nu|, and K_|nu| / I_|nu| ~ pi e^(-2g) (DLMF 10.41.3, 10.41.4).
   const T r = sqrt(nu * nu + x * x);
   const T g = nu < 0 ? T(r + nu * log((r - nu) / x)) : x;
   if (x >= nu && g >= log(16 / tools::epsilon<T>()) / 2)
   {
      const T i = BOOST_MATH_NAMESPACE::cyl_bessel_i(nu, x, pol);
      T m, md;
      if (asymptotic_series(true, nu, x, m, md, T(fabs(i) / 4), pol))
      {
         if (derivative)
         {
            *derivative = BOOST_MATH_NAMESPACE::cyl_bessel_i_prime(nu, x, pol) + md;
         }
         return i + m;
      }
   }
   // The power series. Its terms are all positive for nu > -3/2; for lower orders only the first few alternate,
   // and for -20 <= nu < -3/2 the errors measured against MPFR are within 2 ulps.
   return power_series(true, nu, x, derivative, static_cast<T*>(nullptr), pol);
}

template <class T, class Policy>
T struve_imp(bool modified, const T& nu, const T& x, T* derivative, const Policy& pol)
{
   BOOST_MATH_STD_USING
   const char* function = modified ? "boost::math::struve_l<%1%>(%1%, %1%)" : "boost::math::struve_h<%1%>(%1%, %1%)";
   if (!(BOOST_MATH_NAMESPACE::isfinite)(nu))
   {
      return policies::raise_domain_error<T>(function, "Order must be finite, but got nu = %1%", nu, pol);
   }
   if (!(BOOST_MATH_NAMESPACE::isfinite)(x))
   {
      return policies::raise_domain_error<T>(function, "Argument must be finite, but got x = %1%", x, pol);
   }
   if (x < 0)
   {
      // H_n(-x) = (-1)^(n+1) H_n(x), H'_n(-x) = (-1)^n H'_n(x), and likewise for L
      if (floor(nu) != nu)
      {
         return policies::raise_domain_error<T>(function, "Argument must be >= 0 for non-integer order, but got x = %1%", x, pol);
      }
      const bool odd = fmod(nu, T(2)) != 0;
      const T v = struve_imp(modified, nu, T(-x), derivative, pol);
      if (derivative && odd)
      {
         *derivative = -*derivative;
      }
      return odd ? v : T(-v);
   }
   if (x == 0)
   {
      return value_at_zero(nu, derivative, function, pol);
   }
   return modified ? struve_l_imp(nu, x, derivative, pol) : struve_h_imp(nu, x, derivative, pol);
}

template <class T1, class T2, class Policy>
inline typename tools::promote_args<T1, T2>::type struve_api(bool modified, T1 nu, T2 x,
   typename tools::promote_args<T1, T2>::type* derivative, const Policy&)
{
   typedef typename tools::promote_args<T1, T2>::type result_type;
   typedef typename policies::evaluation<result_type, Policy>::type value_type;
   typedef typename policies::normalise<
      Policy,
      policies::promote_float<false>,
      policies::promote_double<false>,
      policies::discrete_quantile<>,
      policies::assert_undefined<> >::type forwarding_policy;
   const char* function = modified ? "boost::math::struve_l<%1%>(%1%, %1%)" : "boost::math::struve_h<%1%>(%1%, %1%)";
   value_type d = 0;
   const value_type v = struve_imp(modified, static_cast<value_type>(nu), static_cast<value_type>(x), derivative ? &d : static_cast<value_type*>(nullptr), forwarding_policy());
   if (derivative)
   {
      *derivative = policies::checked_narrowing_cast<result_type, forwarding_policy>(d, function);
   }
   return policies::checked_narrowing_cast<result_type, forwarding_policy>(v, function);
}

}} // namespace detail::struve

// Struve function H_nu(x), DLMF 11.2.1, and optionally its derivative with respect to x.
BOOST_MATH_EXPORT template <class T1, class T2, class Policy>
inline typename tools::promote_args<T1, T2>::type struve_h(T1 nu, T2 x, typename tools::promote_args<T1, T2>::type* derivative, const Policy& pol)
{
   return detail::struve::struve_api(false, nu, x, derivative, pol);
}

BOOST_MATH_EXPORT template <class T1, class T2>
inline typename tools::promote_args<T1, T2>::type struve_h(T1 nu, T2 x, typename tools::promote_args<T1, T2>::type* derivative = nullptr)
{
   return detail::struve::struve_api(false, nu, x, derivative, policies::policy<>());
}

BOOST_MATH_EXPORT template <class T1, class T2, class Policy, BOOST_MATH_NAMESPACE::enable_if_t<policies::is_policy_v<Policy>, bool> = true>
inline typename tools::promote_args<T1, T2>::type struve_h(T1 nu, T2 x, const Policy& pol)
{
   return detail::struve::struve_api(false, nu, x, static_cast<typename tools::promote_args<T1, T2>::type*>(nullptr), pol);
}

// Modified Struve function L_nu(x), DLMF 11.2.2, and optionally its derivative with respect to x.
BOOST_MATH_EXPORT template <class T1, class T2, class Policy>
inline typename tools::promote_args<T1, T2>::type struve_l(T1 nu, T2 x, typename tools::promote_args<T1, T2>::type* derivative, const Policy& pol)
{
   return detail::struve::struve_api(true, nu, x, derivative, pol);
}

BOOST_MATH_EXPORT template <class T1, class T2>
inline typename tools::promote_args<T1, T2>::type struve_l(T1 nu, T2 x, typename tools::promote_args<T1, T2>::type* derivative = nullptr)
{
   return detail::struve::struve_api(true, nu, x, derivative, policies::policy<>());
}

BOOST_MATH_EXPORT template <class T1, class T2, class Policy, BOOST_MATH_NAMESPACE::enable_if_t<policies::is_policy_v<Policy>, bool> = true>
inline typename tools::promote_args<T1, T2>::type struve_l(T1 nu, T2 x, const Policy& pol)
{
   return detail::struve::struve_api(true, nu, x, static_cast<typename tools::promote_args<T1, T2>::type*>(nullptr), pol);
}

BOOST_MATH_NAMESPACE_END

#endif // BOOST_MATH_SPECIAL_STRUVE_HPP
