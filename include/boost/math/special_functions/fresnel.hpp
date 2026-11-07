//  (C) Copyright Nick Thompson 2026.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_SPECIAL_FRESNEL_HPP
#define BOOST_MATH_SPECIAL_FRESNEL_HPP

// The Fresnel integrals C(x) = int_0^x cos(pi t^2/2) dt and S(x) = int_0^x sin(pi t^2/2) dt (DLMF 7.2.7, 7.2.8).

#include <boost/math/tools/config.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/math/tools/promotion.hpp>
#include <boost/math/tools/precision.hpp>
#include <boost/math/tools/compensated_arithmetic.hpp>
#include <boost/math/tools/detail/complex_from_real.hpp>
#include <boost/math/policies/error_handling.hpp>
#include <boost/math/special_functions/math_fwd.hpp>
#include <boost/math/special_functions/fpclassify.hpp>
#include <boost/math/special_functions/sin_pi.hpp>
#include <boost/math/special_functions/cos_pi.hpp>
#ifndef BOOST_MATH_BUILD_MODULE
#include <cmath>
#include <limits>
#include <type_traits>
#include <vector>
#endif

BOOST_MATH_NAMESPACE_BEGIN

namespace detail { namespace fresnel {

enum which { want_c = 1, want_s = 2, want_both = 3 };

// D, the precision in bits, is fixed at compile time (D > 0), so that what depends on it can be computed once; for
// variable-precision types (D = 0) it is recomputed at the current precision on every call.
template <class T, class Policy>
using precision_key = std::integral_constant<int,
   (std::numeric_limits<T>::is_specialized && std::numeric_limits<T>::radix == 2 && std::numeric_limits<T>::digits > 0
    && std::numeric_limits<T>::digits != (std::numeric_limits<int>::max)()) ? policies::precision<T, Policy>::type::value : 0>;

// pi - hi = sin(hi) for hi = pi rounded, so pi = hi + sin(hi) to twice working precision.
template <class T>
T make_pi_lo()
{
   BOOST_MATH_STD_USING
   return sin(constants::pi<T>());
}

template <class T, int D>
T pi_lo(std::integral_constant<int, D>)
{
   static const T lo = make_pi_lo<T>();
   return lo;
}

template <class T>
T pi_lo(std::integral_constant<int, 0>)
{
   return make_pi_lo<T>();
}

// pi x^2 in double-word arithmetic, x^2 exactly.
template <class T, class Key>
double_word<T> pi_x_squared(const T& x, Key key)
{
   double_word<T> x2;
   two_product<T>()(x, x, x2.hi, x2.lo);
   return double_word_arithmetic<T>().multiply(x2, double_word<T>{ constants::pi<T>(), pi_lo<T>(key) });
}

// The coefficients a_k = i^k / (k! (2k + 1)) of C + iS = x sum_k a_k theta^k, with the sign of i^k, in double-word
// arithmetic: as many as theta <= theta1 = digits / 10 can need, and one to spare.
template <class T>
std::vector<double_word<T> > make_series_coefficients(int digits)
{
   BOOST_MATH_STD_USING
   const double_word_arithmetic<T> dw;
   const T theta1 = T(digits) / 10;
   const T eps = ldexp(T(1), 1 - digits);
   std::vector<double_word<T> > v;
   double_word<T> u = { T(1), T(0) };   // 1 / k!
   T r = 1;                             // theta1^k / k!
   for (unsigned k = 0; k <= theta1 + 1 || r > eps * eps; ++k)
   {
      if (k > 0)
      {
         u = dw.divide(u, T(k));
         r *= theta1 / k;
      }
      const double_word<T> t = dw.divide(u, T(2 * k + 1));
      v.push_back(k % 4 > 1 ? double_word_arithmetic<T>::negate(t) : t);
   }
   return v;
}

template <class T, int D>
const std::vector<double_word<T> >& series_coefficients(int digits, std::integral_constant<int, D>)
{
   static const std::vector<double_word<T> > v = make_series_coefficients<T>(digits);
   return v;
}

template <class T>
std::vector<double_word<T> > series_coefficients(int digits, std::integral_constant<int, 0>)
{
   return make_series_coefficients<T>(digits);
}

// C and S for theta = pi x^2 / 2 <= digits / 10 from the Maclaurin series, by Horner's rule in theta^2 in double-word
// arithmetic: the terms add up to about e^theta / (2 theta) times the result, and C' = cos(theta) amplifies an error
// in theta by x^2 / 2.
template <class T, class Key>
void series(const T& x, T& c, T& s, int digits, const T& eps, Key key)
{
   typedef double_word_arithmetic<T> dw_type;
   typedef typename dw_type::value_type dw;
   const dw_type a;
   const std::vector<dw>& coef = series_coefficients<T>(digits, key);
   const dw pix2 = pi_x_squared(x, key);
   const dw theta = { T(pix2.hi / 2), T(pix2.lo / 2) };
   const dw theta2 = a.multiply(theta, theta);
   // The last term: x theta^n / (n! (2n + 1)) below eps / 16 of C and S, which exceed min(x theta / 4, 1/4).
   const T tol = eps / 16 * (std::min)(T(theta.hi / 4), T(T(0.25) / x));
   std::size_t n = 1;
   for (T r = theta.hi; n + 2 < coef.size() && (n <= theta.hi || r / (2 * n + 1) > tol); ++n)
      r *= theta.hi / (n + 1);
   n += n % 2 == 0;   // odd, so that C and S have (n + 1) / 2 terms each
   dw sc = coef[n - 1], ss = coef[n];
   for (std::size_t k = n; k >= 3; k -= 2)
   {
      sc = dw_type::add(a.multiply(sc, theta2), coef[k - 3]);
      ss = dw_type::add(a.multiply(ss, theta2), coef[k - 2]);
   }
   c = a.multiply(sc, x).hi;
   s = a.multiply(a.multiply(ss, theta), x).hi;
}

// g + i f = x / t_0 for theta > digits / 10 (DLMF 7.5.3, 7.5.4), with t_0 = 1 - 2i theta - 1*2/(5 - 2i theta - 3*4/(9 - ...)),
// the even part of Laplace's continued fraction for erfc((1 - i) sqrt(pi) x / 2) (DLMF 7.9.2), evaluated backwards.
// Measured, n terms give the full precision once n > 0.1 digits^2 / theta, and n > digits / (1.6 log2(theta)) for
// large theta.
template <class T, class Key>
void auxiliary(const T& x, T& f, T& g, int digits, const T& eps, Key key)
{
   BOOST_MATH_STD_USING
   if (x * x > 1 / eps)
   {
      // f = (1 - 3/(pi x^2)^2 + ...)/(pi x), g = (1 - 15/(pi x^2)^2 + ...)/(pi^2 x^3): the corrections are below eps^2.
      f = 1 / (constants::pi<T>() * x);
      g = constants::pi<T>() * f * f * f;
      return;
   }
   typedef double_word_arithmetic<T> dw_type;
   typedef typename dw_type::value_type dw;
   const dw_type a;
   const dw beta = dw_type::negate(pi_x_squared(x, key));   // -2 theta, the imaginary part of every b_n = 4n + 1 + i beta
   const T p = T(digits);
   const T n_max = ceil((std::max)(T(-T(0.22) * p * p / beta.hi),
                                   T(p * constants::ln_two<T>() / (T(1.6) * log(-beta.hi / 2))))) + 2;
   T tr = 4 * n_max + 1, ti = beta.hi;
   for (T n = n_max; n > 1; --n)
   {
      // t = b_{n-1} + a_n / t, a_n = -(2n - 1)(2n)
      const T q = -(2 * n - 1) * (2 * n) / (tr * tr + ti * ti);
      tr = (4 * n - 3) + q * tr;
      ti = beta.hi - q * ti;
   }
   // Nearly all the rounding error is in the last step, through beta, and in x / t_0: those in double-word arithmetic.
   const T q = -2 / (tr * tr + ti * ti);
   const T t0r = 1 + q * tr;
   const dw t0i = dw_type::fast_two_sum(beta.hi, T(beta.lo - q * ti));
   T d, de;
   two_product<T>()(t0i.hi, t0i.hi, d, de);
   const dw m = a.divide(dw{ x, T(0) }, dw_type::fast_two_sum(d, T(de + (2 * t0i.hi * t0i.lo + t0r * t0r))));
   g = a.multiply(m, t0r).hi;
   f = -a.multiply(m, t0i).hi;
}

// sin(pi x^2 / 2) and cos(pi x^2 / 2), with x^2 / 2 = h + l exactly, h reduced mod 2 exactly and l added by the
// angle-addition formulas: rounding h + l would cost up to pi/2 ulps.
template <class T, class Policy>
void phase(const T& x, T& sn, T& cs, const Policy& pol)
{
   BOOST_MATH_STD_USING
   if (x >= 2 / tools::epsilon<T>())
   {
      // x is an even integer, so x^2 / 2 is too (and x^2 may overflow).
      sn = 0;
      cs = 1;
      return;
   }
   T p, e;
   two_product<T>()(x, x, p, e);
   T h = fmod(T(p / 2), T(2));
   if (h > 1)
      h -= 2;
   const T l = e / 2;
   const T sh = BOOST_MATH_NAMESPACE::sin_pi(h, pol);
   const T ch = BOOST_MATH_NAMESPACE::cos_pi(h, pol);
   // For |l| < sqrt(eps) / 8, sin(pi l) = pi l and cos(pi l) = 1 to within eps / 6.
   const bool small = fabs(l) < tools::root_epsilon<T>() / 8;
   const T sl = small ? T(constants::pi<T>() * l) : T(BOOST_MATH_NAMESPACE::sin_pi(l, pol));
   const T cl = small ? T(1) : T(BOOST_MATH_NAMESPACE::cos_pi(l, pol));
   sn = sh * cl + ch * sl;
   cs = ch * cl - sh * sl;
}

// 1/2 + a b - c d, a b and c d exactly, the sums in double-word arithmetic: C and S fall to 1/3 near x = 2.
template <class T>
T half_plus(const T& a, const T& b, const T& c, const T& d)
{
   typedef double_word_arithmetic<T> dw_type;
   typedef typename dw_type::value_type dw;
   const two_product<T> product;
   dw ab, cd;
   product(a, b, ab.hi, ab.lo);
   product(c, d, cd.hi, cd.lo);
   return dw_type::add(dw_type::add(dw{ T(0.5), T(0) }, ab), dw_type::negate(cd)).hi;
}

template <class T, class Policy>
void fresnel_imp(const T& x, T& c, T& s, which w, const Policy& pol)
{
   BOOST_MATH_STD_USING
   if (x < 0)
   {
      fresnel_imp(T(-x), c, s, w, pol);
      c = -c;
      s = -s;
      return;
   }
   if ((BOOST_MATH_NAMESPACE::isinf)(x))
   {
      c = s = T(0.5);
      return;
   }
   const precision_key<T, Policy> key;
   const int digits = policies::digits<T, Policy>();
   const T eps = policies::get_epsilon<T, Policy>();
   if (constants::half_pi<T>() * x * x <= T(digits) / 10)
   {
      series(x, c, s, digits, eps, key);
      return;
   }
   T f, g, sn, cs;
   auxiliary(x, f, g, digits, eps, key);
   phase(x, sn, cs, pol);
   if (w & want_c)
      c = half_plus(f, sn, g, cs);
   if (w & want_s)
      s = half_plus(T(-f), cs, g, sn);
}

// C(x) and S(x) in the result type, of which only those in w are computed and narrowed.
template <class T, class Policy>
void fresnel_api(T x, typename tools::promote_args<T>::type& c, typename tools::promote_args<T>::type& s, which w,
                 const char* function, const Policy&)
{
   typedef typename tools::promote_args<T>::type result_type;
   typedef typename policies::evaluation<result_type, Policy>::type value_type;
   typedef typename policies::normalise<
      Policy,
      policies::promote_float<false>,
      policies::promote_double<false>,
      policies::discrete_quantile<>,
      policies::assert_undefined<> >::type forwarding_policy;
   const value_type z = static_cast<value_type>(x);
   if ((BOOST_MATH_NAMESPACE::isnan)(z))
   {
      c = s = policies::raise_domain_error<result_type>(function, "Argument x must not be NaN, but got %1%", static_cast<result_type>(x), forwarding_policy());
      return;
   }
   value_type cv = 0, sv = 0;
   fresnel_imp(z, cv, sv, w, forwarding_policy());
   if (w & want_c)
      c = policies::checked_narrowing_cast<result_type, forwarding_policy>(cv, function);
   if (w & want_s)
      s = policies::checked_narrowing_cast<result_type, forwarding_policy>(sv, function);
}

}} // namespace detail::fresnel

// C(x) = int_0^x cos(pi t^2/2) dt, DLMF 7.2.7.
BOOST_MATH_EXPORT template <class T, class Policy>
inline typename tools::promote_args<T>::type fresnel_c(T x, const Policy& pol)
{
   typename tools::promote_args<T>::type c = 0, s = 0;
   detail::fresnel::fresnel_api(x, c, s, detail::fresnel::want_c, "boost::math::fresnel_c<%1%>(%1%)", pol);
   return c;
}

BOOST_MATH_EXPORT template <class T>
inline typename tools::promote_args<T>::type fresnel_c(T x)
{
   return fresnel_c(x, policies::policy<>());
}

// S(x) = int_0^x sin(pi t^2/2) dt, DLMF 7.2.8.
BOOST_MATH_EXPORT template <class T, class Policy>
inline typename tools::promote_args<T>::type fresnel_s(T x, const Policy& pol)
{
   typename tools::promote_args<T>::type c = 0, s = 0;
   detail::fresnel::fresnel_api(x, c, s, detail::fresnel::want_s, "boost::math::fresnel_s<%1%>(%1%)", pol);
   return s;
}

BOOST_MATH_EXPORT template <class T>
inline typename tools::promote_args<T>::type fresnel_s(T x)
{
   return fresnel_s(x, policies::policy<>());
}

// C(x) + i S(x) = int_0^x cis(pi t^2/2) dt, in std::complex, or for multiprecision types the complex type of
// tools::detail::complex_from_real (whose header the caller includes).
BOOST_MATH_EXPORT template <class T, class Policy>
inline tools::detail::complex_from_real_t<typename tools::promote_args<T>::type> fresnel_cis(T x, const Policy& pol)
{
   typedef tools::detail::complex_from_real_t<typename tools::promote_args<T>::type> result_type;
   typename tools::promote_args<T>::type c = 0, s = 0;
   detail::fresnel::fresnel_api(x, c, s, detail::fresnel::want_both, "boost::math::fresnel_cis<%1%>(%1%)", pol);
   return result_type(c, s);
}

BOOST_MATH_EXPORT template <class T>
inline tools::detail::complex_from_real_t<typename tools::promote_args<T>::type> fresnel_cis(T x)
{
   return fresnel_cis(x, policies::policy<>());
}

BOOST_MATH_NAMESPACE_END

#endif // BOOST_MATH_SPECIAL_FRESNEL_HPP
