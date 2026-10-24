
//  (C) Copyright John Maddock 2006.
//  (C) Copyright Jacob Hass 2026.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_SPECIAL_LAGUERRE_HPP
#define BOOST_MATH_SPECIAL_LAGUERRE_HPP

#ifdef _MSC_VER
#pragma once
#endif

#include <boost/math/special_functions/math_fwd.hpp>
#include <boost/math/tools/config.hpp>
#include <boost/math/policies/error_handling.hpp>
#include <boost/math/special_functions/detail/orthogonal_polynomial.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/special_functions/digamma.hpp>
#include <boost/math/tools/compensated_arithmetic.hpp>

BOOST_MATH_NAMESPACE_BEGIN

// Recurrence relation for Laguerre polynomials:
BOOST_MATH_EXPORT template <class T1, class T2, class T3>
inline typename tools::promote_args<T1, T2, T3>::type
   laguerre_next(unsigned n, T1 x, T2 Ln, T3 Lnm1)
{
   typedef typename tools::promote_args<T1, T2, T3>::type result_type;
   return ((2 * n + 1 - result_type(x)) * result_type(Ln) - n * result_type(Lnm1)) / (n + 1);
}

namespace detail{

//
// Implement Laguerre polynomials via recurrence.  On its own the recurrence loses up to n^1.5 epsilon or so
// where L_n oscillates, as both of its solutions are of similar size there and the rounding errors of each step
// accumulate.  So we find the rounding error of each step almost exactly, and carry it forward with the same
// recurrence: the result is then within the conditioning of L_n.
//
template <class T>
T laguerre_imp(unsigned n, T x)
{
   T p0 = 1;
   T p1 = 1 - x;

   if(n == 0)
      return p0;

   // Errors in p0 and p1, starting with that of 1 - x:
   T e0 = 0;
   T e1 = two_sum_error(T(1), T(-x), p1);
   detail::two_product<T> product;

   for(unsigned c = 1; c < n; ++c)
   {
      // a + ea = 2c + 1 - x exactly:
      T a = (2 * c + 1) - x;
      T ea = two_sum_error(T(2 * c + 1), T(-x), a);
      // s1 + r1 = a p1 and s2 + r2 = c p0 exactly:
      T s1, r1, s2, r2;
      product(a, p1, s1, r1);
      product(T(c), p0, s2, r2);
      // d + ed = s1 - s2 exactly:
      T d = s1 - s2;
      T ed = two_sum_error(s1, T(-s2), d);
      T p2 = d / (c + 1);
      // d - (c + 1) p2, where d - q is exact as q is so close to d:
      T q, rq;
      product(T(c + 1), p2, q, rq);
      T rd = (d - q) - rq;
      T e2 = (a * e1 - c * e0 + rd + ed + r1 - r2 + ea * p1) / (c + 1);
      p0 = p1;
      p1 = p2;
      e0 = e1;
      e1 = e2;
   }
   return p1 + e1;
}

template <class T, class Policy>
inline typename tools::promote_args<T>::type
laguerre(unsigned n, T x, const Policy&, const std::true_type&)
{
   typedef typename tools::promote_args<T>::type result_type;
   typedef typename policies::evaluation<result_type, Policy>::type value_type;
   return policies::checked_narrowing_cast<result_type, Policy>(detail::laguerre_imp(n, static_cast<value_type>(x)), "boost::math::laguerre<%1%>(unsigned, %1%)");
}

template <class T>
inline typename tools::promote_args<T>::type
   laguerre(unsigned n, unsigned m, T x, const std::false_type&)
{
   return BOOST_MATH_NAMESPACE::laguerre(n, m, x, policies::policy<>());
}

} // namespace detail

BOOST_MATH_EXPORT template <class T>
inline typename tools::promote_args<T>::type
   laguerre(unsigned n, T x)
{
   return laguerre(n, x, policies::policy<>());
}

// Recurrence for associated polynomials:
BOOST_MATH_EXPORT template <class T1, class T2, class T3>
inline typename tools::promote_args<T1, T2, T3>::type
   laguerre_next(unsigned n, unsigned l, T1 x, T2 Pl, T3 Plm1)
{
   typedef typename tools::promote_args<T1, T2, T3>::type result_type;
   return ((2 * n + l + 1 - result_type(x)) * result_type(Pl) - (n + l) * result_type(Plm1)) / (n+1);
}

namespace detail{
// Laguerre Associated Polynomial:
template <class T, class Policy>
T laguerre_imp(unsigned n, unsigned m, T x, const Policy& pol)
{
   // Special cases:
   if(m == 0)
      return BOOST_MATH_NAMESPACE::laguerre(n, x, pol);

   T p0 = 1;

   if(n == 0)
      return p0;

   T p1 = m + 1 - x;

   unsigned c = 1;

   while(c < n)
   {
      std::swap(p0, p1);
      p1 = static_cast<T>(laguerre_next(c, m, x, p0, p1));
      ++c;
   }
   return p1;
}

}

BOOST_MATH_EXPORT template <class T, class Policy>
inline typename tools::promote_args<T>::type
   laguerre(unsigned n, unsigned m, T x, const Policy& pol)
{
   typedef typename tools::promote_args<T>::type result_type;
   typedef typename policies::evaluation<result_type, Policy>::type value_type;
   return policies::checked_narrowing_cast<result_type, Policy>(detail::laguerre_imp(n, m, static_cast<value_type>(x), pol), "boost::math::laguerre<%1%>(unsigned, unsigned, %1%)");
}

BOOST_MATH_EXPORT template <class T1, class T2>
inline typename laguerre_result<T1, T2>::type
   laguerre(unsigned n, T1 m, T2 x)
{
   typedef typename policies::is_policy<T2>::type tag_type;
   return detail::laguerre(n, m, x, tag_type());
}

namespace detail{

// The generalized Laguerre polynomials L_n^(alpha), orthogonal for x^alpha exp(-x) on [0, infinity).
template <class Real>
class laguerre_family
{
public:
   laguerre_family(unsigned n, Real alpha_ = 0) : N(n), alpha(alpha_) {}

   Real p0() const
   {
      return 1;
   }

   Real lower_bound() const
   {
      return 0;
   }

   Real upper_bound() const
   {  // Gershgorin's theorem for the Jacobi matrix, diagonal 2k + alpha + 1 and off-diagonal sqrt(k (k + alpha)),
      // bounds the zeros by 4N + 2 alpha; (3 + sqrt 2) N + 3 alpha + 3 also covers the anti-Gauss rule.
      using std::sqrt;
      return (3 + sqrt(Real(2))) * static_cast<Real>(N) + 3 * alpha + 3;
   }

   unsigned num_roots() const
   {
      return N;
   }

   // Gamma(alpha + 1), with the rounding error e of alpha + 1 = s + e to first order: Gamma(s) (1 + psi(s) e).
   Real area() const
   {
      const Real s = alpha + 1;
      return BOOST_MATH_NAMESPACE::tgamma(s) * (1 + BOOST_MATH_NAMESPACE::digamma(s) * BOOST_MATH_NAMESPACE::detail::two_sum_error(alpha, Real(1), s));
   }

   unsigned weight_multiplicity(unsigned) const
   {
      return 1;
   }

   Real next(const Real& x, const Real& p, const Real& p_previous, unsigned k) const
   {
      const Real degree = static_cast<Real>(k);
      // x enters as x p, not through 2k + 1 + alpha - x, which would round away the low bits of the small zeros.
      return ((Real(2) * degree + 1 + alpha) * p - (degree + alpha) * p_previous - x * p) / (degree + 1);
   }

   // p_n(x), or Laurie's p~_n(x) with the coefficient of p_{n-2} doubled in the last step, to about twice the
   // working precision: each step carries the rounding errors of its products, sums and quotient forward (a
   // compensated recurrence). Near the small zeros the terms of the recurrence are about n times p_n, so in
   // working precision those zeros, and the weights from p_{n-1} there, are only good to about n eps / x.
   Real compensated_value(const Real& x, unsigned n, bool anti, Real& previous) const
   {
      const BOOST_MATH_NAMESPACE::detail::two_product<Real> product;
      Real p = p0(), e = 0, p_previous = 0, e_previous = 0;
      for (unsigned k = 0; k < n; ++k)
      {
         const Real c = (anti && (k + 1 == n)) ? Real(2) : Real(1);
         const Real d = static_cast<Real>(k + 1);
         const Real a = static_cast<Real>(2 * k + 1) + alpha;
         const Real ea = BOOST_MATH_NAMESPACE::detail::two_sum_error<Real>(static_cast<Real>(2 * k + 1), alpha, a);
         Real b = static_cast<Real>(k) + alpha;
         const Real eb = c * BOOST_MATH_NAMESPACE::detail::two_sum_error<Real>(static_cast<Real>(k), alpha, b);
         b *= c;
         Real ap, eap, bp, ebp, xp, exp;
         product(a, p, ap, eap);
         product(b, p_previous, bp, ebp);
         product(x, p, xp, exp);
         const Real s1 = ap - bp;
         const Real s = s1 - xp;
         const Real err = BOOST_MATH_NAMESPACE::detail::two_sum_error<Real>(ap, -bp, s1) + BOOST_MATH_NAMESPACE::detail::two_sum_error<Real>(s1, -xp, s)
            + eap - ebp - exp + (a - x) * e + ea * p - b * e_previous - eb * p_previous;
         const Real q = s / d;
         Real qd, eqd;
         product(q, d, qd, eqd);
         p_previous = p;
         e_previous = e;
         p = q;
         e = (((s - qd) - eqd) + err) / d;
      }
      previous = p_previous + e_previous;
      return p + e;
   }

   Real slope(unsigned k) const
   {  // The coefficient of x p in next()
      return -1 / static_cast<Real>(k + 1);
   }

   Real derivative(const Real& x, const Real& p, const Real& p_previous) const
   {
      return (static_cast<Real>(N) * p - (static_cast<Real>(N) + alpha) * p_previous) / x;
   }

   Real second_derivative(const Real& x, const Real& p, const Real& p_prime) const
   {
      return ((x - alpha - 1) * p_prime - static_cast<Real>(N) * p) / x;
   }

   Real second_derivative_ratio(const Real& x, const Real& p, const Real& p_prime) const
   {
      return (x - alpha - 1) / x - static_cast<Real>(N) / x * p / p_prime;
   }

   Real weight(const Real& x, const Real& p, const Real& p_prime) const
   {  // Up to the factor Gamma(N + alpha + 1) / N!, which the scaling to area() supplies.
      Real delta = -p / p_prime;
      Real p_prime_at_root = p_prime + second_derivative(x, p, p_prime) * delta;
      return Real(1) / (p_prime_at_root * p_prime_at_root) / x;
   }

   // For glaser_liu_rokhlin: u = x^((alpha + 1) / 2) exp(-x / 2) p_N satisfies x^2 u'' + g(x) u = 0 with
   // g = (1 - alpha^2) / 4 + (N + (alpha + 1) / 2) x - x^2 / 4; the coefficients of s(t + h) and g(t + h) in h.
   static constexpr unsigned ends = 1;

   void normal_form(const double_word_arithmetic<Real>& dw, const double_word<Real>& t, double_word<Real> (&s)[5], double_word<Real> (&g)[3], unsigned m) const
   {
      using word = double_word<Real>;
      const word c = dw.add(word{ Real(m), Real(0) }, dw.divide(dw.two_sum(alpha, Real(1)), Real(2)));
      const word a = dw.divide(dw.add(dw.negate(dw.multiply(word{ alpha, Real(0) }, alpha)), Real(1)), Real(4));
      const word square = dw.multiply(t, t);
      s[0] = square;
      s[1] = word{ 2 * t.hi, 2 * t.lo };
      s[2] = word{ Real(1), Real(0) };
      g[0] = dw.add(dw.add(word{ -square.hi / 4, -square.lo / 4 }, dw.multiply(c, t)), a);
      g[1] = dw.add(c, word{ -t.hi / 2, -t.lo / 2 });
      g[2] = word{ Real(-0.25), Real(0) };
   }

   Real singular_distance(const Real& x) const
   {
      return x;
   }

   // x^alpha exp(-x) at z = hi + lo; past where exp(-x) alone underflows, the square of its square root.
   Real weight_function(const double_word<Real>& z) const
   {
      using std::exp;
      using std::pow;
      const Real correction = 1 + (alpha / z.hi - 1) * z.lo;
      if (z.hi < -tools::log_min_value<Real>() / 2)
         return pow(z.hi, alpha) * exp(-z.hi) * correction;
      const Real root = pow(z.hi, alpha / 2) * exp(-z.hi / 2);
      return root * root * correction;
   }

   // L_N^alpha(lambda w) = L_N^alpha(0) sum_k c_k w^k with c_0 = 1, lambda a power of 2 near 1 / (4N), so that the
   // smallest zero is at w of order 1 (Gatteschi: x_1 ~ j_{alpha,1}^2 / (4N + 2 alpha + 2)).
   Real series_origin(unsigned) const
   {
      return 0;
   }

   Real series_step(unsigned) const
   {
      using std::frexp;
      using std::ldexp;
      int exponent;
      frexp(4 * Real(N) + 2 * alpha + 2, &exponent);
      return ldexp(Real(1), -exponent);
   }

   double_word<Real> series_ratio(const double_word_arithmetic<Real>& dw, unsigned k, unsigned end, unsigned m) const
   {  // c_{k+1} / c_k = lambda (k - m) / ((k + 1) (k + 1 + alpha)) for L_m^alpha
      const double_word<Real> denominator = dw.multiply(dw.two_sum(Real(k + 1), alpha), Real(k + 1));
      return dw.divide(double_word<Real>{ (Real(k) - Real(m)) * series_step(end), Real(0) }, denominator);
   }

   // L_N = (A x + B) L_{N-1} - C L_{N-2} with A = -1 / N, B = (2N - 1 + alpha) / N, C = (N - 1 + alpha) / N; C_2 = -2 C.
   void anti_gauss_recurrence(const double_word_arithmetic<Real>& dw, double_word<Real>& A, double_word<Real>& B, double_word<Real>& C2) const
   {
      A = dw.divide(double_word<Real>{ Real(-1), Real(0) }, Real(N));
      B = dw.divide(dw.two_sum(2 * Real(N) - 1, alpha), Real(N));
      C2 = dw.divide(dw.two_sum(-2 * Real(N - 1), -2 * alpha), Real(N));
   }

   // L_{N-1}^alpha(0) / L_{N-2}^alpha(0) = (N - 1 + alpha) / (N - 1)
   double_word<Real> anti_gauss_ratio(const double_word_arithmetic<Real>& dw, unsigned) const
   {
      return dw.divide(dw.two_sum(Real(N - 1), alpha), Real(N - 1));
   }

   // sigma' / sigma = (alpha + 1) / (2 x) - 1 / 2
   double_word<Real> scale_log_derivative(const double_word_arithmetic<Real>& dw, const double_word<Real>& x) const
   {
      return dw.add(dw.divide(dw.two_sum(alpha, Real(1)), double_word<Real>{ 2 * x.hi, 2 * x.lo }), Real(-0.5));
   }

   // The leading coefficient x of the Sturm-Liouville form (x^(alpha + 1) exp(-x) p')' = -N x^alpha exp(-x) p.
   Real sturm_liouville(const double_word<Real>& z) const
   {
      return z.hi + z.lo;
   }

private:
   unsigned N;
   Real alpha;
};
}

BOOST_MATH_EXPORT template <class T, class Policy = policies::policy<>>
inline std::vector<T> laguerre_zeros(unsigned n)
{
   static const char* function = "boost::math::laguerre_zeros(unsigned %1%)";
   if (n == 0)
   {
      policies::raise_domain_error(function, "Laguerre polynomial has no roots for n = 0.", n, Policy());
      return { std::numeric_limits<T>::quiet_NaN() };
   }
   if (n == 1)
   {
      return { static_cast<T>(1) };
   }
   return detail::polynomial_zeros<T>(n, detail::laguerre_family<T>(n), function, Policy());
}

BOOST_MATH_NAMESPACE_END

#endif // BOOST_MATH_SPECIAL_LAGUERRE_HPP



