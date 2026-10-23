
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

   Real area() const
   {
      return BOOST_MATH_NAMESPACE::tgamma(alpha + 1);
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

private:
   unsigned N;
   Real alpha;
};
}

BOOST_MATH_EXPORT template <class T, class Policy = policies::policy<>>
inline std::vector<T> laguerre_zeros(unsigned n)
{
   static const char* function = "boost::math::laguerre_zeros(unsigned %1%)";
   if ((!(BOOST_MATH_NAMESPACE::isfinite)(n)))
   {
      policies::raise_domain_error(function, "The argument n to the laguerre_zeros function must be finite (got n=%1%).", n, Policy());
      return { std::numeric_limits<T>::quiet_NaN() };
   }
   if (n == 0)
   {
      policies::raise_domain_error(function, "Laguerre polynomial has no roots for n = 0.", n, Policy());
      return { std::numeric_limits<T>::quiet_NaN() };
   }
   if (n == 1)
   {
      return { static_cast<T>(1) };
   }

   detail::orthogonal_polynomial<T, detail::laguerre_family<T> > evaluate(n);
   std::vector<T> roots = evaluate.abscissa();

   if (!(BOOST_MATH_NAMESPACE::isfinite)(roots[0]))
   {
      // In this case, the recurrence has overflowed and returned NaNs. Throw evaluation error
      // or return vector of nans.
      policies::raise_evaluation_error(function, "End points of the interval is not finite. Probably due to n being too large.", n, Policy());
      return roots;
   }

   return roots;
}

BOOST_MATH_NAMESPACE_END

#endif // BOOST_MATH_SPECIAL_LAGUERRE_HPP



