
//  (C) Copyright John Maddock 2006.
//  (C) Copyright Matt Borland 2024.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_SPECIAL_HERMITE_HPP
#define BOOST_MATH_SPECIAL_HERMITE_HPP

#ifdef _MSC_VER
#pragma once
#endif

#include <boost/math/tools/config.hpp>
#include <boost/math/tools/promotion.hpp>
#include <boost/math/special_functions/math_fwd.hpp>
#include <boost/math/policies/error_handling.hpp>
#ifndef BOOST_MATH_HAS_GPU_SUPPORT
#ifndef BOOST_MATH_BUILD_MODULE
#include <vector>
#endif
#include <boost/math/special_functions/detail/orthogonal_polynomial.hpp>
#endif

BOOST_MATH_NAMESPACE_BEGIN

// Recurrence relation for Hermite polynomials:
BOOST_MATH_EXPORT template <class T1, class T2, class T3>
BOOST_MATH_GPU_ENABLED inline typename tools::promote_args<T1, T2, T3>::type 
   hermite_next(unsigned n, T1 x, T2 Hn, T3 Hnm1)
{
   using promoted_type = tools::promote_args_t<T1, T2, T3>;
   return (2 * promoted_type(x) * promoted_type(Hn) - 2 * n * promoted_type(Hnm1));
}

namespace detail{

// Implement Hermite polynomials via recurrence:
template <class T>
BOOST_MATH_GPU_ENABLED T hermite_imp(unsigned n, T x)
{
   T p0 = 1;
   T p1 = 2 * x;

   if(n == 0)
      return p0;

   unsigned c = 1;

   while(c < n)
   {
      BOOST_MATH_GPU_SAFE_SWAP(p0, p1);
      p1 = static_cast<T>(hermite_next(c, x, p0, p1));
      ++c;
   }
   return p1;
}

} // namespace detail

BOOST_MATH_EXPORT template <class T, class Policy>
BOOST_MATH_GPU_ENABLED inline typename tools::promote_args<T>::type 
   hermite(unsigned n, T x, const Policy&)
{
   typedef typename tools::promote_args<T>::type result_type;
   typedef typename policies::evaluation<result_type, Policy>::type value_type;
   return policies::checked_narrowing_cast<result_type, Policy>(detail::hermite_imp(n, static_cast<value_type>(x)), "boost::math::hermite<%1%>(unsigned, %1%)");
}

BOOST_MATH_EXPORT template <class T>
BOOST_MATH_GPU_ENABLED inline typename tools::promote_args<T>::type 
   hermite(unsigned n, T x)
{
   return BOOST_MATH_NAMESPACE::hermite(n, x, policies::policy<>());
}

#ifndef BOOST_MATH_HAS_GPU_SUPPORT

namespace detail {

template <class Real>
class hermite_family
{
public:
   hermite_family(unsigned n) : N(n), a_values(N), b_values(N)
   {
      using std::sqrt;
      for (unsigned k = 0; k < N; ++k)
      {
         a_values[k] = sqrt(Real(2) / Real(k + 1));
         b_values[k] = -sqrt(Real(k) / Real(k + 1));
      }
      p0_value = 1 / sqrt(sqrt(BOOST_MATH_NAMESPACE::constants::pi<Real>()));
      derivative_scale = sqrt(Real(2) * Real(N));
   }

   const Real& a(unsigned k) const
   {
      return a_values[k];
   }

   const Real& b(unsigned k) const
   {
      return b_values[k];
   }

   const Real& p0() const
   {
      return p0_value;
   }

   unsigned num_roots() const
   {  // Need to account for odd N having root at 0
      return N / 2 + N % 2;
   }

   Real lower_bound() const
   {
      return 0;
   }

   Real upper_bound() const
   {
      using std::sqrt;
      return sqrt(Real(2) * Real(N) + Real(2));
   }

   Real area() const
   {
      return BOOST_MATH_NAMESPACE::constants::root_pi<Real>();
   }

   unsigned weight_multiplicity(unsigned i) const
   {
      return (N % 2 && i == 0) ? 1 : 2;
   }

   Real next(const Real& x, const Real& p, const Real& p_previous, unsigned k) const
   {
      return a(k) * x * p + b(k) * p_previous;
   }

   Real slope(unsigned k) const
   {  // The coefficient of x p in next()
      return a(k);
   }

   Real derivative(const Real&, const Real&, const Real& p_previous) const
   {
      return derivative_scale * p_previous;
   }

   Real second_derivative(const Real& x, const Real& p, const Real& p_prime) const
   {
      return Real(2) * x * p_prime - Real(2) * Real(N) * p;
   }

   Real second_derivative_ratio(const Real& x, const Real& p, const Real& p_prime) const
   {
      return Real(2) * x - Real(2) * Real(N) * p / p_prime;
   }

   Real weight(const Real& z, const Real& p, const Real& p_prime) const
   {
      Real delta = -p / p_prime;
      Real p_prime_at_root = p_prime + second_derivative(z, p, p_prime) * delta;
      return Real(2) / (p_prime_at_root * p_prime_at_root);
   }

   // For glaser_liu_rokhlin: u = exp(-x^2 / 2) p_N satisfies u'' + (2N + 1 - x^2) u = 0; the coefficients of
   // s(t + h) = 1 and g(t + h) as polynomials in h. Its solutions are entire.
   static constexpr unsigned ends = 0;

   void normal_form(const double_word_arithmetic<Real>& dw, const double_word<Real>& t, double_word<Real> (&s)[5], double_word<Real> (&g)[3], unsigned m) const
   {
      s[0] = double_word<Real>{ Real(1), Real(0) };
      g[0] = dw.add(dw.negate(dw.multiply(t, t)), Real(2) * Real(m) + 1);
      g[1] = double_word<Real>{ -2 * t.hi, -2 * t.lo };
      g[2] = double_word<Real>{ Real(-1), Real(0) };
   }

   Real singular_distance(const Real&) const
   {
      return tools::max_value<Real>();
   }

   Real weight_function(const double_word<Real>& z) const
   {
      using std::exp;
      const double_word<Real> square = double_word_arithmetic<Real>().multiply(z, z);
      return exp(-square.hi) * (1 - square.lo);
   }

   Real sturm_liouville(const double_word<Real>&) const
   {
      return 1;
   }

   // p_N = (A x + B) p_{N-1} - C p_{N-2}: A = sqrt(2 / N), B = 0, and C_2 = -2 C = -2 sqrt((N - 1) / N), in double words.
   void anti_gauss_recurrence(const double_word_arithmetic<Real>& dw, double_word<Real>& A, double_word<Real>& B, double_word<Real>& C2) const
   {
      A = square_root(dw, dw.divide(double_word<Real>{ Real(2), Real(0) }, Real(N)));
      B = double_word<Real>{ Real(0), Real(0) };
      const double_word<Real> C = square_root(dw, dw.divide(double_word<Real>{ Real(N - 1), Real(0) }, Real(N)));
      C2 = double_word<Real>{ -2 * C.hi, -2 * C.lo };
   }

   // One step of Newton's method for the square root, in double words.
   static double_word<Real> square_root(const double_word_arithmetic<Real>& dw, const double_word<Real>& x)
   {
      using std::sqrt;
      const Real root = sqrt(x.hi);
      Real p, e;
      two_product<Real>()(root, root, p, e);
      return dw.add(double_word<Real>{ root, Real(0) }, (((x.hi - p) - e) + x.lo) / (2 * root));
   }

   // For the anti-Gauss polynomial: p_{N-1}(0) / p_{N-2}'(0) = -1 / sqrt(2 (N - 1)) or p_{N-1}'(0) / p_{N-2}(0) = sqrt(2 (N - 1)),
   // from p_k' = sqrt(2 k) p_{k-1} and p_{N-1}(0) = b(N - 2) p_{N-3}(0).
   double_word<Real> anti_gauss_ratio(const double_word_arithmetic<Real>& dw, unsigned) const
   {
      const double_word<Real> r = square_root(dw, double_word<Real>{ 2 * Real(N - 1), Real(0) });
      return (N - 1) % 2 ? r : dw.divide(double_word<Real>{ Real(-1), Real(0) }, r);
   }

private:
   unsigned N;
   std::vector<Real> a_values;
   std::vector<Real> b_values;
   Real p0_value;
   Real derivative_scale;
};

} // namespace detail

// The nonnegative zeros of H_n, in increasing order; the others are their negatives.
BOOST_MATH_EXPORT template <class T, class Policy>
inline std::vector<T> hermite_zeros(unsigned n, const Policy& pol)
{
   using value_type = typename policies::evaluation<T, Policy>::type;
   return detail::orthogonal_polynomial_zeros<T, value_type>(n, detail::hermite_family<value_type>(n), "boost::math::hermite_zeros<%1%>(unsigned)", pol);
}

BOOST_MATH_EXPORT template <class T>
inline std::vector<T> hermite_zeros(unsigned n)
{
   return BOOST_MATH_NAMESPACE::hermite_zeros<T>(n, policies::policy<>());
}

#endif // BOOST_MATH_HAS_GPU_SUPPORT

BOOST_MATH_NAMESPACE_END

#endif // BOOST_MATH_SPECIAL_HERMITE_HPP



