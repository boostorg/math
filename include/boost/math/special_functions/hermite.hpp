
//  (C) Copyright John Maddock 2006.
//  (C) Copyright Matt Borland 2024.
//  (C) Copyright Jacob Hass 2026.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_SPECIAL_HERMITE_HPP
#define BOOST_MATH_SPECIAL_HERMITE_HPP

#ifdef _MSC_VER
#pragma once
#endif

#ifndef BOOST_MATH_BUILD_MODULE
#include <vector>
#endif
#include <boost/math/tools/config.hpp>
#include <boost/math/tools/promotion.hpp>
#include <boost/math/special_functions/math_fwd.hpp>
#include <boost/math/policies/error_handling.hpp>
#include <boost/math/special_functions/detail/orthogonal_polynomial.hpp>

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
      p0_value = 1 / sqrt(sqrt(boost::math::constants::pi<Real>()));
      derivative_scale = sqrt(Real(2 * N));
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
      return (N+1) / 2;
   }

   Real lower_bound() const
   {
      return 0;
   }

   Real upper_bound() const
   {
      return sqrt(Real(2 * N + 2));
   }

   Real next(const Real& x, const Real& p, const Real& p_previous, unsigned k) const
   {
      return a(k) * x * p + b(k) * p_previous;
   }

   Real derivative(const Real& x, const Real& p, const Real& p_previous) const
   {
      return derivative_scale * p_previous;
   }

   Real second_derivative(const Real& x, const Real& p, const Real& p_prime) const
   {
      return 2 * x * p_prime - 2 * Real(N) * p;
   }

   Real second_derivative_ratio(const Real& x, const Real& p, const Real& p_prime) const
   {
      return 2 * x - 2 * Real(N) * p / p_prime;
   }

   Real weight(const Real& z, const Real& p, const Real& p_prime) const
   {
      Real delta = -p / p_prime;
      Real p_prime_at_root = p_prime + second_derivative(z, p, p_prime) * delta;
      return 2 / (p_prime_at_root * p_prime_at_root);
   }

private:
   unsigned N;
   std::vector<Real> a_values;
   std::vector<Real> b_values;
   Real p0_value;
   Real derivative_scale;
};

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

// Todo: add policy for nan handling
BOOST_MATH_EXPORT template <class T>
inline std::vector<T> hermite_zeros(unsigned n)
{
   detail::orthogonal_polynomial<T, detail::hermite_family<T> > evaluate(n);
   std::vector<T> roots = evaluate.abscissa();
   return roots;
}

BOOST_MATH_NAMESPACE_END

#endif // BOOST_MATH_SPECIAL_HERMITE_HPP



