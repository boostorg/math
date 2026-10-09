
//  (C) Copyright John Maddock 2006.
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

// Implement Laguerre polynomials via recurrence:
template <class T>
T laguerre_imp(unsigned n, T x)
{
   T p0 = 1;
   T p1 = 1 - x;

   if(n == 0)
      return p0;

   unsigned c = 1;

   while(c < n)
   {
      std::swap(p0, p1);
      p1 = laguerre_next(c, x, p0, p1);
      ++c;
   }
   return p1;
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

template <class Real>
class laguerre_family
{
public:
   laguerre_family(unsigned n) : N(n) {}

   Real p0() const
   {
      return 1;
   }

   Real lower_bound() const
   {
      return 0;
   }

   Real upper_bound() const
   { // See https://mathoverflow.net/questions/251607/zeroes-of-laguerre-polynomials
      return N + (N-1) * sqrt(N) + 1;
   }

   unsigned num_roots() const
   {
      return N;
   }

   Real next(const Real& x, const Real& p, const Real& p_previous, unsigned k) const
   {
      return ((2*k + 1 - x) * p - k * p_previous) / (k + 1);
   }

   Real derivative(const Real& x, const Real& p, const Real& p_previous) const
   {
      return N / x * (p - p_previous);
   }

   Real second_derivative(const Real& x, const Real& p, const Real& p_prime) const
   {
      return ((x-1) * p_prime - N * p) / x;;
   }

   Real second_derivative_ratio(const Real& x, const Real& p, const Real& p_prime) const
   {
      return (x-1) / x - N / x * p / p_prime;
   }

   Real weight(const Real& x, const Real& p, const Real& p_prime) const
   {
      Real delta = -p / p_prime;
      Real p_prime_at_root = p_prime + second_derivative(x, p, p_prime) * delta;
      return 1 / (p_prime_at_root * p_prime_at_root) / x;
   }

private:
   unsigned N;
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



