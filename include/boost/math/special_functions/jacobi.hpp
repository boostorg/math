//  (C) Copyright Nick Thompson 2019.
//  (C) Copyright Jacob Hass 2026.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_SPECIAL_JACOBI_HPP
#define BOOST_MATH_SPECIAL_JACOBI_HPP

#ifndef BOOST_MATH_BUILD_MODULE
#include <cmath>
#include <limits>
#include <stdexcept>
#endif

#include <boost/math/tools/config.hpp>
#include <boost/math/special_functions/beta.hpp>
#include <boost/math/special_functions/detail/orthogonal_polynomial.hpp>

BOOST_MATH_NAMESPACE_BEGIN

BOOST_MATH_EXPORT template<typename Real>
Real jacobi(unsigned n, Real alpha, Real beta, Real x)
{
    static_assert(!std::is_integral<Real>::value, "Jacobi polynomials do not work with integer arguments.");

    if (n == 0) {
        return Real(1);
    }
    Real y0 = 1;
    Real y1 = (alpha+1) + (alpha+beta+2)*(x-1)/Real(2);

    Real yk = y1;
    Real k = 2;
    Real k_max = n*(1+std::numeric_limits<Real>::epsilon());
    while(k < k_max)
    {
        // Hoping for lots of common subexpression elimination by the compiler:
        Real denom = 2*k*(k+alpha+beta)*(2*k+alpha+beta-2);
        Real gamma1 = (2*k+alpha+beta-1)*( (2*k+alpha+beta)*(2*k+alpha+beta-2)*x + alpha*alpha -beta*beta);
        Real gamma0 = -2*(k+alpha-1)*(k+beta-1)*(2*k+alpha+beta);
        yk = (gamma1*y1 + gamma0*y0)/denom;
        y0 = y1;
        y1 = yk;
        k += 1;
    }
    return yk;
}

BOOST_MATH_EXPORT template<typename Real>
Real jacobi_derivative(unsigned n, Real alpha, Real beta, Real x, unsigned k)
{
    if (k > n) {
        return Real(0);
    }
    Real scale = 1;
    for(unsigned j = 1; j <= k; ++j) {
        scale *= (alpha + beta + n + j)/2;
    }

    return scale*jacobi<Real>(n-k, alpha + k, beta+k, x);
}

BOOST_MATH_EXPORT template<typename Real>
Real jacobi_prime(unsigned n, Real alpha, Real beta, Real x)
{
    return jacobi_derivative<Real>(n, alpha, beta, x, 1);
}

BOOST_MATH_EXPORT template<typename Real>
Real jacobi_double_prime(unsigned n, Real alpha, Real beta, Real x)
{
    return jacobi_derivative<Real>(n, alpha, beta, x, 2);
}

namespace detail {

template <class Real>
class jacobi_family
{
public:
   jacobi_family(unsigned n, Real alpha_, Real beta_) : N(n), alpha(alpha_), beta(beta_) {}

   Real p0() const
   {
      return 1;
   }

   Real lower_bound() const
   {
      return -1;
   }

   Real upper_bound() const
   {
      return 1;
   }

   unsigned num_roots() const
   {
      return N;
   }

   Real area() const
   {
      using std::pow;
      return pow(Real(2), alpha + beta + 1) * BOOST_MATH_NAMESPACE::beta(alpha + 1, beta + 1);
   }

   unsigned weight_multiplicity(unsigned) const
   {
      return 1;
   }

   Real next(const Real& x, const Real& p, const Real& p_previous, unsigned k) const
   {  // See https://dlmf.nist.gov/18.9
      if (k == 0)
      {
         return ((alpha - beta) + (alpha + beta + 2) * x) / 2;
      }

      const Real degree = static_cast<Real>(k);
      const Real sum = alpha + beta;
      const Real twice_degree_sum = 2 * degree + sum;
      const Real first = (twice_degree_sum + 1) *
         (twice_degree_sum * (twice_degree_sum + 2) * x + alpha * alpha - beta * beta) * p;
      const Real second = 2 * (degree + alpha) * (degree + beta) * (twice_degree_sum + 2) * p_previous;
      const Real denominator = 2 * (degree + 1) * (degree + sum + 1) * twice_degree_sum;
      return (first - second) / denominator;
   }

   Real derivative(const Real& x, const Real& p, const Real& p_previous) const
   {  // See https://dlmf.nist.gov/18.9#E17
      const Real degree = static_cast<Real>(N);
      const Real sum = alpha + beta;
      if (x == 1)
      {
         return p * degree * (degree + sum + 1) / (2 * (alpha + 1));
      }
      if (x == -1)
      {
         return -p * degree * (degree + sum + 1) / (2 * (beta + 1));
      }
      return (degree * (alpha - beta - (2 * degree + sum) * x) * p +
         2 * (degree + alpha) * (degree + beta) * p_previous) / (1 - x * x);
   }

   Real second_derivative(const Real& x, const Real& p, const Real& p_prime) const
   {
      const Real degree = static_cast<Real>(N);
      const Real sum = alpha + beta;
      return ((alpha - beta + (sum + 2) * x) * p_prime -
         degree * (degree + sum + 1) * p) / (1 - x * x);
   }

   Real second_derivative_ratio(const Real& x, const Real& p, const Real& p_prime) const
   {
      const Real degree = static_cast<Real>(N);
      const Real sum = alpha + beta;
      return (alpha - beta + (sum + 2) * x - degree * (degree + sum + 1) * p / p_prime) / (1 - x * x);
   }

   Real weight(const Real& x, const Real& p, const Real& p_prime) const
   {
      Real delta = -p / p_prime;
      Real p_prime_at_root = p_prime + second_derivative(x, p, p_prime) * delta;
      return 1 / ((1 - x * x) * p_prime_at_root * p_prime_at_root);
   }

private:
   unsigned N;
   Real alpha;
   Real beta;
};
}

BOOST_MATH_EXPORT template <class T, class Policy = boost::math::policies::policy<> >
inline std::vector<T> jacobi_zeros(unsigned n, T alpha, T beta)
{
   static const char* function = "boost::math::jacobi_zeros(unsigned %1%)";
   if ((!(BOOST_MATH_NAMESPACE::isfinite)(n)))
   {
      policies::raise_domain_error(function, "The argument n to the jacobi_zeros function must be finite (got n=%1%).", n, Policy());
      return { std::numeric_limits<T>::quiet_NaN() };
   }
   if ((alpha <= -1) || (beta <=-1))
   {
      policies::raise_domain_error(function, "Parameters alpha and beta must be > -1.", n, Policy());
      return { std::numeric_limits<T>::quiet_NaN() };
   }
   if ((!(BOOST_MATH_NAMESPACE::isfinite)(alpha)) || (!(BOOST_MATH_NAMESPACE::isfinite)(beta)))
   {
      policies::raise_domain_error(function, "Parameters alpha and beta must be finite.", n, Policy());
      return { std::numeric_limits<T>::quiet_NaN() };
   }

   if (n == 0)
   {
      policies::raise_domain_error(function, "Jacob polynomial has no roots for n = 0.", n, Policy());
      return { std::numeric_limits<T>::quiet_NaN() };
   }
   if (n == 1)
   {
      return { static_cast<T>(1) - static_cast<T>(2) * (alpha + static_cast<T>(1)) / (alpha + beta + static_cast<T>(2)) };
   }
   detail::jacobi_family<T> family(n, alpha, beta);
   detail::orthogonal_polynomial<T, detail::jacobi_family<T> > evaluate(n, family);
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
#endif
