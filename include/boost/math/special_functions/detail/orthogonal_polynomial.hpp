//  Copyright Jacob Hass 2026.
//  Copyright Nick Thompson 2026.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_DETAIL_ORTHOGONAL_POLYNOMIAL_HPP
#define BOOST_MATH_DETAIL_ORTHOGONAL_POLYNOMIAL_HPP

#ifndef BOOST_MATH_BUILD_MODULE
#include <cmath>
#include <vector>
#include <limits>
#include <stdexcept>
#include <algorithm>
#include <utility>
#endif
#include <boost/math/constants/constants.hpp>
#include <boost/math/special_functions/fpclassify.hpp>

namespace boost { namespace math { namespace detail {

template <class Real, class Family>
class orthogonal_polynomial
{
public:
orthogonal_polynomial(unsigned n) : N(n)
{
   if (N == 0)
   {
      throw std::domain_error("Gauss quadrature needs at least one point.");
   }

   values = calculate_values();
}
const std::vector<Real>& abscissa() const
{
   return values.first;
}
const std::vector<Real>& weights() const
{
   return values.second;
}

private:
unsigned N;
std::pair<std::vector<Real>, std::vector<Real> > values;
struct recurrence
{
   Real p;
   Real p_prime;
   unsigned sign_changes;
};

class orthonormal_hermite
{
public:
   orthonormal_hermite(unsigned n) : N(n), family(N) {}

   recurrence operator()(const Real& x) const
   {
      Real p = family.p0();
      Real p_previous = 0;
      bool last_negative = false;
      unsigned sign_changes = 0;
      for (unsigned k = 0; k < N; ++k)
      {
         Real p_next = family.next(x, p, p_previous, k);
         p_previous = p;
         p = p_next;
         // A zero p_k(x) lies between values of opposite sign, so skipping it keeps the count right.
         if ((p != 0) && ((p < 0) != last_negative))
         {
            ++sign_changes;
            last_negative = !last_negative;
         }
      }
      return recurrence{ p, family.derivative(x, p, p_previous), sign_changes };
   }

   Real second_derivative_ratio(const Real& x, const Real& p, const Real& p_prime) const
   {
      return family.second_derivative_ratio(x, p, p_prime);
   }

   Real weight(const Real& z, const Real& p, const Real& p_prime) const
   {
      return family.weight(z, p, p_prime);
   }

   unsigned num_roots() const
   {
      return family.num_roots();
   }

   Real lower_bound() const
   {
      return family.lower_bound();
   }

   Real upper_bound() const
   {
      return family.upper_bound();
   }

private:
   unsigned N;
   Family family;
};

// The nonnegative zeros of the orthogonal polynomial H_N, in increasing order, and their weights.
// Everything is evaluated with the orthonormal recurrence
// whose values stay representable long after H_N itself overflows. Evaluating the polynomial from its monomial
// coefficients instead would be hopelessly ill-conditioned.
//
// By Sturm's theorem for orthogonal polynomials, the number of sign changes in p_0(x), ..., p_N(x)
// changes by one at each zero of p_N. Its direction depends on the family. Bisecting on that count
// isolates each zero, and Halley's method then converges to it inside its bracket. Asymptotic initial guesses alone are not enough:
// for N = 215, even Halley's method from those of Numerical Recipes' gauher finds the wrong zero.
std::pair<std::vector<Real>, std::vector<Real> > calculate_values()
{
   using std::abs;
   using std::sqrt;
   const orthonormal_hermite evaluate(N);

   const unsigned num_roots = evaluate.num_roots();
   unsigned interior_roots = num_roots;
   std::vector<Real> x(num_roots), w(num_roots);
   // All zeros lie below upper_bound(), so none lie above `upper`. Every p_k grows with x beyond
   // its zeros, so if the recurrence overflows anywhere we need it, it overflows here: report
   // that with NaNs, which gauss_hermite turns into an evaluation error.
   Real upper = evaluate.upper_bound();
   Real lower = evaluate.lower_bound();
   recurrence top = evaluate(upper);
   if (!(boost::math::isfinite)(top.p) || !(boost::math::isfinite)(top.p_prime))
   {
      std::fill(x.begin(), x.end(), std::numeric_limits<Real>::quiet_NaN());
      std::fill(w.begin(), w.end(), std::numeric_limits<Real>::quiet_NaN());
      return std::make_pair(x, w);
   }
   // Check if lower bound is 0. If so, we have an odd number of roots, and the middle one is at 0.
   const recurrence bottom = evaluate(lower);
   if (bottom.p == 0)
   {
      x[0] = lower;
      w[0] = evaluate.weight(lower, bottom.p, bottom.p_prime);
      interior_roots--;
   }
   const bool count_increases = top.sign_changes > bottom.sign_changes;
   // Find the zeros from the largest down, following the variation-count direction of this family.
   for (unsigned t = 1; t <= interior_roots; ++t)
   {
      Real low = lower;
      Real high = upper;
      const unsigned target = count_increases
         ? top.sign_changes - t + 1
         : top.sign_changes + t;
      for (int i = 0; i < tools::digits<Real>(); ++i)
      {
         Real mid = (low + high) / 2;
         unsigned variations = evaluate(mid).sign_changes;
         if (count_increases)
         {
            if (variations >= target)
               high = mid;
            else
               low = mid;
         }
         else
         {
            if (variations >= target)
               low = mid;
            else
               high = mid;
         }
      }
      // Halley's method, falling back to bisection whenever a step would leave the bracket.
      // With u = p/p' and v = p''/p' = 2z - 2N u, the step is u/(1 - uv/2); the ratios cannot
      // overflow even where p'^2 would.
      const bool lower_negative = evaluate(low).p < 0;
      Real z = (low + high) / 2;
      recurrence r = evaluate(z);
      for (int i = 0; (r.p != 0) && (i < 4 * tools::digits<Real>()); ++i)
      {
         if ((r.p < 0) == lower_negative)
            low = z;
         else
            high = z;
         Real u = r.p / r.p_prime;
         Real v = evaluate.second_derivative_ratio(z, r.p, r.p_prime);
         Real next = z - u / (1 - u * v / 2);
         if (!((next > low) && (next < high)))
            next = (low + high) / 2;
         const bool converged = abs(next - z) <= 2 * tools::epsilon<Real>() * abs(next);
         z = next;
         r = evaluate(z);
         if (converged)
            break;
      }
      x[num_roots - t] = z;
      w[num_roots - t] = evaluate.weight(z, r.p, r.p_prime);
      upper = low;
   }
   return std::make_pair(x, w);
}

};

} // namespace detail
} // namespace math
} // namespace boost

#endif // BOOST_MATH_DETAIL_ORTHOGONAL_POLYNOMIAL_HPP
