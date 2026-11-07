//  (C) Copyright Nick Thompson 2026.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// Newton's method and Cuyt's multivariate Halley method for a system F(x) = 0, x in R^n, globalized by a
// backtracking line search on |F|^2/2.
//
// Halley: with the Newton step a = -J^-1 F and b = J^-1 F''[a, a], x_i += a_i^2 / (a_i + b_i/2).
// A. Cuyt and P. Van der Cruyssen, Abstract Pade approximants for the solution of a system of nonlinear
//   equations, Comput. Math. Appl. 9 (1983) 617-624.
// A. Cuyt and L. B. Rall, Computational implementation of the multivariate Halley method for solving
//   nonlinear systems of equations, ACM TOMS 11 (1985) 20-36.

#ifndef BOOST_MATH_TOOLS_MULTIVARIATE_ROOTS_HPP
#define BOOST_MATH_TOOLS_MULTIVARIATE_ROOTS_HPP

#include <boost/math/tools/config.hpp>
#include <boost/math/policies/policy.hpp>
#include <boost/math/policies/error_handling.hpp>
#include <boost/math/special_functions/fpclassify.hpp>
#include <boost/math/tools/precision.hpp>
#ifndef BOOST_MATH_BUILD_MODULE
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>
#endif

BOOST_MATH_NAMESPACE_BEGIN
namespace tools {
namespace detail {

// J[i][j] for a container of rows, J[i*n + j] for a row-major array.
template <class M, class = void>
struct is_nested_matrix : std::false_type {};

template <class M>
struct is_nested_matrix<M, decltype(void(std::declval<M const&>()[0][0]))> : std::true_type {};

template <class Real, class M>
bool copy_jacobian(M const& J, std::vector<Real>& lu, std::size_t n, std::true_type)
{
   if (J.size() != n)
      return false;
   for (std::size_t i = 0; i < n; ++i)
   {
      if (J[i].size() != n)
         return false;
      for (std::size_t j = 0; j < n; ++j)
         lu[i * n + j] = J[i][j];
   }
   return true;
}

template <class Real, class M>
bool copy_jacobian(M const& J, std::vector<Real>& lu, std::size_t n, std::false_type)
{
   if (J.size() != n * n)
      return false;
   for (std::size_t k = 0; k < n * n; ++k)
      lu[k] = J[k];
   return true;
}

// LU with partial pivoting of the row-major n x n matrix a, in place. False when a pivot is not above n eps times
// the largest entry of its column, so the test does not depend on the scaling of the variables.
template <class Real>
bool lu_factor(std::vector<Real>& a, std::vector<std::size_t>& piv, std::vector<Real>& colmax, std::size_t n)
{
   using std::abs;
   for (std::size_t j = 0; j < n; ++j)
   {
      colmax[j] = 0;
      for (std::size_t i = 0; i < n; ++i)
         colmax[j] = (std::max)(colmax[j], Real(abs(a[i * n + j])));
   }
   Real const tiny = n * tools::epsilon<Real>();
   for (std::size_t k = 0; k < n; ++k)
   {
      std::size_t p = k;
      for (std::size_t i = k + 1; i < n; ++i)
         if (abs(a[i * n + k]) > abs(a[p * n + k]))
            p = i;
      piv[k] = p;
      if (p != k)
         for (std::size_t j = 0; j < n; ++j)
            std::swap(a[k * n + j], a[p * n + j]);
      Real const d = a[k * n + k];
      if (!(abs(d) > tiny * colmax[k]) || !(BOOST_MATH_NAMESPACE::isfinite)(d))
         return false;
      for (std::size_t i = k + 1; i < n; ++i)
      {
         Real const l = a[i * n + k] /= d;
         for (std::size_t j = k + 1; j < n; ++j)
            a[i * n + j] -= l * a[k * n + j];
      }
   }
   return true;
}

template <class Real, class V>
void lu_solve(std::vector<Real> const& a, std::vector<std::size_t> const& piv, std::size_t n, V& b)
{
   for (std::size_t k = 0; k < n; ++k)
      std::swap(b[k], b[piv[k]]);
   for (std::size_t i = 1; i < n; ++i)
      for (std::size_t j = 0; j < i; ++j)
         b[i] -= a[i * n + j] * b[j];
   for (std::size_t i = n; i-- > 0;)
   {
      for (std::size_t j = i + 1; j < n; ++j)
         b[i] -= a[i * n + j] * b[j];
      b[i] /= a[i * n + i];
   }
}

template <class V>
auto inf_norm(V const& v) -> typename std::decay<decltype(v[0])>::type
{
   using std::abs;
   typename std::decay<decltype(v[0])>::type r = 0;
   for (std::size_t i = 0; i < v.size(); ++i)
      r = (std::max)(r, decltype(r)(abs(v[i])));
   return r;
}

// The Halley component a^2 / (a + b/2) = a / q, q = 1 + b/(2a). Cuyt and Rall (1985), eq. (7.2), take 0/0 as 0.
// Where q < 1/2 the Pade step points against the Newton step or is more than twice as long, and the iteration is
// far from the cubic regime, where q -> 1: take the Newton component a instead.
template <class Real>
Real halley_component(Real const& a, Real const& b)
{
   if (a == 0)
      return 0;
   Real const q = 1 + b / (2 * a);
   return q >= Real(0.5) ? Real(a / q) : a;
}

// b = F''[a, a], from the closure f returned as its third element; nothing for Newton. False if it has the wrong size.
template <class R, class V>
bool second_derivative(std::true_type, R& r, V const& a, V& b)
{
   auto const& w = std::get<2>(r)(a);
   if (w.size() != a.size())
      return false;
   for (std::size_t i = 0; i < a.size(); ++i)
      b[i] = w[i];
   return true;
}

template <class R, class V>
bool second_derivative(std::false_type, R&, V const&, V&)
{
   return true;
}

// x with every component the value of the evaluation error, NaN unless the policy's handler says otherwise.
template <class Container>
Container fail(Container x, char const* function, char const* message)
{
   using Real = typename std::decay<decltype(x[0])>::type;
   Real const v = policies::raise_evaluation_error<Real>(function, message, std::numeric_limits<Real>::quiet_NaN(),
                                                         policies::policy<>());
   for (std::size_t i = 0; i < x.size(); ++i)
      x[i] = v;
   return x;
}

template <bool Halley, class F, class Container>
Container multivariate_root_finder(F& f, Container x, int digits, std::uintmax_t& max_iter, char const* function)
{
   using std::sqrt;
   using std::ldexp;
   using Real = typename std::decay<decltype(x[0])>::type;
   using R = typename std::decay<decltype(f(x))>::type;
   using nested = is_nested_matrix<typename std::decay<decltype(std::get<1>(std::declval<R&>()))>::type>;
   std::size_t const n = x.size();
   std::uintmax_t const limit = max_iter;
   max_iter = 0;
   if (n == 0)
   {
      policies::raise_domain_error<Real>(function, "x must not be empty", Real(0), policies::policy<>());
      return x;
   }
   if (digits > tools::digits<Real>() || digits < 1)
      digits = tools::digits<Real>();
   Real const tol = ldexp(Real(1), 1 - digits);
   Real const eps = tools::epsilon<Real>();
   using std::cbrt;
   // A Newton step this small that cannot decrease |F| is rounding noise; see the line search below.
   Real const noise = cbrt(eps);
   Real const armijo = Real(1) / 10000;

   std::vector<Real> lu(n * n), colmax(n), g(n);
   std::vector<std::size_t> piv(n);
   Container fx(x), a(x), b(x), xt(x);
   // F into fx, J into lu, |F|^2/2 into phi, or infinity if F has a non-finite entry, and whether F = 0 exactly into
   // zero (phi can underflow). False if F or J has the wrong shape.
   // Every trial point is unpacked, so fx and lu hold the accepted point's values only once it is accepted; within an
   // iteration only a, b and g, already computed, are used.
   bool zero = false;
   auto unpack = [&](R const& r, Real& phi) {
      auto const& Fv = std::get<0>(r);
      if (Fv.size() != n || !copy_jacobian<Real>(std::get<1>(r), lu, n, nested()))
         return false;
      phi = 0;
      zero = true;
      for (std::size_t i = 0; i < n; ++i)
      {
         fx[i] = Fv[i];
         phi += fx[i] * fx[i] / 2;
         zero = zero && fx[i] == 0;
      }
      if (!(BOOST_MATH_NAMESPACE::isfinite)(phi))
         phi = std::numeric_limits<Real>::infinity();
      return true;
   };

   // The tuple holds a closure for F''[v, v] for Halley, which can be copied but not assigned.
   std::unique_ptr<R> r(new R(f(x)));
   Real phi;
   Real last = std::numeric_limits<Real>::infinity();
   if (!unpack(*r, phi))
      return fail(x, function, "F(x) must have x.size() entries and J(x) x.size() rows of x.size()");
   if (phi == std::numeric_limits<Real>::infinity())
      return fail(x, function, "F is not finite at the initial guess");
   while (max_iter < limit && !zero)
   {
      // The gradient g = J^T F of |F|^2/2, and the Cauchy step -(|g|^2/|Jg|^2) g, the least |F + J p| along -g,
      // before the LU overwrites J.
      Real gg = 0, jgjg = 0;
      for (std::size_t j = 0; j < n; ++j)
      {
         g[j] = 0;
         for (std::size_t i = 0; i < n; ++i)
            g[j] += lu[i * n + j] * fx[i];
         gg += g[j] * g[j];
      }
      for (std::size_t i = 0; i < n; ++i)
      {
         Real jg = 0;
         for (std::size_t j = 0; j < n; ++j)
            jg += lu[i * n + j] * g[j];
         jgjg += jg * jg;
      }
      Real const cauchy = gg / jgjg;
      Real const cnorm = cauchy * inf_norm(g);
      if (!lu_factor(lu, piv, colmax, n))
         return fail(x, function, "The Jacobian is singular");
      for (std::size_t i = 0; i < n; ++i)
         a[i] = -fx[i];
      lu_solve(lu, piv, n, a);
      ++max_iter;
      Real const xnorm = inf_norm(x);
      Real const anorm = inf_norm(a);
      if (anorm <= tol * xnorm)
      {
         for (std::size_t i = 0; i < n; ++i)
            x[i] += a[i];
         break;
      }
      // Steps at the rounding level of x that have stopped shrinking are noise: x is as good as it gets. (At a
      // multiple root the steps shrink linearly, so the iteration goes on.)
      if (anorm <= sqrt(eps) * xnorm && anorm >= last)
         break;
      last = anorm;
      Real phit = 0;
      bool newton = true;
      BOOST_MATH_IF_CONSTEXPR (Halley)
      {
         if (!second_derivative(std::integral_constant<bool, Halley>(), *r, a, b))
            return fail(x, function, "d2(v) must have x.size() entries");
         lu_solve(lu, piv, n, b);
         for (std::size_t i = 0; i < n; ++i)
         {
            Real const c = halley_component(a[i], b[i]);
            newton = newton && c == a[i];
            xt[i] = x[i] + c;
         }
         if (!newton)
         {
            std::unique_ptr<R> rt(new R(f(xt)));
            if (!unpack(*rt, phit))
               return fail(x, function, "F(x) must have x.size() entries and J(x) x.size() rows of x.size()");
            if (phit < (1 - 2 * armijo) * phi || phit == 0)
            {
               x = xt;
               r = std::move(rt);
               phi = phit;
               continue;
            }
         }
      }
      // Backtrack along the Newton step until |F|^2/2 decreases enough. On the line x + t a, |F|^2/2 has
      // slope g^T a at t = 0, -2 phi for the Newton step; minimize the quadratic through that and the last trial,
      // kept in [t/10, t/2]. As in a dogleg step, once t a is shorter than the Cauchy step, restart along that.
      Real slope = -2 * phi;
      Real pnorm = anorm;
      bool steep = false;
      Real t = 1;
      for (;;)
      {
         for (std::size_t i = 0; i < n; ++i)
            xt[i] = x[i] + t * a[i];
         std::unique_ptr<R> rt(new R(f(xt)));
         if (!unpack(*rt, phit))
            return fail(x, function, "F(x) must have x.size() entries and J(x) x.size() rows of x.size()");
         // Strict, so that a step too short to change |F| is not accepted; phit = 0 is a decrease unless phi underflows.
         if (phit < phi + armijo * t * slope || phit == 0)
         {
            x = xt;
            r = std::move(rt);
            phi = phit;
            break;
         }
         if (t * pnorm <= tol * xnorm || t < eps)
         {
            // At a local minimum of |F| that is not a root J is singular, so the Newton step is large or undefined.
            // A step below eps^(1/3) |x| that cannot decrease |F| is instead the rounding noise of an accurate x,
            // about cond(J) eps |x|, for cond(J) up to eps^(-2/3), 1e10 in double.
            if (anorm <= noise * xnorm)
               return x;
            return fail(x, function, "The line search cannot decrease |F|: x is near a local minimum of |F| that is not a root, or J is wrong");
         }
         Real s = t / 10;
         if (phit != std::numeric_limits<Real>::infinity())
         {
            Real const curvature = phit - phi - slope * t;
            if (curvature > 0)
               s = -slope * t * t / (2 * curvature);
         }
         t = (std::min)((std::max)(s, t / 10), t / 2);
         if (!steep && t * anorm < cnorm && (BOOST_MATH_NAMESPACE::isfinite)(cnorm))
         {
            steep = true;
            for (std::size_t i = 0; i < n; ++i)
               a[i] = -cauchy * g[i];
            slope = -cauchy * gg;
            pnorm = cnorm;
            t = 1;
         }
      }
   }
   return x;
}

} // namespace detail

// A root of F near x. f(x) returns a pair or tuple (F(x), J(x)), with J a container of rows or a row-major array
// of n*n entries. On exit max_iter is the number of iterations taken.
BOOST_MATH_EXPORT template <class F, class Container>
Container multivariate_newton(F f, Container x, int digits, std::uintmax_t& max_iter)
{
   return detail::multivariate_root_finder<false>(f, std::move(x), digits, max_iter,
      "boost::math::tools::multivariate_newton<%1%>");
}

BOOST_MATH_EXPORT template <class F, class Container>
Container multivariate_newton(F f, Container x, int digits)
{
   std::uintmax_t max_iter = (std::numeric_limits<std::uintmax_t>::max)();
   return multivariate_newton(f, std::move(x), digits, max_iter);
}

// As multivariate_newton, but f(x) returns (F(x), J(x), d2), where d2(v) is F''(x)[v, v], the second derivative
// of t -> F(x + t v) at t = 0.
BOOST_MATH_EXPORT template <class F, class Container>
Container multivariate_halley(F f, Container x, int digits, std::uintmax_t& max_iter)
{
   return detail::multivariate_root_finder<true>(f, std::move(x), digits, max_iter,
      "boost::math::tools::multivariate_halley<%1%>");
}

BOOST_MATH_EXPORT template <class F, class Container>
Container multivariate_halley(F f, Container x, int digits)
{
   std::uintmax_t max_iter = (std::numeric_limits<std::uintmax_t>::max)();
   return multivariate_halley(f, std::move(x), digits, max_iter);
}

} // namespace tools
BOOST_MATH_NAMESPACE_END

#endif // BOOST_MATH_TOOLS_MULTIVARIATE_ROOTS_HPP
