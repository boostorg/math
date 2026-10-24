//  Copyright Jacob Hass 2026.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_QUADRATURE_GAUSS_JACOBI_HPP
#define BOOST_MATH_QUADRATURE_GAUSS_JACOBI_HPP

#ifdef _MSC_VER
#pragma once
#endif

#include <boost/math/tools/config.hpp>

#ifdef BOOST_MATH_NO_CXX17_IF_CONSTEXPR
#error "The header <boost/math/quadrature/gauss_jacobi.hpp> requires C++17 or later."
#endif

#ifndef BOOST_MATH_BUILD_MODULE
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <type_traits>
#include <utility>
#include <vector>
#endif
#include <boost/math/constants/constants.hpp>
#include <boost/math/policies/error_handling.hpp>
#include <boost/math/policies/policy.hpp>
#include <boost/math/special_functions/fpclassify.hpp>
#include <boost/math/tools/big_constant.hpp>
#include <boost/math/tools/precision.hpp>
#include <boost/math/quadrature/detail/quadrature_constant.hpp>
#include <boost/math/special_functions/detail/orthogonal_polynomial.hpp>
#include <boost/math/special_functions/beta.hpp>
#include <boost/math/special_functions/digamma.hpp>
#include <boost/math/tools/compensated_arithmetic.hpp>

BOOST_MATH_NAMESPACE_BEGIN namespace detail {

template <class Real>
class jacobi_family
{
public:
   // The anti-Gauss rule can have nodes outside [-1, 1], but within 3/2: doubling the last coefficient b_n <= 1
   // of the Jacobi matrix, whose norm is at most 1, moves its eigenvalues by at most (sqrt(2) - 1) sqrt(b_n).
   jacobi_family(unsigned n, Real alpha_, Real beta_, bool anti_gauss_ = false) : N(n), alpha(alpha_), beta(beta_), anti_gauss(anti_gauss_) {}

   Real p0() const
   {
      return 1;
   }

   Real lower_bound() const
   {
      return anti_gauss ? Real(-1.5) : Real(-1);
   }

   Real upper_bound() const
   {
      return anti_gauss ? Real(1.5) : Real(1);
   }

   unsigned num_roots() const
   {
      return N;
   }

   // 2^(alpha + beta + 1) B(alpha + 1, beta + 1). The sums alpha + 1, beta + 1 and their sum are rounded, and the gamma
   // functions take that to first order, Gamma(s + e) = Gamma(s) (1 + psi(s) e); else the area, and with it every weight,
   // moves by eps s psi(s), as by 7 eps for alpha = -0.9, beta = 10. The gamma functions are also more accurate than
   // beta(), which is out by 9 eps at (11, 11), where they do not overflow.
   Real area() const
   {
      using std::pow;
      const Real a = alpha + 1;
      const Real b = beta + 1;
      const Real c = a + b;
      const Real ea = BOOST_MATH_NAMESPACE::detail::two_sum_error(alpha, Real(1), a);
      const Real eb = BOOST_MATH_NAMESPACE::detail::two_sum_error(beta, Real(1), b);
      const Real ec = BOOST_MATH_NAMESPACE::detail::two_sum_error(a, b, c) + ea + eb;
      const Real psi_c = BOOST_MATH_NAMESPACE::digamma(c);
      const Real B = c < 100
         ? BOOST_MATH_NAMESPACE::tgamma(a) * BOOST_MATH_NAMESPACE::tgamma(b) / BOOST_MATH_NAMESPACE::tgamma(c) * (1 + BOOST_MATH_NAMESPACE::digamma(a) * ea + BOOST_MATH_NAMESPACE::digamma(b) * eb - psi_c * ec)
         : BOOST_MATH_NAMESPACE::beta(a, b) * (1 + (BOOST_MATH_NAMESPACE::digamma(a) - psi_c) * ea + (BOOST_MATH_NAMESPACE::digamma(b) - psi_c) * eb);
      return 2 * pow(Real(2), alpha) * pow(Real(2), beta) * B;
   }

   unsigned weight_multiplicity(unsigned) const
   {
      return 1;
   }

   Real next(const Real& x, const Real& p, const Real& p_previous, unsigned k) const
   {  // See https://dlmf.nist.gov/18.9
      if (k == 0)
      {
         return ((alpha - beta) + (alpha + beta + 2) * x) / 2 * p;
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

   Real slope(unsigned k) const
   {  // The coefficient of x p in next()
      const Real degree = static_cast<Real>(k);
      const Real sum = alpha + beta;
      if (k == 0)
      {
         return (sum + 2) / 2;
      }
      return (2 * degree + sum + 1) * (2 * degree + sum + 2) / (2 * (degree + 1) * (degree + sum + 1));
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
         2 * (degree + alpha) * (degree + beta) * p_previous) / ((2 * degree + sum) * (1 - x * x));
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

   // For glaser_liu_rokhlin: u = (1 - x)^((alpha + 1) / 2) (1 + x)^((beta + 1) / 2) p_N satisfies s u'' + g u = 0 with
   // s = (1 - x)^2 (1 + x)^2 and g = kappa (1 - x) (1 + x) + A (1 + x)^2 + B (1 - x)^2, where A = (1 - alpha^2) / 4,
   // B = (1 - beta^2) / 4 and kappa = ((2N + alpha + beta + 1)^2 + 1 - alpha^2 - beta^2) / 4. The coefficients of
   // s(t + h) and g(t + h) in h, from d = 1 - t and e = 1 + t, which keep their accuracy near the ends.

   void normal_form(const double_word_arithmetic<Real>& dw, const double_word<Real>& t, double_word<Real> (&s)[5], double_word<Real> (&g)[3], unsigned m) const
   {
      using word = double_word<Real>;
      const word alpha2 = dw.multiply(word{ alpha, Real(0) }, alpha);
      const word beta2 = dw.multiply(word{ beta, Real(0) }, beta);
      const word A = dw.divide(dw.add(dw.negate(alpha2), Real(1)), Real(4));
      const word B = dw.divide(dw.add(dw.negate(beta2), Real(1)), Real(4));
      const word nu = dw.add(dw.two_sum(alpha, beta), Real(2) * Real(m) + 1);
      const word kappa = dw.divide(dw.add(dw.add(dw.multiply(nu, nu), dw.negate(dw.add(alpha2, beta2))), Real(1)), Real(4));
      const word d = dw.add(dw.negate(t), Real(1));
      const word e = dw.add(t, Real(1));
      // (d - h)(e + h) = q0 + q1 h - h^2
      const word q0 = dw.multiply(d, e);
      const word q1 = dw.add(d, dw.negate(e));
      s[0] = dw.multiply(q0, q0);
      s[1] = dw.multiply(q0, word{ 2 * q1.hi, 2 * q1.lo });
      s[2] = dw.add(dw.multiply(q1, q1), word{ -2 * q0.hi, -2 * q0.lo });
      s[3] = word{ -2 * q1.hi, -2 * q1.lo };
      s[4] = word{ Real(1), Real(0) };
      const word Ae = dw.multiply(A, e);
      const word Bd = dw.multiply(B, d);
      g[0] = dw.add(dw.add(dw.multiply(kappa, q0), dw.multiply(Ae, e)), dw.multiply(Bd, d));
      const word difference = dw.add(Ae, dw.negate(Bd));
      g[1] = dw.add(dw.multiply(kappa, q1), word{ 2 * difference.hi, 2 * difference.lo });
      g[2] = dw.add(dw.add(A, B), dw.negate(kappa));
   }

   Real singular_distance(const Real& x) const
   {
      return (std::min)(Real(1 - x), Real(1 + x));
   }

   // (1 - z)^alpha (1 + z)^beta at z = hi + lo.
   Real weight_function(const double_word<Real>& z) const
   {
      using std::pow;
      const double_word<Real> d = double_word_arithmetic<Real>::add(double_word_arithmetic<Real>::negate(z), Real(1));
      const double_word<Real> e = double_word_arithmetic<Real>::add(z, Real(1));
      return pow(d.hi, alpha) * pow(e.hi, beta) * (1 + alpha * d.lo / d.hi + beta * e.lo / e.hi);
   }

   // The series of P_N about either end, in powers of w with c_0 = 1, x = 1 - 2 lambda w or -1 + 2 lambda w: P_N(x) is a
   // multiple of 2F1(-N, N + alpha + beta + 1; alpha + 1; (1 - x) / 2), and of the same with alpha and beta swapped in
   // (1 + x) / 2. The first end is the one where the exponent of the weight is smaller. lambda is a power of 2 near
   // 1 / (2N + alpha + beta + 1)^2, which puts the nearest zero at w of order 1.
   static constexpr unsigned ends = 2;

   Real series_origin(unsigned end) const
   {
      return at_one(end) ? Real(1) : Real(-1);
   }

   Real series_step(unsigned end) const
   {
      using std::frexp;
      using std::ldexp;
      int exponent;
      frexp(2 * Real(N) + alpha + beta + 1, &exponent);
      return ldexp(at_one(end) ? Real(-2) : Real(2), -2 * exponent);
   }

   double_word<Real> series_ratio(const double_word_arithmetic<Real>& dw, unsigned k, unsigned end, unsigned m) const
   {  // c_{k+1} / c_k = lambda (k - m) (k + m + alpha + beta + 1) / ((k + 1) (k + 1 + alpha)) for P_m, alpha and beta swapped at -1
      using std::abs;
      const double_word<Real> numerator = dw.multiply(dw.add(dw.two_sum(alpha, beta), Real(k + m + 1)), (Real(k) - Real(m)) * abs(series_step(end)) / 2);
      const double_word<Real> denominator = dw.multiply(dw.two_sum(Real(k + 1), at_one(end) ? alpha : beta), Real(k + 1));
      return dw.divide(numerator, denominator);
   }

   // P_N = (A x + B) P_{N-1} - C P_{N-2} (DLMF 18.9.2), with k = N - 1, s = alpha + beta and d = 2 (k + 1) (k + s + 1) (2k + s):
   // A = (2k + s + 1)(2k + s)(2k + s + 2) / d, B = (2k + s + 1)(alpha^2 - beta^2) / d, C = 2 (k + alpha)(k + beta)(2k + s + 2) / d.
   void anti_gauss_recurrence(const double_word_arithmetic<Real>& dw, double_word<Real>& A, double_word<Real>& B, double_word<Real>& C2) const
   {
      using word = double_word<Real>;
      const Real k = Real(N - 1);
      const word sum = dw.two_sum(alpha, beta);
      const word m0 = dw.add(sum, 2 * k), m1 = dw.add(sum, 2 * k + 1), m2 = dw.add(sum, 2 * k + 2);
      const word d = dw.multiply(dw.multiply(dw.add(sum, k + 1), m0), 2 * (k + 1));
      A = dw.divide(dw.multiply(dw.multiply(m1, m0), m2), d);
      B = dw.divide(dw.multiply(m1, dw.multiply(dw.two_sum(alpha, -beta), sum)), d);
      const word C = dw.divide(dw.multiply(dw.multiply(dw.two_sum(k, alpha), dw.two_sum(k, beta)), m2), d);
      C2 = word{ -4 * C.hi, -4 * C.lo };
   }

   // P_{N-1}(e) / P_{N-2}(e) = (N - 1 + alpha) / (N - 1) at e = 1, -(N - 1 + beta) / (N - 1) at -1
   double_word<Real> anti_gauss_ratio(const double_word_arithmetic<Real>& dw, unsigned end) const
   {
      const double_word<Real> r = dw.divide(dw.two_sum(Real(N - 1), at_one(end) ? alpha : beta), Real(N - 1));
      return at_one(end) ? r : dw.negate(r);
   }

   // sigma' / sigma = -(alpha + 1) / (2 (1 - x)) + (beta + 1) / (2 (1 + x))
   double_word<Real> scale_log_derivative(const double_word_arithmetic<Real>& dw, const double_word<Real>& x) const
   {
      const double_word<Real> d = dw.add(dw.negate(x), Real(1)), e = dw.add(x, Real(1));
      return dw.add(dw.divide(dw.two_sum(beta, Real(1)), double_word<Real>{ 2 * e.hi, 2 * e.lo }),
                    dw.negate(dw.divide(dw.two_sum(alpha, Real(1)), double_word<Real>{ 2 * d.hi, 2 * d.lo })));
   }

   // The leading coefficient 1 - x^2 of the Sturm-Liouville form ((1 - x)^(alpha + 1) (1 + x)^(beta + 1) p')' = -lambda w p.
   Real sturm_liouville(const double_word<Real>& z) const
   {
      const double_word_arithmetic<Real> dw;
      return dw.multiply(dw.add(dw.negate(z), Real(1)), dw.add(z, Real(1))).hi;
   }

private:
   bool at_one(unsigned end) const
   {
      return (alpha <= beta) == (end == 0);
   }

   unsigned N;
   Real alpha;
   Real beta;
   bool anti_gauss;
};

} // namespace detail
BOOST_MATH_NAMESPACE_END

BOOST_MATH_NAMESPACE_BEGIN namespace quadrature { namespace detail {

template <class Real, unsigned Category>
class gauss_jacobi_detail
{
   unsigned N;
   Real alpha;
   Real beta;
   std::vector<Real> abscissa_vals;
   std::vector<Real> weight_vals;
   bool anti_gauss_internal;
   lazy_rule<Real> anti;
   using polynomial = BOOST_MATH_NAMESPACE::detail::orthogonal_polynomial<Real, BOOST_MATH_NAMESPACE::detail::jacobi_family<Real> >;
public:
   gauss_jacobi_detail(unsigned N_, Real alpha_, Real beta_) : N(N_), alpha(alpha_), beta(beta_)
   {
      // gauss_jacobi raises the domain error.
      if (!in_domain())
      {
         abscissa_vals.assign(N, std::numeric_limits<Real>::quiet_NaN());
         weight_vals.assign(N, std::numeric_limits<Real>::quiet_NaN());
         anti_gauss_internal = false;
         return;
      }
      BOOST_MATH_NAMESPACE::detail::jacobi_family<Real> family(N, alpha, beta);
      polynomial value(N, family);
      abscissa_vals = value.abscissa();
      weight_vals = value.weights();
      // Laurie 1996, Theorem 3: the anti-Gauss rule has a node outside [-1, 1] when alpha or beta is below -1/2
      // and N is large, and for some N when it is near -1/2. At -1/2 it can have nodes on -1 and 1.
      anti_gauss_internal = polynomial::anti_gauss_is_internal(N + 1, family, Real(1)) && polynomial::anti_gauss_is_internal(N + 1, family, Real(-1));
    }
   const std::vector<Real>& abscissa() const
   {
      return abscissa_vals;
   }

   const std::vector<Real>& weights() const
   {
      return weight_vals;
   }

   bool anti_gauss_is_internal() const
   {
      return anti_gauss_internal;
   }

   // Laurie's anti-Gauss rule with N + 1 points, for integrate_averaged; empty unless it is internal.
   const std::vector<Real>& anti_gauss_abscissa() const
   {
      return anti_gauss().x;
   }

   const std::vector<Real>& anti_gauss_weights() const
   {
      return anti_gauss().w;
   }

protected:
   const auto& anti_gauss() const
   {
      return anti.get([this](std::vector<Real>& x, std::vector<Real>& w)
      {
         if (!anti_gauss_internal)
            return;
         const polynomial rule(N + 1, BOOST_MATH_NAMESPACE::detail::jacobi_family<Real>(N + 1, alpha, beta, true), true);
         x = rule.abscissa();
         w = rule.weights();
         // Nodes on the endpoints come out a few eps either side.
         for (Real& z : x)
            z = (std::min)((std::max)(z, Real(-1)), Real(1));
      });
   }

   bool in_domain() const
   {
      return N > 0 && alpha > -1 && beta > -1 && (BOOST_MATH_NAMESPACE::isfinite)(alpha) && (BOOST_MATH_NAMESPACE::isfinite)(beta);
   }
};

} // namespace detail

template <class Real, class Policy = BOOST_MATH_NAMESPACE::policies::policy<> >
class gauss_jacobi : public detail::gauss_jacobi_detail<Real, detail::quadrature_constant_category<Real>::value>
{
   template <class F, class K, class Norm>
   K anti_gauss_average(F& f, K zero, const K& result, Real& L1, Norm norm, Real* error) const
   {
      if (this->in_domain() && !this->anti_gauss_is_internal())
      {
         L1 = policies::raise_evaluation_error("boost::math::quadrature::gauss_jacobi<%1%>::integrate_averaged",
            "The anti-Gauss rule has a node outside [-1, 1] for these alpha and beta, so there is no averaged rule.",
            std::numeric_limits<Real>::quiet_NaN(), Policy());
         if (error)
            *error = L1;
         return zero * L1;
      }
      return detail::anti_gauss_average(f, this->anti_gauss_abscissa(), this->anti_gauss_weights(), false, zero, result, L1, norm, error, "boost::math::quadrature::gauss_jacobi<%1%>::integrate_averaged", Policy());
   }

public:
   gauss_jacobi(unsigned N_, Real alpha_, Real beta_) : detail::gauss_jacobi_detail<Real, detail::quadrature_constant_category<Real>::value>(N_, alpha_, beta_)
   {
      if (N_ == 0)
      {
         policies::raise_domain_error("boost::math::quadrature::gauss_jacobi<%1%>::gauss_jacobi",
            "The number of points must be > 0, but got %1%.", Real(0), Policy());
      }
      if (!(alpha_ > -1) || !(BOOST_MATH_NAMESPACE::isfinite)(alpha_))
      {
         policies::raise_domain_error("boost::math::quadrature::gauss_jacobi<%1%>::gauss_jacobi",
            "The parameter alpha must be finite and > -1, but got %1%.", alpha_, Policy());
      }
      if (!(beta_ > -1) || !(BOOST_MATH_NAMESPACE::isfinite)(beta_))
      {
         policies::raise_domain_error("boost::math::quadrature::gauss_jacobi<%1%>::gauss_jacobi",
            "The parameter beta must be finite and > -1, but got %1%.", beta_, Policy());
      }
   }

   // Integrates f(x) (1 - x)^alpha * (1 + x)^beta over [-1, 1].
   template <class F>
   auto integrate(F f, Real* pL1 = nullptr) const->decltype(std::declval<F>()(std::declval<Real>()))
   {
      using K = decltype(f(Real(0)));
      return integrate(f, K(Real(0)), [](const K& z) { using std::abs; return abs(z); }, pL1);
   }

   // Explicit zero and norm for vector- and matrix-valued integrands.
   template <class F, class Norm>
   auto integrate(F f, const decltype(std::declval<F>()(std::declval<Real>()))& zero, Norm norm, Real* pL1 = nullptr) const
      ->decltype(static_cast<Real>(norm(f(Real(0)))), f(Real(0)))
   {
      using K = decltype(f(Real(0)));
      static_assert(!std::is_integral<K>::value,
                    "The return type cannot be integral, it must be either a real or complex floating point type.");
      Real nan = std::numeric_limits<Real>::quiet_NaN();
      if (!this->in_domain() || detail::overflowed(this->abscissa(), "boost::math::quadrature::gauss_jacobi<%1%>::integrate", pL1, &nan, Policy()))
      {
         if (pL1)
            *pL1 = nan;
         K result = zero * nan;
         return result;
      }
      K result = zero;
      Real L1 = static_cast<Real>(norm(result));
      for (unsigned i = 0; i < this->abscissa().size(); ++i)
      {
         K fp = f(static_cast<Real>(this->abscissa()[i]));
         result += fp * static_cast<Real>(this->weights()[i]);
         L1 += static_cast<Real>(norm(fp)) * static_cast<Real>(this->weights()[i]);
      }
      if (pL1)
         *pL1 = L1;
      return result;
   }
   // Laurie's averaged rule (G + A) / 2: N + 1 more evaluations of f than integrate, exact for polynomials of degree
   // 2N + 1. The error is half the difference of G and A, a bound on the error when they bracket the integral.
   template <class F>
   auto integrate_averaged(F f, Real* error = nullptr, Real* pL1 = nullptr) const->decltype(std::declval<F>()(std::declval<Real>()))
   {
      using K = decltype(f(Real(0)));
      return integrate_averaged(f, K(Real(0)), [](const K& z) { using std::abs; return abs(z); }, error, pL1);
   }

   template <class F, class Norm>
   auto integrate_averaged(F f, const decltype(std::declval<F>()(std::declval<Real>()))& zero, Norm norm, Real* error = nullptr, Real* pL1 = nullptr) const
      ->decltype(static_cast<Real>(norm(f(Real(0)))), f(Real(0)))
   {
      Real L1;
      auto result = integrate(f, zero, norm, &L1);
      result = anti_gauss_average(f, zero, result, L1, norm, error);
      if (pL1)
         *pL1 = L1;
      return result;
   }
};
} // namespace quadrature

BOOST_MATH_NAMESPACE_END

#endif // BOOST_MATH_QUADRATURE_GAUSS_JACOBI_HPP
