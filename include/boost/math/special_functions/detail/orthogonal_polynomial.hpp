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
#include <boost/math/special_functions/trunc.hpp>
#include <boost/math/tools/compensated_arithmetic.hpp>
#include <boost/math/policies/error_handling.hpp>
#include <boost/math/policies/policy.hpp>

BOOST_MATH_NAMESPACE_BEGIN namespace detail {

// The zeros of p_N and their Gauss weights in O(N) operations, after Glaser, Liu & Rokhlin, SIAM J. Sci. Comput. 29
// (2007), or those of Laurie's anti-Gauss polynomial p~_N = (A x + B) p_{N-1} - 2 C p_{N-2} (where p_N = (A x + B) p_{N-1}
// - C p_{N-2}) and its weights. The family supplies the normal form s u'' + g u = 0 of the equation of u = sigma p_m, with
// sigma^2 the weight function W times the leading coefficient A_SL of the Sturm-Liouville form: u is bounded where p_m
// oscillates, unlike p_m. The zeros are those of F = (B + A x) u_0 + C_2 u_1: F = u_0 = sigma p_N for the Gauss rule, and
// u_0 = sigma p_{N-1}, u_1 = sigma p_{N-2}, C_2 = -2 C for the anti-Gauss rule.
//
// From a zero, a Pruefer transform of the equation of u_0 predicts the next zero of p_N, or brackets that of p~_N between
// the next two of p_{N-1}, which interlace with them (Laurie 1996). Newton's method on the Taylor series of F, whose
// coefficients follow from the equations, finds it, and the series carry u and u' there, where the Gauss weight is
// W / F'^2 and the anti-Gauss weight W A_SL / (u_0 F'), up to a constant. The series are summed in double words, as the
// rounding errors of the steps add up, and through W move the weights by thousands of eps. They converge out to the ends
// of a finite interval, and near those the zeros come instead from the series of the polynomials about the end.
//
// solve() returns false where it cannot vouch for the result, and the caller falls back on bisection.
template <class Real, class Family>
class glaser_liu_rokhlin
{
   using arithmetic = double_word_arithmetic<Real>;
   using word = double_word<Real>;
   // u = sigma p_m at the anchor t, and its scaled Taylor coefficients a_k = H^k u^(k)(t) / k!, b_k = k (k - 1) a_k there.
   struct solution
   {
      unsigned degree;
      word u, u_prime, g[3];
      std::vector<word> a, b;
      std::size_t terms;
   };

public:
   glaser_liu_rokhlin(const Family& family_, unsigned n, bool anti_gauss) : family(family_), N(n), anti(anti_gauss)
   {
      using std::ldexp;
      tolerance = ldexp(tools::epsilon<Real>() * tools::epsilon<Real>(), -2);
      max_terms = 8 * static_cast<unsigned>(tools::digits<Real>()) + 64;
      u[0].degree = anti ? N - 1 : N;
      u[1].degree = N - 2;
      u[1].u = u[1].u_prime = word{ Real(0), Real(0) };
      if (anti)
         family.anti_gauss_recurrence(dw, A, B, C2);
   }

   // The zeros in increasing order, the nonnegative ones only for a symmetric family, and their weights up to a constant factor.
   bool solve(std::vector<Real>& x, std::vector<Real>& w)
   {
      BOOST_MATH_STD_USING
      const unsigned num_roots = family.num_roots();
      x.clear();
      w.clear();
      word t;
      int direction = 1;
      bool from_zero = true;
      // The zeros found from the far end, from the end inwards, and their weights up to another constant factor.
      std::vector<word> far_zeros;
      std::vector<Real> far_weights;
      if (!begin(std::integral_constant<bool, Family::ends == 0>(), t, direction, from_zero, x, w, far_zeros, far_weights))
         return false;
      // The innermost of the far zeros is found twice, and its two weights give the ratio of the constants.
      unsigned count = num_roots - static_cast<unsigned>(x.size() + far_zeros.size()) + (far_zeros.empty() ? 0 : 1);
      for (; count > 0; --count)
      {
         if (!step(t, direction, from_zero))
            return false;
         x.push_back(t.hi);
         w.push_back(weight(t));
      }
      if (!far_zeros.empty())
      {
         if (!(abs(x.back() - far_zeros.back().hi) <= 4 * tools::epsilon<Real>() * (std::max)(Real(abs(x.back())), Real(1))))
            return false;
         const Real scale = w.back() / far_weights.back();
         for (std::size_t i = far_zeros.size() - 1; i-- > 0;)
         {
            x.push_back(far_zeros[i].hi);
            w.push_back(far_weights[i] * scale);
         }
      }
      if (direction < 0)
      {
         std::reverse(x.begin(), x.end());
         std::reverse(w.begin(), w.end());
      }
      for (std::size_t i = 0; i < x.size(); ++i)
      {
         if (!(BOOST_MATH_NAMESPACE::isfinite)(x[i]) || !(BOOST_MATH_NAMESPACE::isfinite)(w[i]) || !(w[i] >= 0)
            || ((i > 0) && !(x[i] > x[i - 1])) || (x[i] < family.lower_bound()) || (x[i] > family.upper_bound()))
            return false;
      }
      return true;
   }

private:
   const Family& family;
   unsigned N;
   bool anti;
   word A{ Real(0), Real(0) }, B{ Real(1), Real(0) }, C2{ Real(0), Real(0) };
   arithmetic dw;
   Real tolerance;
   unsigned max_terms;
   solution u[2];
   std::vector<word> reciprocals;
   // s = the product of (x - e)^2 over the ends e of the interval
   static constexpr unsigned s_degree = 2 * Family::ends;
   word s[5];
   Real last_step = 0;

   unsigned solutions() const
   {
      return anti ? 2 : 1;
   }

   // B + A x
   word factor(const word& x) const
   {
      return arithmetic::add(dw.multiply(x, A), B);
   }

   // The weight at a zero z of F, where u and u' are known.
   Real weight(const word& z) const
   {
      BOOST_MATH_STD_USING
      const Real F_prime = A.hi * u[0].u.hi + factor(z).hi * u[0].u_prime.hi + C2.hi * u[1].u_prime.hi;
      // The anti-Gauss weights share the sign of a constant factor.
      return anti ? family.weight_function(z) * family.sturm_liouville(z) / abs(u[0].u.hi * F_prime)
                  : family.weight_function(z) / (F_prime * F_prime);
   }

   // A symmetric family: p_m is even or odd, and so is u = sigma p_m. Start at 0, where F = 0 when it is odd.
   bool begin(std::true_type, word& t, int& direction, bool& from_zero, std::vector<Real>& x, std::vector<Real>& w,
              std::vector<word>&, std::vector<Real>&)
   {
      t = word{ Real(0), Real(0) };
      direction = 1;
      // p_{N-1} and p_{N-2} with the ratio of their first nonzero derivatives at 0.
      const word rho = anti ? family.anti_gauss_ratio(dw, 0) : word{ Real(1), Real(0) };
      for (unsigned i = 0; i < solutions(); ++i)
      {
         const word one = i == 0 ? rho : word{ Real(1), Real(0) };
         const bool odd = u[i].degree % 2;
         u[i].u = odd ? word{ Real(0), Real(0) } : one;
         u[i].u_prime = odd ? one : word{ Real(0), Real(0) };
      }
      from_zero = N % 2;
      if (from_zero)
      {
         x.push_back(Real(0));
         w.push_back(weight(t));
      }
      return true;
   }

   // Otherwise the zeros nearest each end come from the series about it, and the march starts from the last of those at
   // the first end, with sigma = 1 there. Then the weights of the zeros from the series, a constant over A_SL p_N'^2 or
   // p_{N-1} p~_N', are W A_SL at the last zero times those over A_SL u_0'^2 or u_0 F'.
   bool begin(std::false_type, word& t, int& direction, bool& from_zero, std::vector<Real>& x, std::vector<Real>& w,
              std::vector<word>& far_zeros, std::vector<Real>& far_weights)
   {
      std::vector<word> zeros;
      std::vector<Real> weights;
      if (!series_zeros(0, N, zeros, weights, true))
         return false;
      t = zeros.back();
      const Real scale = zeros.size() < N ? family.weight_function(t) * family.sturm_liouville(t) : Real(1);
      for (std::size_t i = 0; i < zeros.size(); ++i)
      {
         x.push_back(zeros[i].hi);
         w.push_back(scale * weights[i]);
      }
      direction = family.series_step(0) > 0 ? 1 : -1;
      from_zero = true;
      // At least one zero from the far end, and at most up to the last from the first; where they overlap, the march has
      // nothing to add, and the last zero from the first end is the one found twice.
      if ((Family::ends == 2) && (zeros.size() < N))
         series_zeros(1, N - static_cast<unsigned>(zeros.size()) + 1, far_zeros, far_weights, false);
      return true;
   }

   // The zeros of F nearest the given end, from the series of the p_m about it, in double words, and their weights up to a
   // constant factor, at least one and at most max_count. Laguerre's method from the end, with the zeros already found
   // divided out, climbs monotonically to the next, as F has N real zeros, all on one side; Newton's method finishes it.
   // The series cancel more the further out the zero: the march takes over where they would cost more than a few eps,
   // or where its steps need not be shortened, the zeros spaced by less than a third of their distance from the end (a
   // quarter for the anti-Gauss rule, whose steps span about one and a half spacings).
   // With anchor, u and u' at the last zero are set for the march.
   bool series_zeros(unsigned end, unsigned max_count, std::vector<word>& zeros, std::vector<Real>& weights, bool anchor)
   {
      BOOST_MATH_STD_USING
      const Real eps = tools::epsilon<Real>();
      const Real c_step = family.series_step(end);
      const Real origin = family.series_origin(end);
      // F = (B + A x) rho y_0 + C_2 y_1 with y_m = p_m / p_m(end), x = origin + c_step v.
      const word rho = anti ? family.anti_gauss_ratio(dw, end) : word{ Real(1), Real(0) };
      // The ratios r_k = c_{k+1} / c_k of the coefficients of each y in v: the terms t_k = c_k v^k are summed from them,
      // as the coefficients themselves fall out of range long before the terms do at the outer zeros.
      std::vector<word> r[2];
      auto ratio = [&](unsigned i, std::size_t k) -> word
      {
         while (r[i].size() <= k)
            r[i].push_back(family.series_ratio(dw, static_cast<unsigned>(r[i].size()), end, u[i].degree));
         return r[i][k];
      };
      // y, y' and y'' in v at v, until two terms in a row fall below the double-word precision, then F, dF/dv and d2F/dv2,
      // and the sum of the magnitudes of the terms of F; in double words, or for words = false in working precision, which
      // is enough for Laguerre's method except where F, as for the anti-Gauss rule, is a small difference of its terms.
      word y[2][3], F[3];
      Real magnitude;
      auto evaluate_series = [&](const word& v, bool words) -> bool
      {
         magnitude = 0;
         for (unsigned i = 0; i < solutions(); ++i)
         {
            const word scale = i == 0 ? rho : C2;
            word term{ Real(1), Real(0) }, sum[3] = {};
            Real big = 0;
            for (std::size_t k = 0, small = 0; ; ++k)
            {
               for (unsigned d = 0; d < 3; ++d)
               {
                  const Real c = d == 0 ? Real(1) : d == 1 ? Real(k) : Real(k) * (Real(k) - 1);
                  // y'' only steers Laguerre's method, which needs it in double words only for the anti-Gauss rule.
                  sum[d] = (words && ((d < 2) || anti)) ? arithmetic::add(sum[d], dw.multiply(term, c)) : word{ sum[d].hi + c * term.hi, Real(0) };
               }
               magnitude += abs(term.hi * scale.hi);
               big = (std::max)(big, abs(term.hi));
               small = ((k > 1) && (abs(term.hi) <= tolerance * big)) ? small + 1 : 0;
               if ((small == 2) || (k > u[i].degree))
                  break;
               if (k >= max_terms)
                  return false;
               term = words ? dw.multiply(dw.multiply(term, ratio(i, k)), v) : word{ term.hi * ratio(i, k).hi * v.hi, Real(0) };
            }
            // y' and y'' in v from v y' and v^2 y'', or at v = 0 from the coefficients; y_0 times rho.
            y[i][0] = sum[0];
            y[i][1] = v.hi == 0 ? ratio(i, 0) : dw.divide(sum[1], v);
            y[i][2] = v.hi == 0 ? dw.multiply(dw.multiply(ratio(i, 0), ratio(i, 1)), Real(2)) : dw.divide(dw.divide(sum[2], v), v);
            for (unsigned d = 0; (i == 0) && (d < 3); ++d)
               y[0][d] = dw.multiply(y[0][d], rho);
         }
         const word b = factor(arithmetic::add(dw.multiply(v, c_step), origin));
         const word a{ A.hi * c_step, A.lo * c_step };
         F[0] = dw.multiply(y[0][0], b);
         F[1] = arithmetic::add(dw.multiply(y[0][0], a), dw.multiply(y[0][1], b));
         F[2] = arithmetic::add(dw.multiply(y[0][1], word{ 2 * a.hi, 2 * a.lo }), dw.multiply(y[0][2], b));
         for (unsigned d = 0; anti && (d < 3); ++d)
            F[d] = arithmetic::add(F[d], dw.multiply(y[1][d], C2));
         return true;
      };
      zeros.clear();
      weights.clear();
      std::vector<Real> found;
      for (unsigned count = 0; count < max_count; ++count)
      {
         // From the end, then from just past the last zero, clear of the next and of the zeros divided out (v is scaled to
         // put the first zeros of p_N about 1 apart).
         const Real gap = count == 0 ? Real(0) : found.back() - (count > 1 ? found[count - 2] : Real(0));
         word v{ count == 0 ? Real(0) : found.back() + (std::max)(gap, Real(1)) / 64, Real(0) };
         const Real n = static_cast<Real>(N - count);
         // To about sqrt(eps), and one step more, which would take it to eps but for rounding errors where the terms
         // cancel, as for the anti-Gauss rule; Newton's method finishes.
         for (int iteration = 0, close = 0; (iteration < 100) && (close < (anti ? 1 : 2)); ++iteration)
         {
            // Just past a zero, where v starts, the terms of F cancel, and the steps need it to a few digits whatever the
            // precision: in float they would jump past the next zero. Use double words where those digits are lost.
            if (!evaluate_series(v, anti) || (!anti && (iteration == 0) && !(magnitude * eps <= abs(F[0].hi) / 65536) && !evaluate_series(v, true)))
               return !zeros.empty();
            if (F[0].hi == 0)
               break;
            Real G = F[1].hi / F[0].hi;
            Real G2 = G * G - F[2].hi / F[0].hi;
            for (const Real& z : found)
            {
               G -= 1 / (v.hi - z);
               G2 -= 1 / ((v.hi - z) * (v.hi - z));
            }
            const Real d = (std::max)(Real((n - 1) * (n * G2 - G * G)), Real(0));
            const Real next = v.hi - n / (G < 0 ? G - sqrt(d) : G + sqrt(d));
            close += (close > 0) || (abs(next - v.hi) <= sqrt(eps) * abs(next));
            v = word{ next, Real(0) };
         }
         // A zero on the end itself, as anti-Gauss rules have for an exponent -1/2 there.
         const bool on_end = (count == 0) && (abs(v.hi) <= 16 * eps);
         if (on_end)
            v = word{ Real(0), Real(0) };
         else if (!(v.hi > (found.empty() ? Real(0) : found.back())))
            return !zeros.empty();
         // Newton's method in double words, until a step of a few eps, after which y, y' and F' follow to first order.
         bool converged = on_end && evaluate_series(v, true);
         for (int iteration = 0; (iteration < 6) && !converged; ++iteration)
         {
            // The terms of the series cancel; their rounding errors must not reach the working precision.
            if (!evaluate_series(v, true) || !(magnitude * eps <= abs(F[1].hi * v.hi) / 64))
               return !zeros.empty();
            const Real delta = -(F[0].hi + F[0].lo) / F[1].hi;
            v = arithmetic::add(v, delta);
            converged = abs(delta) <= 16 * eps * v.hi;
            for (unsigned i = 0; converged && (i < solutions()); ++i)
            {
               y[i][0] = arithmetic::add(y[i][0], dw.multiply(y[i][1], delta));
               y[i][1] = arithmetic::add(y[i][1], dw.multiply(y[i][2], delta));
            }
            if (converged)
               F[1] = arithmetic::add(F[1], dw.multiply(F[2], delta));
         }
         if (!converged)
            return !zeros.empty();
         const word z = arithmetic::add(word{ c_step * v.hi, c_step * v.lo }, origin);
         // F' and y' in x.
         const Real F_prime = F[1].hi / c_step;
         zeros.push_back(z);
         weights.push_back(1 / abs(F_prime * (anti ? y[0][0].hi : family.sturm_liouville(z) * F_prime)));
         found.push_back(v.hi);
         if (anchor)
         {
            // u = sigma y with sigma = 1 at z: u' = y' + (sigma' / sigma) y.
            const word r_sigma = family.scale_log_derivative(dw, z);
            for (unsigned i = 0; i < solutions(); ++i)
            {
               u[i].u = y[i][0];
               u[i].u_prime = arithmetic::add(dw.divide(y[i][1], c_step), dw.multiply(y[i][0], r_sigma));
            }
            last_step = abs(c_step * (v.hi - (count > 0 ? found[count - 1] : Real(0))));
         }
         if ((count > 0) && ((anti ? 4 : 3) * (v.hi - found[count - 1]) <= v.hi))
            break;
      }
      return true;
   }

   // f = g / s for u_0 and f' at t + h, in working precision.
   void f(const Real& h, Real& value, Real& derivative) const
   {
      Real sv = 0, sd = 0, gv = 0, gd = 0;
      for (unsigned i = s_degree + 1; i-- > 0;)
      {
         sd = sd * h + sv;
         sv = sv * h + s[i].hi;
      }
      for (unsigned i = 3; i-- > 0;)
      {
         gd = gd * h + gv;
         gv = gv * h + u[0].g[i].hi;
      }
      value = gv / sv;
      derivative = (gd - value * sd) / sv;
   }

   // The step h from the anchor to where the Pruefer angle of u_0, tan(theta) = sqrt(f) u / u' (f = g / s), reaches
   // theta1 from theta0: dx/dtheta = 1 / (sqrt(f) + f' sin(2 theta) / (4 f)) is nearly constant where f varies slowly on the
   // scale of the oscillation, and one step of Runge-Kutta per half turn predicts the zero to a fraction of a percent.
   // NaN past a turning point, where there is no next zero.
   Real predict(const Real& theta0, const Real& theta1) const
   {
      BOOST_MATH_STD_USING
      const int steps = (std::max)(BOOST_MATH_NAMESPACE::itrunc(ceil(abs(theta1 - theta0) / constants::pi<Real>() - Real(0.01))), 1);
      const Real d = (theta1 - theta0) / steps;
      auto slope = [&](const Real& h, const Real& theta) -> Real
      {
         Real value, f_prime;
         f(h, value, f_prime);
         const Real denominator = sqrt(value) + f_prime / (4 * value) * sin(2 * theta);
         return denominator > 0 ? 1 / denominator : std::numeric_limits<Real>::quiet_NaN();
      };
      Real h = 0, theta = theta0;
      for (int i = 0; i < steps; ++i)
      {
         const Real k1 = slope(h, theta);
         const Real k2 = slope(h + d / 2 * k1, theta + d / 2);
         const Real k3 = slope(h + d / 2 * k2, theta + d / 2);
         const Real k4 = slope(h + d * k3, theta + d);
         h += d / 6 * (k1 + 2 * k2 + 2 * k3 + k4);
         theta += d;
      }
      return h;
   }

   // The scaled Taylor coefficients of a solution, until the terms fall below the double-word precision for
   // |tau| <= tau_max: from s u'' + g u = 0,
   // b_{m+2} = -(sum_{i>=1} s_i H^i b_{m+2-i} + sum_i g_i H^{i+2} a_{m-i}) / s_0, a_{m+2} = b_{m+2} / ((m + 2) (m + 1)).
   // Returns false if they do not fall.
   bool taylor(solution& v, const Real& H, const Real& tau_max)
   {
      BOOST_MATH_STD_USING
      word st[5], gt[3];
      const word inverse = dw.divide(word{ Real(1), Real(0) }, s[0]);
      Real power = H;
      for (unsigned i = 1; i <= s_degree; ++i, power *= H)
         st[i] = dw.multiply(word{ s[i].hi * power, s[i].lo * power }, inverse);
      power = H * H;
      for (unsigned i = 0; i < 3; ++i, power *= H)
         gt[i] = dw.multiply(word{ v.g[i].hi * power, v.g[i].lo * power }, inverse);
      std::vector<word>& a = v.a;
      std::vector<word>& b = v.b;
      a.resize((std::max)(a.size(), std::size_t(64)));
      b.resize(a.size());
      a[0] = v.u;
      a[1] = word{ v.u_prime.hi * H, v.u_prime.lo * H };
      b[0] = b[1] = word{ Real(0), Real(0) };
      Real big = (std::max)(abs(a[0].hi), abs(a[1].hi) * tau_max);
      power = tau_max;
      unsigned small = 0;
      for (unsigned m = 0; m + 2 < max_terms; ++m)
      {
         if (m + 2 >= a.size())
         {
            a.resize(2 * (m + 2));
            b.resize(2 * (m + 2));
         }
         if (m >= reciprocals.size())
            reciprocals.push_back(dw.divide(word{ Real(1), Real(0) }, static_cast<Real>((m + 2) * (m + 1))));
         // From the oldest coefficients to the newest, which is then the only one the sum waits for.
         Real sum = 0, error = 0;
         for (unsigned i = (std::min)(2u, m) + 1; i-- > 0;)
            dw.multiply_add(sum, error, gt[i], a[m - i]);
         for (unsigned i = (std::min)(s_degree, m + 2); i >= 1; --i)
            dw.multiply_add(sum, error, st[i], b[m + 2 - i]);
         // Left unnormalized, which the products above take in their stride, to keep the renormalizations off the
         // dependency chain from one coefficient to the next.
         b[m + 2] = word{ -sum, -error };
         Real product = 0, product_error = 0;
         dw.multiply_add(product, product_error, b[m + 2], reciprocals[m]);
         a[m + 2] = word{ product, product_error };
         power *= tau_max;
         const Real term = abs(a[m + 2].hi) * power;
         big = (std::max)(big, term);
         // Coefficients can vanish, as the odd ones do for an even u about 0: two in a row must be small.
         small = (term <= tolerance * big) ? small + 1 : 0;
         if (small == 2)
         {
            v.terms = m + 3;
            return true;
         }
      }
      return false;
   }

   // A series at tau and its first derivative in double words, by compensated Horner, its second in working precision.
   void evaluate(const solution& v, const Real& tau, word& value, word& slope, Real& curvature) const
   {
      Real y = 0, d = 0, y_error = 0, d_error = 0;
      curvature = 0;
      for (std::size_t k = v.terms; k > 0; --k)
      {
         curvature = curvature * tau + 2 * d;
         dw.horner_step(d, d_error, tau, word{ y, y_error });
         dw.horner_step(y, y_error, tau, v.a[k - 1]);
      }
      value = arithmetic::two_sum(y, y_error);
      slope = arithmetic::two_sum(d, d_error);
   }

   // F and F' = dF/dtau at t + H tau in working precision.
   void combination(const word& t, const Real& H, const Real& tau, Real& F, Real& F_prime) const
   {
      F = F_prime = 0;
      for (unsigned i = 0; i < solutions(); ++i)
      {
         Real y = 0, d = 0;
         for (std::size_t k = u[i].terms; k-- > 0;)
         {
            d = d * tau + y;
            y = y * tau + u[i].a[k].hi;
         }
         const Real c = i == 0 ? B.hi + A.hi * (t.hi + H * tau) : C2.hi;
         F += c * y;
         F_prime += c * d + (i == 0 ? A.hi * H * y : Real(0));
      }
   }

   // From the anchor t, a zero of F (or for a symmetric family, possibly 0), to the next zero in the direction.
   bool step(word& t, int direction, bool& from_zero)
   {
      BOOST_MATH_STD_USING
      const Real eps = tools::epsilon<Real>();
      const Real pi = constants::pi<Real>();
      for (unsigned i = 0; i < solutions(); ++i)
         family.normal_form(dw, t, s, u[i].g, u[i].degree);
      // The bracket [lo, hi] in h for the next zero, and a first guess.
      Real f0, f0_prime;
      f(Real(0), f0, f0_prime);
      // At a zero of the Gauss F = u_0, theta = 0.
      const Real theta0 = (from_zero && !anti) ? Real(0) : (u[0].u_prime.hi == 0 ? pi / 2 : atan(sqrt(f0) * u[0].u.hi / u[0].u_prime.hi));
      Real lo, hi, guess;
      if (!anti)
      {
         // The next multiple of pi, past t if t is a zero.
         guess = predict(theta0, (from_zero || (theta0 * direction > 0)) ? direction * pi : Real(0));
         lo = 3 * guess / 4;
         hi = 5 * guess / 4;
      }
      else
      {
         // Between the next two zeros of p_{N-1}, the first of them t itself where it is one; past the last of them, as
         // far again as the last step.
         const Real theta_lo = ((u[0].u.hi == 0) || (theta0 * direction <= 0)) ? Real(0) : direction * pi;
         lo = u[0].u.hi == 0 ? Real(0) : predict(theta0, theta_lo);
         hi = predict(theta0, theta_lo + direction * pi);
         guess = predict(theta0, theta_lo + direction * pi / 2);
         // Past the last zero of p_{N-1}, and near the turning point before it, the prediction fails; the last zero of
         // p~_N is then within about a step.
         if (!(BOOST_MATH_NAMESPACE::isfinite)(lo * hi * guess))
         {
            const Real step = last_step > 0 ? last_step : pi / sqrt(f0);
            if (!(BOOST_MATH_NAMESPACE::isfinite)(lo))
               lo = direction * step / 8;
            hi = lo + 3 * direction * step / 2;
            guess = (lo + hi) / 2;
         }
         // Allow for the errors of the prediction.
         const Real margin = abs(hi - lo) / 64;
         lo = direction * (std::max)(Real(lo * direction - margin), Real(0));
         hi += direction * margin;
      }
      if (!(BOOST_MATH_NAMESPACE::isfinite)(lo * hi * guess) || !(lo * direction >= 0) || !(hi * direction > lo * direction))
         return false;
      // Within 5/8 of the distance to the nearest singular point the series converge at least like (5/8)^k; the zeros
      // nearer the ends came from the series about them, so the steps here stay inside.
      if (!(8 * abs(hi) <= 5 * family.singular_distance(t.hi)))
         return false;
      int exponent;
      frexp(hi, &exponent);
      const Real H = direction * ldexp(Real(1), exponent);
      for (unsigned i = 0; i < solutions(); ++i)
         if (!taylor(u[i], H, hi / H))
            return false;
      // Newton's method in working precision; for the anti-Gauss rule, whose guess is rougher, falling back on bisection
      // within the bracket, where F changes sign.
      Real tau_lo = lo / H, tau_hi = hi / H, tau = guess / H, F, F_prime, F_lo = -1, F_hi = 1;
      if (anti)
      {
         combination(t, H, tau_lo, F_lo, F_prime);
         combination(t, H, tau_hi, F_hi, F_prime);
         if (!((F_lo < 0) != (F_hi < 0)))
            return false;
      }
      const bool lo_negative = F_lo < 0;
      bool converged = false;
      for (int i = 0; (i < 100) && !converged; ++i)
      {
         combination(t, H, tau, F, F_prime);
         ((F < 0) == lo_negative ? tau_lo : tau_hi) = tau;
         // Close enough for the step in double words below to square the error; F has a floor of rounding errors above
         // eps where its terms cancel, as for the anti-Gauss rule.
         const Real correction = F / F_prime;
         converged = abs(correction) <= pow(eps, Real(0.75)) * abs(tau);
         if (!converged)
         {
            tau -= correction;
            if (anti && !((tau - tau_lo) * (tau - tau_hi) < 0))
               tau = (tau_lo + tau_hi) / 2;
         }
      }
      // The zero must be the one predicted: the next, not this one or the one after.
      if (!converged || !((tau - lo / H) * (tau - hi / H) < 0))
         return false;
      // F and its slope in double words at tau give the zero to about twice the working precision, and u, u' there to
      // second and first order in the correction delta, of order eps tau.
      word value[2], slope[2];
      Real curvature[2];
      for (unsigned i = 0; i < solutions(); ++i)
         evaluate(u[i], tau, value[i], slope[i], curvature[i]);
      const word c = factor(arithmetic::add(t, word{ H * tau, Real(0) }));
      word F_word = dw.multiply(value[0], c);
      Real F_slope = A.hi * H * value[0].hi + c.hi * slope[0].hi;
      if (anti)
      {
         F_word = arithmetic::add(F_word, dw.multiply(value[1], C2));
         F_slope += C2.hi * slope[1].hi;
      }
      const Real delta = -(F_word.hi + F_word.lo) / F_slope;
      const word tz = arithmetic::fast_two_sum(tau, delta);
      for (unsigned i = 0; i < solutions(); ++i)
      {
         u[i].u = arithmetic::add(arithmetic::add(value[i], dw.multiply(slope[i], delta)), curvature[i] * delta * delta / 2);
         const word d = arithmetic::add(slope[i], curvature[i] * delta);
         u[i].u_prime = word{ d.hi / H, d.lo / H };
      }
      last_step = abs(tz.hi * H);
      t = arithmetic::add(t, word{ tz.hi * H, tz.lo * H });
      from_zero = true;
      return true;
   }
};

template <class Real, class Family, class Policy = BOOST_MATH_NAMESPACE::policies::policy<> >
class orthogonal_polynomial
{
public:
   orthogonal_polynomial(unsigned n) : N(n)
   {
      values = calculate_values(Family(N));
   }

   // With anti_gauss, Laurie's anti-Gauss rule with n nodes, the companion of the (n - 1)-point Gauss rule.
   orthogonal_polynomial(unsigned n, const Family& family, bool anti_gauss = false) : N(n), anti(anti_gauss)
   {
      values = calculate_values(family);
   }

   // Laurie 1996, Theorem 1: the anti-Gauss rule with n nodes has no node at or beyond the end e of the support
   // iff p~_n(e) has the sign of p_n(e), with p~_n defined as in strum_counter. Allow a margin, so that nodes
   // on the endpoint, which Chebyshev weights put there exactly, pass despite rounding; f is defined there.
   // An overflowed recurrence passes too, so that the rule reports the overflow.
   static bool anti_gauss_is_internal(unsigned n, const Family& family, const Real& e)
   {
      Real p = family.p0();
      Real p_previous = 0;
      for (unsigned k = 0; k + 1 < n; ++k)
      {
         Real p_next = family.next(e, p, p_previous, k);
         p_previous = p;
         p = p_next;
      }
      return !(family.next(e, p, 2 * p_previous, n - 1) / family.next(e, p, p_previous, n - 1) <= -4 * static_cast<Real>(n) * tools::epsilon<Real>());
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
   bool anti = false;
   std::pair<std::vector<Real>, std::vector<Real> > values;
   struct recurrence
   {
      Real p;
      Real p_prime;
      unsigned sign_changes;
      Real p_previous;
      Real p_prime_previous;
      Real p_double_prime;
   };

   class strum_counter
   {
      public:
         strum_counter(unsigned n, const Family& family_, bool anti_) : N(n), family(family_), anti(anti_) {}

         recurrence operator()(const Real& x) const
         {
            Real p = family.p0();
            Real p_previous = 0;
            // p~_N satisfies no differential equation, so differentiate the recurrence instead.
            Real p_prime = 0, p_prime_previous = 0, p_double_prime = 0, p_double_prime_previous = 0;
            bool last_negative = false;
            unsigned sign_changes = 0;
            for (unsigned k = 0; k < N; ++k)
            {
               // Laurie's anti-Gauss rule is the Gauss rule of the same recurrence with the coefficient of p_{N-2}
               // in its last step doubled: p~_N = p_N - b_{N-1} p_{N-2} in monic form. The Sturm count still holds.
               const Real c = (anti && (k + 1 == N)) ? Real(2) : Real(1);
               Real p_next = family.next(x, p, c * p_previous, k);
               if (anti)
               {
                  const Real slope = family.slope(k);
                  Real p_double_prime_next = family.next(x, p_double_prime, c * p_double_prime_previous, k) + 2 * slope * p_prime;
                  Real p_prime_next = family.next(x, p_prime, c * p_prime_previous, k) + slope * p;
                  p_double_prime_previous = p_double_prime;
                  p_double_prime = p_double_prime_next;
                  p_prime_previous = p_prime;
                  p_prime = p_prime_next;
               }
               p_previous = p;
               p = p_next;
               // A zero p_k(x) lies between values of opposite sign, so skipping it keeps the count right.
               if ((p != 0) && ((p < 0) != last_negative))
               {
                  ++sign_changes;
                  last_negative = !last_negative;
               }
            }
            if (anti)
               return recurrence{ p, p_prime, sign_changes, p_previous, p_prime_previous, p_double_prime };
            return recurrence{ p, family.derivative(x, p, p_previous), sign_changes, p_previous, 0, 0 };
         }

         Real second_derivative_ratio(const Real& x, const recurrence& r) const
         {
            return anti ? r.p_double_prime / r.p_prime : family.second_derivative_ratio(x, r.p, r.p_prime);
         }

         Real weight(const Real& z, const recurrence& r) const
         {
            // The anti-Gauss rule is the Gauss rule of the functional 2I - G_{N-1}, so for monic P its weights are
            // 2 ||P_{N-1}||^2 / (P_{N-1} P~'_N). The family's normalization changes them by a constant factor,
            // which the scaling to area() below removes. Like the family weights, evaluate them at the exact zero
            // to first order: near the ends of the interval they are sensitive to the rounding of z.
            if (!anti)
               return family.weight(z, r.p, r.p_prime);
            const Real delta = -r.p / r.p_prime;
            return 1 / ((r.p_previous + r.p_prime_previous * delta) * (r.p_prime + r.p_double_prime * delta));
         }

         // Replaces r.p, and for the Gauss rule r.p_previous and r.p_prime, with values to about twice the working
         // precision where the family can provide them; returns false where it cannot.
         template <class F = Family>
         auto compensated(const Real& z, recurrence& r, int) const -> decltype(std::declval<const F&>().compensated_value(z, 0u, false, r.p), bool())
         {
            Real p_previous;
            const Real p = family.compensated_value(z, N, anti, p_previous);
            if (!(BOOST_MATH_NAMESPACE::isfinite)(p) || !(BOOST_MATH_NAMESPACE::isfinite)(p_previous))
               return false;
            r.p = p;
            if (!anti)
            {
               r.p_previous = p_previous;
               r.p_prime = family.derivative(z, p, p_previous);
            }
            return true;
         }

         bool compensated(const Real&, recurrence&, long) const
         {
            return false;
         }

      private:
         unsigned N;
         Family family;
         bool anti;
   };

   // The zeros of p_N in increasing order, the nonnegative ones only for a family symmetric about 0, and their weights;
   // where glaser_liu_rokhlin cannot vouch for them, from the family's recurrence. By Sturm's theorem the number of sign
   // changes in p_0(x), ..., p_N(x) changes by one at each zero of p_N, in a direction that depends on the family:
   // bisecting on that count isolates each zero, and Halley's method then converges to it inside its bracket.
   std::pair<std::vector<Real>, std::vector<Real> > calculate_values(const Family& family)
   {
      if (N == 0)
      {
         static const char* function = "boost::math::detail::orthogonal_polynomial<%1%>::orthogonal_polynomial";
         BOOST_MATH_NAMESPACE::policies::raise_domain_error(function, "The argument n used for Gauss quadrature must be > 0 (got n=%1%).", static_cast<Real>(N), Policy());
         std::vector<Real> x = { std::numeric_limits<Real>::quiet_NaN() };
         std::vector<Real> w = { std::numeric_limits<Real>::quiet_NaN() };
         return std::make_pair(x, w);
      }
      using std::abs;
      using std::sqrt;
      const strum_counter evaluate(N, family, anti);

      const unsigned num_roots = family.num_roots();
      unsigned interior_roots = num_roots;
      std::vector<Real> x(num_roots), w(num_roots);
      {
         std::vector<Real> xg, wg;
         if (glaser_liu_rokhlin<Real, Family>(family, N, anti).solve(xg, wg))
            return scaled(family, xg, wg);
      }

      Real upper = family.upper_bound();
      Real lower = family.lower_bound();
      // All zeros lie below upper_bound(), so none lie above `upper`. Every p_k grows with x beyond
      // its zeros, so if the recurrence overflows anywhere we need it, it overflows here: report
      // that with NaNs, which turns into an evaluation error.
      recurrence top = evaluate(upper);
      if (!(BOOST_MATH_NAMESPACE::isfinite)(top.p) || !(BOOST_MATH_NAMESPACE::isfinite)(top.p_prime))
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
         w[0] = evaluate.weight(lower, bottom);
         interior_roots--;
      }
      // Need to know whether the number of sign changes increases or decreases with x, to know which direction to bisect.
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
         // With u = p/p' and v = p''/p' from the family, the step is u/(1 - uv/2); the ratios cannot
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
            Real v = evaluate.second_derivative_ratio(z, r);
            Real next = z - u / (1 - u * v / 2);
            if (!((next > low) && (next < high)))
               next = (low + high) / 2;
            const bool converged = abs(next - z) <= 2 * tools::epsilon<Real>() * abs(next);
            z = next;
            r = evaluate(z);
            if (converged)
               break;
         }
         // Halley's method has the zero of p_N as rounded; where the family can evaluate p_N more accurately,
         // Newton's method on that value finds the zero itself, and the weight uses those values too.
         if (evaluate.compensated(z, r, 0))
         {
            for (int i = 0; i < 2; ++i)
            {
               // Not the bracket: it came from the signs of the rounded p_N, which are wrong this close to the zero.
               const Real next = z - r.p / r.p_prime;
               if (!(abs(next - z) <= sqrt(tools::epsilon<Real>()) * abs(z)) || (next == z))
                  break;
               z = next;
               r = evaluate(z);
               evaluate.compensated(z, r, 0);
            }
         }
         x[num_roots - t] = z;
         w[num_roots - t] = evaluate.weight(z, r);
         upper = low;
      }
      return scaled(family, x, w);
   }

   // The weights share a relative error of a few eps from the rounded coefficients, growing with N. The exact weights
   // sum to area(), so scale them to that: weights with relative errors e_i then have e_i - e, e their weighted mean,
   // plus two roundings. That removes the shared part, which is what integrals of smooth functions see, and moves no
   // weight by more than |e| <= max |e_i|. Sum from the smallest weight, with Neumaier's compensation:
   static std::pair<std::vector<Real>, std::vector<Real> > scaled(const Family& family, std::vector<Real>& x, std::vector<Real>& w)
   {
      BOOST_MATH_NAMESPACE::detail::compensated_sum<Real> sum;
      for (std::size_t i = w.size(); i-- > 0;)
      {
         sum.add(Real(family.weight_multiplicity(static_cast<unsigned>(i))) * w[i]);
      }
      const Real scale = family.area() / sum.value();
      for (Real& weight_i : w)
      {
         weight_i *= scale;
      }
      return std::make_pair(x, w);
   }
};

} // namespace detail
BOOST_MATH_NAMESPACE_END

#endif // BOOST_MATH_DETAIL_ORTHOGONAL_POLYNOMIAL_HPP
