//  Copyright Nick Thompson 2026.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_QUADRATURE_LEJA_HPP
#define BOOST_MATH_QUADRATURE_LEJA_HPP

#ifdef _MSC_VER
#pragma once
#endif

#include <boost/math/tools/config.hpp>

#ifdef BOOST_MATH_NO_CXX17_IF_CONSTEXPR
#error "The header <boost/math/quadrature/leja.hpp> requires C++17 or later."
#endif

#ifndef BOOST_MATH_BUILD_MODULE
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <type_traits>
#include <utility>
#include <vector>
#ifdef BOOST_MATH_HAS_THREADS
#include <atomic>
#include <mutex>
#endif
#endif
#include <boost/math/constants/constants.hpp>
#include <boost/math/policies/error_handling.hpp>
#include <boost/math/policies/policy.hpp>
#include <boost/math/quadrature/gauss_jacobi.hpp>
#include <boost/math/quadrature/gauss_laguerre.hpp>
#include <boost/math/special_functions/fpclassify.hpp>
#include <boost/math/tools/compensated_arithmetic.hpp>
#include <boost/math/tools/precision.hpp>

BOOST_MATH_NAMESPACE_BEGIN namespace quadrature { namespace detail {

// The most points a rule computes. It keeps every rule, O(N^2) numbers, and takes O(N^2) operations to compute them:
// about 60 ms for 1023 in double precision, and 0.4 s for 127 in cpp_bin_float_50.
template <class Real>
constexpr unsigned leja_max_points = std::is_floating_point<Real>::value ? 1023 : 127;

// hi + lo, to twice the precision of Real, for the few quantities the weights cannot afford to round.
template <class Real>
using leja_double = BOOST_MATH_NAMESPACE::detail::double_word<Real>;
template <class Real>
using leja_arithmetic = BOOST_MATH_NAMESPACE::detail::double_word_arithmetic<Real>;

template <class Real>
leja_double<Real> leja_ldexp(const leja_double<Real>& a, int e)
{
   using std::ldexp;
   return {ldexp(a.hi, e), ldexp(a.lo, e)};
}

// A weight w on [lower(), upper()] is described by the total mass mu0, the coefficients of the recurrence
// x pi_k = pi_{k+1} + alpha_k pi_k + beta_k^2 pi_{k-1} of its monic orthogonal polynomials, to twice the precision of Real,
// and log sqrt(u) and its derivatives, where u places the points. u = w, except that a negative exponent of w at an end of the support, where
// sqrt(w) has no maximum, is 0 in u. Every exponent must be finite and > -1; exponent() is one that is not, if there is one.
// The points start from the mean of w, the one-point Gauss rule, not from 0, where x and x^2 vanish and a sparse grid
// of x_1 x_2 would see only zeros.
template <class Real>
struct leja_legendre_weight
{
   static constexpr const char* function = "boost::math::quadrature::leja_legendre<%1%>";
   static Real lower() { return -1; }
   static Real upper() { return 1; }
   static Real exponent() { return 0; }
   static std::vector<Real> start() { return {Real(0), Real(-1), Real(1)}; }
   static Real mu0() { return 2; }
   static Real half_log(const Real&) { return 0; }
   static Real half_log_derivative(const Real&) { return 0; }
   static Real half_log_second_derivative(const Real&) { return 0; }
   static leja_double<Real> alpha(unsigned, const leja_arithmetic<Real>&) { return {0, 0}; }
   static leja_double<Real> beta_squared(unsigned k, const leja_arithmetic<Real>& op)
   {
      return op.divide(leja_double<Real>{Real(k) * k, 0}, leja_double<Real>{Real(4) * k * k - 1, 0});
   }
   // The factor h |h|^(alpha + beta) by which x = c + h t takes the integral over [-1, 1] to one over [c - h, c + h].
   static Real jacobian(const Real& h) { return h; }
};

template <class Real>
struct leja_hermite_weight
{
   static constexpr const char* function = "boost::math::quadrature::leja_hermite<%1%>";
   static Real lower() { return -std::numeric_limits<Real>::infinity(); }
   static Real upper() { return std::numeric_limits<Real>::infinity(); }
   static Real exponent() { return 0; }
   static std::vector<Real> start() { return {Real(0)}; }
   static Real mu0() { return BOOST_MATH_NAMESPACE::constants::root_pi<Real>(); }
   static Real half_log(const Real& x) { return -x * x / 2; }
   static Real half_log_derivative(const Real& x) { return -x; }
   static Real half_log_second_derivative(const Real&) { return -1; }
   static leja_double<Real> alpha(unsigned, const leja_arithmetic<Real>&) { return {0, 0}; }
   static leja_double<Real> beta_squared(unsigned k, const leja_arithmetic<Real>&) { return {Real(k) / 2, 0}; }
};

// x^a exp(-x) on [0, infinity).
template <class Real>
struct leja_laguerre_weight
{
   static constexpr const char* function = "boost::math::quadrature::leja_laguerre<%1%>";
   explicit leja_laguerre_weight(const Real& alpha_ = 0) : a(alpha_), u((std::max)(alpha_, Real(0))) {}
   static Real lower() { return 0; }
   static Real upper() { return std::numeric_limits<Real>::infinity(); }
   Real exponent() const { return a; }
   std::vector<Real> start() const { return {a + 1}; }
   Real mu0() const { return BOOST_MATH_NAMESPACE::detail::laguerre_family<Real>(0, a).area(); }
   Real half_log(const Real& x) const
   {
      using std::log;
      return (u > 0 ? u * log(x) / 2 : Real(0)) - x / 2;
   }
   Real half_log_derivative(const Real& x) const { return (u > 0 ? u / (2 * x) : Real(0)) - Real(1) / 2; }
   Real half_log_second_derivative(const Real& x) const { return u > 0 ? -u / (2 * x * x) : Real(0); }
   leja_double<Real> alpha(unsigned k, const leja_arithmetic<Real>& op) const { return op.add(D{Real(2 * k + 1), 0}, a); }
   leja_double<Real> beta_squared(unsigned k, const leja_arithmetic<Real>& op) const { return op.multiply(op.add(D{Real(k), 0}, a), Real(k)); }
   using D = leja_double<Real>;

   Real a, u;
};

// (1 - x)^a (1 + x)^b on [-1, 1].
template <class Real>
struct leja_jacobi_weight
{
   static constexpr const char* function = "boost::math::quadrature::leja_jacobi<%1%>";
   leja_jacobi_weight(const Real& alpha_, const Real& beta_)
      : a(alpha_), b(beta_), ua((std::max)(alpha_, Real(0))), ub((std::max)(beta_, Real(0))) {}
   static Real lower() { return -1; }
   static Real upper() { return 1; }
   Real exponent() const { return a > -1 && (BOOST_MATH_NAMESPACE::isfinite)(a) ? b : a; }
   Real jacobian(const Real& h) const
   {
      using std::abs;
      using std::pow;
      return h * pow(abs(h), a + b);
   }
   std::vector<Real> start() const { return {(b - a) / (a + b + 2)}; }
   Real mu0() const { return BOOST_MATH_NAMESPACE::detail::jacobi_family<Real>(0, a, b).area(); }
   Real half_log(const Real& x) const
   {
      using std::log;
      return (ua > 0 ? ua * log(1 - x) / 2 : Real(0)) + (ub > 0 ? ub * log(1 + x) / 2 : Real(0));
   }
   Real half_log_derivative(const Real& x) const
   {
      return (ub > 0 ? ub / (2 * (1 + x)) : Real(0)) - (ua > 0 ? ua / (2 * (1 - x)) : Real(0));
   }
   Real half_log_second_derivative(const Real& x) const
   {
      return -(ua > 0 ? ua / (2 * (1 - x) * (1 - x)) : Real(0)) - (ub > 0 ? ub / (2 * (1 + x) * (1 + x)) : Real(0));
   }
   // https://dlmf.nist.gov/18.9; k = 0 and 1 apart, where 2k + a + b or 2k + a + b - 1 can vanish.
   leja_double<Real> alpha(unsigned k, const leja_arithmetic<Real>& op) const
   {
      const D b_minus_a = op.add(D{b, 0}, -a), s = op.add(op.add(D{Real(2 * k), 0}, a), b), s2 = op.add(s, Real(2));
      if (k == 0)
         return op.divide(b_minus_a, s2);
      return op.divide(op.multiply(b_minus_a, op.add(D{b, 0}, a)), op.multiply(s, s2));
   }
   leja_double<Real> beta_squared(unsigned k, const leja_arithmetic<Real>& op) const
   {
      const D s = op.add(op.add(D{Real(2 * k), 0}, a), b), s_s1 = op.multiply(op.multiply(s, s), op.add(s, Real(1)));
      if (k == 1)
         return op.divide(op.multiply(op.multiply(op.add(D{1, 0}, a), Real(4)), op.add(D{1, 0}, b)), s_s1);
      const D ka = op.add(D{Real(k), 0}, a), kb = op.add(D{Real(k), 0}, b);
      const D numerator = op.multiply(op.multiply(op.multiply(ka, Real(4 * k)), kb), op.add(ka, b));
      return op.divide(numerator, op.multiply(s_s1, op.add(s, Real(-1))));
   }
   using D = leja_double<Real>;

   Real a, b, ua, ub;
};

// A vector-valued integrand returns a random-access container of reals or complex numbers, such as std::array or std::vector.
template <class K, class = void>
struct leja_is_vector : std::false_type {};
template <class K>
struct leja_is_vector<K, std::void_t<decltype(std::declval<K&>()[0]), decltype(std::declval<K&>().size())>> : std::true_type {};

// |a - b|, or |a| without b, in the max norm for a vector; NaN if any component is NaN.
template <class Real, class K>
Real leja_distance(const K& a, const K* b = nullptr)
{
   using std::abs;
   if constexpr (leja_is_vector<K>::value)
   {
      Real m = 0;
      for (std::size_t j = 0; j < a.size(); ++j)
      {
         const Real d = b ? Real(abs(a[j] - (*b)[j])) : Real(abs(a[j]));
         if (d != d)
            return d;
         m = (std::max)(m, d);
      }
      return m;
   }
   else
      return b ? Real(abs(a - *b)) : Real(abs(a));
}

// NaN, in every component of a default-constructed vector (none for std::vector).
template <class K, class Real>
K leja_nan()
{
   if constexpr (leja_is_vector<K>::value)
   {
      K y{};
      for (std::size_t j = 0; j < y.size(); ++j)
         y[j] = std::numeric_limits<Real>::quiet_NaN();
      return y;
   }
   else
      return K(std::numeric_limits<Real>::quiet_NaN());
}

// y = a x if first, else y += a x; a vector y takes its size from x.
template <class K, class Real>
void leja_accumulate(K& y, const Real& a, const K& x, bool first)
{
   if constexpr (leja_is_vector<K>::value)
   {
      if (first)
         y = x;
      for (std::size_t j = 0; j < y.size(); ++j)
      {
         if (first)
            y[j] = a * x[j];
         else
            y[j] += a * x[j];
      }
   }
   else if (first)
      y = a * x;
   else
      y += a * x;
}

template <class Real, class Weight, class Policy>
class leja_detail
{
public:
   using value_type = Real;
   using policy_type = Policy;

   // The arguments are the weight's: alpha for leja_laguerre, alpha and beta for leja_jacobi. The points and weights
   // are computed as they are needed, by any thread, up to max_points().
   template <class... Args, typename std::enable_if<std::is_constructible<Weight, const Args&...>::value, bool>::type = true>
   explicit leja_detail(const Args&... args) : m_weight(args...)
   {
      reserve();
      if (!(m_weight.exponent() > -1 && (BOOST_MATH_NAMESPACE::isfinite)(m_weight.exponent())))
         policies::raise_domain_error(Weight::function, "The exponents of the weight must be finite and > -1, but got %1%.", m_weight.exponent(), Policy());
   }

   leja_detail(const leja_detail& rule) : m_weight(rule.m_weight)
   {
#ifdef BOOST_MATH_HAS_THREADS
      std::lock_guard<std::mutex> lock(rule.m_mutex);
#endif
      reserve();
      const unsigned size = rule.m_size;
      m_abscissa.assign(rule.m_abscissa.begin(), rule.m_abscissa.begin() + size);
      m_weights.assign(rule.m_weights.begin(), rule.m_weights.begin() + size + 1);
      m_state = rule.m_state;
      m_size = size;
   }

   leja_detail& operator=(const leja_detail& rule)
   {
      leja_detail copy(rule);
      m_weight = copy.m_weight;
      m_abscissa.swap(copy.m_abscissa);
      m_weights.swap(copy.m_weights);
      std::swap(m_state, copy.m_state);
      m_size = static_cast<unsigned>(copy.m_size);
      return *this;
   }

   // The points computed so far, in Leja order; the n-point rule uses the first n, which weights(n) computes. Those
   // another thread is computing are NaN until it is done.
   const std::vector<Real>& abscissa() const
   {
      return m_abscissa;
   }

   const std::vector<Real>& weights(unsigned n) const
   {
      if (n == 0 || n > max_points())
      {
         policies::raise_domain_error(Weight::function, "The number of points must be in [1, max_points()], but got %1%.", Real(n), Policy());
         return m_weights[0];
      }
      if (n > m_size)
         extend(n);
      // Empty for a weight with an invalid exponent, after a domain error the policy ignored.
      return n <= m_size ? m_weights[n] : m_weights[0];
   }

   static constexpr unsigned max_points()
   {
      return leja_max_points<Real>;
   }

   // Adds one point at a time until |Q_m - Q_{m - max(1, m/4)}| is small for every m in [n - n/4, n]: a single difference,
   // and |Q_n - Q_{n-1}| even more so, often underestimates the error of slowly converging integrands. Like tanh_sinh,
   // it returns its best result if max_points() is exhausted first.
   template <class F>
   auto integrate(F f, Real tolerance = tools::root_epsilon<Real>(), Real* error = nullptr, Real* pL1 = nullptr) const
      ->decltype(std::declval<F>()(std::declval<Real>()))
   {
      using K = decltype(f(Real(0)));
      static_assert(!std::is_integral<K>::value,
                    "The return type cannot be integral, it must be a real or complex floating point type, or a vector of them.");
      std::vector<K> y, Q;
      // NaN for an empty rule.
      K result = leja_nan<K, Real>();
      Real L1 = std::numeric_limits<Real>::quiet_NaN();
      Real err = L1;
      std::vector<Real> e; // e[n - 1] = |Q_n - Q_{n - max(1, n/4)}|
      for (unsigned n = 1; n <= max_points() && !weights(n).empty(); ++n)
      {
         y.push_back(f(m_abscissa[n - 1]));
         result = sum(y, n, L1);
         Q.push_back(result);
         e.push_back(n > 1 ? leja_distance<Real>(result, &Q[n - 1 - (std::max)(1u, n / 4)]) : L1);
         err = *std::max_element(e.end() - (n / 4 + 1), e.end());
         // A NaN from f: max_element need not return it.
         if (L1 != L1)
         {
            err = L1;
            break;
         }
         if (n >= 4 && err <= tolerance * L1)
            break;
      }
      if (error)
         *error = err;
      if (pL1)
         *pL1 = L1;
      return result;
   }

   // The integral over [a, b], through x = (a + b)/2 + h t, h = (b - a)/2: of f for leja_legendre, and of
   // f(x) |b - x|^alpha |x - a|^beta for leja_jacobi. For a > b it is minus the integral over [b, a].
   template <class F, class W = Weight, class = decltype(std::declval<const W&>().jacobian(Real(0)))>
   auto integrate(F f, const Real& a, const Real& b, Real tolerance = tools::root_epsilon<Real>(), Real* error = nullptr,
                  Real* pL1 = nullptr) const->decltype(std::declval<F>()(std::declval<Real>()))
   {
      using std::abs;
      using K = decltype(f(Real(0)));
      const Real h = (b - a) / 2;
      if (!((BOOST_MATH_NAMESPACE::isfinite)(a) && (BOOST_MATH_NAMESPACE::isfinite)(b)) || h == 0)
      {
         const Real r = h == 0 ? Real(0) : policies::raise_domain_error(Weight::function,
            "The bounds must be finite, but got %1%; leja_hermite and leja_laguerre integrate over infinite intervals.",
            (BOOST_MATH_NAMESPACE::isfinite)(a) ? b : a, Policy());
         if (error)
            *error = r;
         if (pL1)
            *pL1 = r;
         // 0, of the shape of f(a), over an empty interval.
         K result = leja_nan<K, Real>();
         if (h == 0)
            leja_accumulate(result, Real(0), f(a), true);
         return result;
      }
      const Real c = (a + b) / 2, J = jacobian(h);
      K result;
      leja_accumulate(result, J, integrate([&](const Real& t) { return f(c + h * t); }, tolerance, error, pL1), true);
      if (error)
         *error *= abs(J);
      if (pL1)
         *pL1 *= abs(J);
      return result;
   }

   // The factor by which the map from [-1, 1] to an interval of half-width h, of the sign of its orientation, scales an
   // integral; leja_legendre and leja_jacobi only.
   template <class W = Weight>
   auto jacobian(const Real& h) const->decltype(std::declval<const W&>().jacobian(h))
   {
      return m_weight.jacobian(h);
   }

private:
   void reserve()
   {
      m_abscissa.reserve(max_points());
      m_weights.reserve(max_points() + 1);
      m_weights.resize(1);
   }

   // Computes the points and weights of the rules up to n points, which the other threads see once m_size is n. The
   // vectors never reallocate, so the references the others hold stay valid.
   void extend(unsigned n) const
   {
#ifdef BOOST_MATH_HAS_THREADS
      std::lock_guard<std::mutex> lock(m_mutex);
#endif
      const unsigned size = m_size;
      if (n <= size || !(m_weight.exponent() > -1 && (BOOST_MATH_NAMESPACE::isfinite)(m_weight.exponent())))
         return;
      std::vector<Real> x(m_abscissa.begin(), m_abscissa.begin() + size);
      m_abscissa.resize(n, std::numeric_limits<Real>::quiet_NaN());
      compute_abscissa(x, n);
      std::copy(x.begin() + size, x.end(), m_abscissa.begin() + size);
      compute_weights(size, n);
      m_size = n;
   }

   template <class K>
   K sum(const std::vector<K>& y, unsigned n, Real& L1) const
   {
      using std::abs;
      K result{};
      L1 = 0;
      for (unsigned i = 0; i < n; ++i)
      {
         leja_accumulate(result, m_weights[n][i], y[i], i == 0);
         L1 += abs(m_weights[n][i]) * leja_distance<Real>(y[i]);
      }
      return result;
   }

   // F(t) = log sqrt(u(t)) + sum_k log|t - x_k| is concave between consecutive points, so each gap has one maximum,
   // found by Newton's method on F' safeguarded by bisection.
   Real derivatives(const std::vector<Real>& points, const Real& t, Real& second) const
   {
      Real first = m_weight.half_log_derivative(t);
      second = m_weight.half_log_second_derivative(t);
      for (const Real& x : points)
      {
         const Real r = 1 / (t - x);
         first += r;
         second -= r * r;
      }
      return first;
   }

   // The points between l and r, or an end of the support where end_l or end_r. lower = F(t) is at most the maximum of
   // F there, and upper at least; they are equal, at the maximum t, after a search.
   struct gap
   {
      Real l, r, t, lower, upper, scale;
      bool end_l, end_r;
   };

   // Sets t to the maximum of F on the gap, starting from t0.
   void search(const std::vector<Real>& points, gap& g, Real t) const
   {
      using std::abs;
      using std::frexp;
      using std::log;
      const Real inf = std::numeric_limits<Real>::infinity();
      Real l = g.l, r = g.r, second;
      // At a finite end of the support, u > 0 can make the end itself the maximum.
      const bool at_l = g.end_l && l > -inf && derivatives(points, l, second) <= 0;
      if (at_l || (g.end_r && r < inf && derivatives(points, r, second) >= 0))
         t = at_l ? l : r;
      else
      {
         if (l == -inf || r == inf)
         {  // Double the ray's length until F' changes sign.
            const Real end = l == -inf ? r : l;
            const Real direction = l == -inf ? -1 : 1;
            Real step = 1;
            while (direction * derivatives(points, end + direction * step, second) > 0)
               step *= 2;
            l = end + direction * step;
            r = end + direction * step / 2;
            if (step == 1)
               r = end;
            if (direction > 0)
               std::swap(l, r);
         }
         // Relative to the gap, not to t, which can be 0.
         const Real tolerance = 2 * tools::epsilon<Real>() * (std::max)(abs(l), abs(r));
         if (!(t > l && t < r))
            t = (l + r) / 2;
         for (unsigned i = 0; i < 4 * tools::digits<Real>(); ++i)
         {
            const Real first = derivatives(points, t, second);
            (first > 0 ? l : r) = t;
            Real next = t - first / second;
            if (abs(next - t) <= tolerance)
            {
               t = next;
               break;
            }
            if (!(next > l && next < r))
               next = (l + r) / 2;
            t = next;
            if (r - l <= 2 * tolerance)
               break;
         }
      }
      // F(t), with one logarithm: the product of the |t - x_k| is P 2^e, renormalized every four factors. Each factor
      // adds at most eps to its rounding error, and log rounds |log P + e ln 2| and |log sqrt(u)| with relative error eps.
      Real P = 1;
      int e = 0;
      for (std::size_t k = 0; k < points.size(); k += 4)
      {
         Real block = abs(t - points[k]);
         for (std::size_t j = k + 1; j < (std::min)(k + 4, points.size()); ++j)
            block *= abs(t - points[j]);
         int exponent;
         P = frexp(P * block, &exponent);
         e += exponent;
      }
      const Real half_log = m_weight.half_log(t);
      const Real product = log(P) + e * constants::ln_two<Real>();
      g.t = t;
      g.lower = g.upper = half_log + product;
      g.scale = abs(half_log) + abs(product) + 2 * points.size();
   }

   // Adds points to points until it has N. A new point x adds log|t - x| to F, so a gap's maximum grows by at most the
   // largest log|t - x| over the gap, and by at least log|t - x| at the last maximum: only the gaps whose upper bound
   // reaches the largest lower bound need a search, about 8 of them for each point up to 1023.
   void compute_abscissa(std::vector<Real>& points, unsigned N) const
   {
      using std::abs;
      using std::log;
      const Real eps = tools::epsilon<Real>();
      const Real inf = std::numeric_limits<Real>::infinity();
      const std::vector<Real> start = m_weight.start();
      while (points.size() < (std::min)(static_cast<std::size_t>(N), start.size()))
         points.push_back(start[points.size()]);
      std::vector<gap>& gaps = m_state.gaps;
      if (points.size() < N && gaps.empty())
      {
         std::vector<Real> sorted = points;
         std::sort(sorted.begin(), sorted.end());
         for (std::size_t g = 0; g <= sorted.size(); ++g)
         {
            gap G{g == 0 ? m_weight.lower() : sorted[g - 1], g == sorted.size() ? m_weight.upper() : sorted[g], 0, 0, 0, 0,
                  g == 0, g == sorted.size()};
            if (G.l < G.r)
            {
               search(points, G, (G.l + G.r) / 2);
               gaps.push_back(G);
            }
         }
      }
      while (points.size() < N)
      {
         // The bounds gain a rounding error with each point.
         Real best = -inf, scale = 0;
         for (const gap& g : gaps)
         {
            best = (std::max)(best, g.lower);
            scale = (std::max)(scale, g.scale);
         }
         const Real slack = 4 * (points.size() + 8) * eps * (abs(best) + scale);
         std::vector<std::size_t> searched;
         for (std::size_t i = 0; i < gaps.size(); ++i)
         {
            if (gaps[i].upper < best - slack)
               continue;
            if (gaps[i].upper > gaps[i].lower)
               search(points, gaps[i], gaps[i].t);
            searched.push_back(i);
         }
         // Maxima within 8 eps scale, a few times the rounding error of F, are tied, and the leftmost wins. Exact ties
         // are common (+-1/sqrt(3) for Legendre, +-1 for Hermite), and this way every precision picks the same point.
         std::size_t top = searched[0];
         for (std::size_t i : searched)
            if (gaps[i].lower > gaps[top].lower)
               top = i;
         std::size_t c = 0;
         while (std::find(searched.begin(), searched.end(), c) == searched.end()
                || gaps[c].lower < gaps[top].lower - 8 * eps * gaps[top].scale)
            ++c;
         const Real x = gaps[c].t;
         for (gap& g : gaps)
         {
            g.lower += log(abs(g.t - x));
            g.upper += g.l == -inf || g.r == inf ? inf : log((std::max)(abs(g.l - x), abs(g.r - x)));
         }
         points.push_back(x);
         gap left{gaps[c].l, x, 0, 0, 0, 0, gaps[c].end_l, false}, right{x, gaps[c].r, 0, 0, 0, 0, false, gaps[c].end_r};
         gaps.erase(gaps.begin() + c);
         for (gap* g : {&right, &left})
         {
            if (g->l < g->r)
            {
               search(points, *g, (g->l + g->r) / 2);
               gaps.insert(gaps.begin() + c, *g);
            }
         }
      }
   }

   // The rules in Newton form: Q_{n+1} f = Q_n f + I[omega_n] f[x_0, ..., x_n], where omega_n = prod_{k<n} (x - x_k) and
   // f[x_0, ..., x_n] = sum_i f(x_i) / d_i, d_i = prod_{j<=n, j!=i} (x_i - x_j), so each rule costs O(n). I[omega_n] is
   // mu0 c_0, from the coefficients c of omega_n in the monic orthogonal polynomials scaled by 2^-e_k, e_k ~ log2 of the
   // norm of pi_k, which keeps them in range. It cancels by up to 1e16, so c, d and the weights carry twice the precision
   // of Real; c, and each d_i, have a power of 2 of their own.
   void compute_weights(unsigned old, unsigned N) const
   {
      using std::abs;
      using std::floor;
      using std::frexp;
      using std::log;
      using D = leja_double<Real>;
      const leja_arithmetic<Real> op;
      state& f = m_state;
      const Real mu0 = m_weight.mu0();
      auto normalize = [](D& v, int& e) {
         int k;
         frexp(v.hi, &k);
         v = leja_ldexp(v, -k);
         e += k;
      };
      for (unsigned n = old; n < N; ++n)
      {
         const Real& x = m_abscissa[n];
         if (n == 0)
         {
            f.c = {D{1, 0}};
            f.e = {0};
            f.d = {D{1, 0}};
            f.d_exponent = {0};
            f.w = {D{mu0, 0}};
            m_weights.push_back({mu0});
            continue;
         }
         // omega_n = (x - x_{n-1}) omega_{n-1}.
         f.log2_norm += log(m_weight.beta_squared(n, op).hi) / (2 * constants::ln_two<Real>());
         f.e.push_back(static_cast<int>(floor(f.log2_norm + Real(0.5))));
         std::vector<D> c(n + 1, D{0, 0});
         Real largest = 0;
         for (unsigned k = 0; k <= n; ++k)
         {
            if (k > 0)
               c[k] = leja_ldexp(f.c[k - 1], f.e[k] - f.e[k - 1]);
            if (k < n)
               c[k] = op.add(c[k], op.multiply(op.add(m_weight.alpha(k, op), -m_abscissa[n - 1]), f.c[k]));
            if (k + 1 < n)
               c[k] = op.add(c[k], leja_ldexp(op.multiply(m_weight.beta_squared(k + 1, op), f.c[k + 1]), f.e[k] - f.e[k + 1]));
            largest = (std::max)(largest, Real(abs(c[k].hi)));
         }
         int e;
         frexp(largest, &e);
         for (D& v : c)
            v = leja_ldexp(v, -e);
         f.c_exponent += e;
         f.c.swap(c);
         // d_i gains the factor x_i - x, and d_n is new.
         D d{1, 0};
         int d_exponent = 0;
         for (unsigned i = 0; i < n; ++i)
         {
            f.d[i] = op.multiply(f.d[i], op.two_sum(m_abscissa[i], -x));
            normalize(f.d[i], f.d_exponent[i]);
            d = op.multiply(d, op.two_sum(x, -m_abscissa[i]));
            normalize(d, d_exponent);
         }
         f.d.push_back(d);
         f.d_exponent.push_back(d_exponent);
         f.w.push_back(D{0, 0});
         const D I = op.multiply(f.c[0], mu0);
         std::vector<Real> w(n + 1);
         for (unsigned i = 0; i <= n; ++i)
         {
            f.w[i] = op.add(f.w[i], leja_ldexp(op.divide(I, f.d[i]), f.c_exponent - f.d_exponent[i]));
            w[i] = f.w[i].hi + f.w[i].lo;
         }
         m_weights.push_back(std::move(w));
      }
   }

   // The Newton form of the rules so far, and the gaps between the points, so that a rule computes the same points and
   // weights however it grows.
   struct state
   {
      std::vector<leja_double<Real>> c, d, w;
      std::vector<int> e, d_exponent;
      int c_exponent = 0;
      Real log2_norm = 0;
      std::vector<gap> gaps;
   };

   Weight m_weight;
   mutable std::vector<Real> m_abscissa;
   mutable std::vector<std::vector<Real>> m_weights; // m_weights[n] is the n-point rule; m_weights[0] is empty.
   mutable state m_state;
#ifdef BOOST_MATH_HAS_THREADS
   mutable std::mutex m_mutex;
   mutable std::atomic<unsigned> m_size{0};
#else
   mutable unsigned m_size = 0;
#endif
};

} // namespace detail

// Nested Leja quadrature for the integral of f over [-1, 1], the Legendre weight w = 1.
template <class Real, class Policy = BOOST_MATH_NAMESPACE::policies::policy<> >
using leja_legendre = detail::leja_detail<Real, detail::leja_legendre_weight<Real>, Policy>;

// Nested weighted Leja quadrature for the integral of f(x) exp(-x^2) over the real line.
template <class Real, class Policy = BOOST_MATH_NAMESPACE::policies::policy<> >
using leja_hermite = detail::leja_detail<Real, detail::leja_hermite_weight<Real>, Policy>;

// Nested weighted Leja quadrature for the integral of f(x) x^alpha exp(-x) over [0, infinity), alpha > -1:
// leja_laguerre<Real>(alpha = 0).
template <class Real, class Policy = BOOST_MATH_NAMESPACE::policies::policy<> >
using leja_laguerre = detail::leja_detail<Real, detail::leja_laguerre_weight<Real>, Policy>;

// Nested weighted Leja quadrature for the integral of f(x) (1 - x)^alpha (1 + x)^beta over [-1, 1], alpha, beta > -1:
// leja_jacobi<Real>(alpha, beta).
template <class Real, class Policy = BOOST_MATH_NAMESPACE::policies::policy<> >
using leja_jacobi = detail::leja_detail<Real, detail::leja_jacobi_weight<Real>, Policy>;

} // namespace quadrature
BOOST_MATH_NAMESPACE_END

#endif // BOOST_MATH_QUADRATURE_LEJA_HPP
