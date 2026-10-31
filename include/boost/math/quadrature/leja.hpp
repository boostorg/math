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
#include <boost/math/tools/compensated_arithmetic.hpp>
#include <boost/math/tools/precision.hpp>

BOOST_MATH_NAMESPACE_BEGIN namespace quadrature { namespace detail {

// The most points a rule computes. The weights are accurate well beyond, but the rules take O(N^3) to compute: about
// 3 s for 1023 in double precision, and 2.4 s for 127 in cpp_bin_float_50.
template <class Real>
constexpr unsigned leja_max_points = std::is_floating_point<Real>::value ? 1023 : 127;

// A weight w on [lower(), upper()] is described by the total mass mu0, the coefficients of the recurrence
// x p_k = beta_{k+1} p_{k+1} + alpha_k p_k + beta_k p_{k-1} of its orthonormal polynomials, and log sqrt(u) and its
// derivatives, where u places the points. u = w, except that a negative exponent of w at an end of the support, where
// sqrt(w) has no maximum, is 0 in u. Every exponent must be > -1; exponent() is one that is not, if there is one.
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
   static Real alpha(unsigned) { return 0; }
   static Real beta(unsigned k)
   {
      using std::sqrt;
      return k == 0 ? Real(0) : Real(k) / sqrt(Real(4) * k * k - 1);
   }
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
   static Real alpha(unsigned) { return 0; }
   static Real beta(unsigned k)
   {
      using std::sqrt;
      return sqrt(Real(k) / 2);
   }
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
   Real alpha(unsigned k) const { return 2 * k + 1 + a; }
   Real beta(unsigned k) const
   {
      using std::sqrt;
      return sqrt(k * (k + a));
   }

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
   Real exponent() const { return a > -1 ? b : a; }
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
   // https://dlmf.nist.gov/18.9, orthonormal; k = 0 and 1 apart, where 2k + a + b or 2k + a + b - 1 can vanish.
   Real alpha(unsigned k) const
   {
      const Real s = 2 * k + a + b;
      return k == 0 ? (b - a) / (a + b + 2) : (b - a) * (b + a) / (s * (s + 2));
   }
   Real beta(unsigned k) const
   {
      using std::sqrt;
      const Real s = 2 * k + a + b;
      if (k <= 1)
         return k == 0 ? Real(0) : 2 * sqrt((1 + a) * (1 + b) / (s + 1)) / s;
      return 2 * sqrt(k * (k + a) * (k + b) * (k + a + b) / ((s + 1) * (s - 1))) / s;
   }

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
      if (!(m_weight.exponent() > -1))
         policies::raise_domain_error(Weight::function, "The exponents of the weight must be > -1, but got %1%.", m_weight.exponent(), Policy());
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
      m_factorization = rule.m_factorization;
      m_size = size;
   }

   leja_detail& operator=(const leja_detail& rule)
   {
      leja_detail copy(rule);
      m_weight = copy.m_weight;
      m_abscissa.swap(copy.m_abscissa);
      m_weights.swap(copy.m_weights);
      std::swap(m_factorization, copy.m_factorization);
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

   // The n-point rule; integral n so that integrate(f, 10) does not compete with the adaptive overload.
   template <class F, class Integer, typename std::enable_if<std::is_integral<Integer>::value, bool>::type = true>
   auto integrate(F f, Integer n, Real* pL1 = nullptr) const->decltype(std::declval<F>()(std::declval<Real>()))
   {
      using K = decltype(f(Real(0)));
      static_assert(!std::is_integral<K>::value,
                    "The return type cannot be integral, it must be a real or complex floating point type, or a vector of them.");
      // Empty, after a domain error the policy ignored, unless 1 <= n <= max_points() and the weight is valid.
      const unsigned m = static_cast<unsigned>(weights(n > 0 ? static_cast<unsigned>(n) : 0u).size());
      std::vector<K> y;
      for (unsigned i = 0; i < m; ++i)
         y.push_back(f(m_abscissa[i]));
      Real L1 = std::numeric_limits<Real>::quiet_NaN();
      K result = m ? sum(y, m, L1) : leja_nan<K, Real>();
      if (pL1)
         *pL1 = L1;
      return result;
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
      if (n <= size || !(m_weight.exponent() > -1))
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

   // Adds points to points until it has N.
   void compute_abscissa(std::vector<Real>& points, unsigned N) const
   {
      using std::abs;
      using std::frexp;
      using std::log;
      const Real eps = tools::epsilon<Real>();
      const std::vector<Real> start = m_weight.start();
      while (points.size() < (std::min)(static_cast<std::size_t>(N), start.size()))
         points.push_back(start[points.size()]);
      std::vector<Real> sorted, argmax, value, scale, curvature;
      std::vector<Real>& previous = m_factorization.argmax;
      std::vector<Real>& previous_curvature = m_factorization.curvature;
      std::size_t& chosen = m_factorization.chosen;
      // F(t), with one logarithm: the product of the |t - x_k| is P 2^e, renormalized every four factors. Each factor
      // adds at most eps to its rounding error, and log rounds |log P + e ln 2| and |log sqrt(u)| with relative error eps.
      auto candidate = [&](const Real& t, const Real& second) {
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
         argmax.push_back(t);
         curvature.push_back(second);
         value.push_back(half_log + product);
         scale.push_back(abs(half_log) + abs(product) + 2 * points.size());
      };
      while (points.size() < N)
      {
         sorted = points;
         std::sort(sorted.begin(), sorted.end());
         argmax.clear();
         curvature.clear();
         value.clear();
         scale.clear();
         // Gap g lies between sorted[g - 1] and sorted[g]; gaps 0 and sorted.size() end at the ends of the support, and
         // are empty where a point is there.
         for (std::size_t g = 0; g <= sorted.size(); ++g)
         {
            Real l = g == 0 ? m_weight.lower() : sorted[g - 1];
            Real r = g == sorted.size() ? m_weight.upper() : sorted[g];
            Real second;
            if (!(l < r))
            {
               argmax.push_back(l);
               curvature.push_back(0);
               value.push_back(-std::numeric_limits<Real>::infinity());
               scale.push_back(0);
               continue;
            }
            // At a finite end of the support, u > 0 can make the end itself the maximum.
            if ((g == 0 && l > -std::numeric_limits<Real>::infinity() && derivatives(points, l, second) <= 0)
                || (g == sorted.size() && r < std::numeric_limits<Real>::infinity() && derivatives(points, r, second) >= 0))
            {
               candidate(g == 0 ? l : r, Real(0));
               continue;
            }
            if (l == -std::numeric_limits<Real>::infinity() || r == std::numeric_limits<Real>::infinity())
            {  // Double the ray's length until F' changes sign.
               const Real end = g == 0 ? sorted.front() : sorted.back();
               const Real direction = g == 0 ? -1 : 1;
               Real step = 1;
               while (direction * derivatives(points, end + direction * step, second) > 0)
                  step *= 2;
               l = end + direction * step;
               r = end + direction * step / 2;
               if (step == 1)
                  r = end;
               if (g != 0)
                  std::swap(l, r);
            }
            // Relative to the gap, not to t, which can be 0.
            const Real tolerance = 2 * eps * (std::max)(abs(l), abs(r));
            // The new point x split gap chosen in two. In the others, F gained log|t - x|, and one Newton step from the
            // last maximum s, with F'(s) = 0 and F''(s) known, costs nothing.
            Real t = (l + r) / 2;
            const std::size_t old = g < chosen ? g : g - 1;
            if (!previous.empty() && g != chosen && g != chosen + 1 && previous_curvature[old] < 0)
            {
               const Real d = 1 / (previous[old] - points.back());
               const Real next = previous[old] - d / (previous_curvature[old] - d * d);
               if (next > l && next < r)
                  t = next;
            }
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
            candidate(t, second);
         }
         // Maxima within 8 eps scale, a few times the rounding error of F, are tied, and the
         // leftmost wins. Exact ties are common (+-1/sqrt(3) for Legendre, +-1 for Hermite), and this way every
         // precision picks the same point.
         const std::size_t best = std::max_element(value.begin(), value.end()) - value.begin();
         std::size_t g = 0;
         while (value[g] < value[best] - 8 * eps * scale[best])
            ++g;
         points.push_back(argmax[g]);
         chosen = g;
         previous.swap(argmax);
         previous_curvature.swap(curvature);
      }
   }

   // The weights solve sum_i w_i p_k(x_i) = sqrt(mu0) delta_k0, k < n. Scaling column i by sqrt(u(x_i)) makes the
   // matrix A_ki = sqrt(u(x_i)) p_k(x_i) well conditioned where the Vandermonde matrix p_k(x_i) is not. The rule with n
   // points has the leading n x n block of A, so one LU factorization without pivoting serves every n: the Leja order
   // pivots it, since each point maximizes the next pivot. Each point adds a row of L and a column of U, and each rule a
   // back substitution, in O(N^3) for all of them.
   void compute_weights(unsigned old, unsigned N) const
   {
      using std::abs;
      using std::exp;
      using std::log;
      using std::sqrt;
      const Real mu0 = m_weight.mu0();
      const Real eps = tools::epsilon<Real>();
      factorization& f = m_factorization;
      const BOOST_MATH_NAMESPACE::detail::two_product<Real> product;
      // A_ki of point i, whose recurrence is at degree k, and the recurrence to k + 1. A_ki = p_k(x) exp(log_scale):
      // sqrt(u(x)) underflows far out, where p_k(x) is large.
      auto next = [&](unsigned i, unsigned k) {
         const Real a = f.p[i] * exp(f.log_scale[i]);
         const Real p = ((m_abscissa[i] - m_weight.alpha(k)) * f.p[i] - m_weight.beta(k) * f.p_previous[i]) / m_weight.beta(k + 1);
         f.p_previous[i] = f.p[i];
         f.p[i] = p;
         if (abs(p) > 1 / eps)
         {
            f.p[i] *= eps;
            f.p_previous[i] *= eps;
            f.log_scale[i] -= log(eps);
         }
         return a;
      };
      std::vector<Real> r;
      for (unsigned m = old; m < N; ++m)
      {
         // Point m and degree m: column m and row m of A.
         f.p.push_back(1 / sqrt(mu0));
         f.p_previous.push_back(0);
         f.log_scale.push_back(m_weight.half_log(m_abscissa[m]));
         for (unsigned k = 0; k < m; ++k)
            f.A[k].push_back(next(m, k));
         f.A.emplace_back();
         for (unsigned i = 0; i <= m; ++i)
            f.A[m].push_back(next(i, m));
         // Column m of U, row m of L, and y = L^-1 sqrt(mu0) e_0, whose first n entries are those of every n.
         std::vector<Real> u(m + 1), l(m);
         for (unsigned i = 0; i < m; ++i)
         {
            u[i] = f.A[i][m];
            for (unsigned t = 0; t < i; ++t)
               u[i] -= f.L[i][t] * u[t];
         }
         for (unsigned c = 0; c < m; ++c)
         {
            l[c] = f.A[m][c];
            for (unsigned t = 0; t < c; ++t)
               l[c] -= l[t] * f.U[c][t];
            l[c] /= f.U[c][c];
         }
         u[m] = f.A[m][m];
         Real y = m == 0 ? sqrt(mu0) : Real(0);
         for (unsigned t = 0; t < m; ++t)
         {
            u[m] -= l[t] * u[t];
            y -= l[t] * f.y[t];
         }
         f.L.push_back(std::move(l));
         f.U.push_back(std::move(u));
         f.y.push_back(y);

         const unsigned n = m + 1;
         auto back_substitute = [&](std::vector<Real>& v) {
            for (unsigned i = n; i-- > 0;)
            {
               for (unsigned j = i + 1; j < n; ++j)
                  v[i] -= f.U[j][i] * v[j];
               v[i] /= f.U[i][i];
            }
         };
         std::vector<Real> w(f.y);
         back_substitute(w);
         // One step of iterative refinement, with a compensated residual (Ogita, Rump and Oishi's Dot2): the solve loses
         // far more than the rounding of A, as much as n^2 eps for Laguerre.
         r.resize(n);
         for (unsigned k = 0; k < n; ++k)
         {
            Real sum = k == 0 ? sqrt(mu0) : Real(0);
            Real error = 0;
            for (unsigned i = 0; i < n; ++i)
            {
               Real p, e;
               product(-f.A[k][i], w[i], p, e);
               const Real t = sum + p;
               error += BOOST_MATH_NAMESPACE::detail::two_sum_error(sum, p, t) + e;
               sum = t;
            }
            r[k] = sum + error;
            for (unsigned t = 0; t < k; ++t)
               r[k] -= f.L[k][t] * r[t];
         }
         back_substitute(r);
         for (unsigned i = 0; i < n; ++i)
            w[i] = (w[i] + r[i]) * exp(m_weight.half_log(m_abscissa[i]));
         m_weights.push_back(std::move(w));
      }
   }

   // A_ki, k, i <= m, for the points so far; their recurrences, at degree m + 1; and A = L U, L unit lower triangular,
   // L[m][c] for c < m, U[m][i] = U_im for i <= m, and y = L^-1 sqrt(mu0) e_0. Also the last search: the maximum of F
   // in each gap, F'' there, and the gap of the last point, so that a rule computes the same points however it grows.
   struct factorization
   {
      std::vector<std::vector<Real>> A, L, U;
      std::vector<Real> p, p_previous, log_scale, y, argmax, curvature;
      std::size_t chosen = 0;
   };

   Weight m_weight;
   mutable std::vector<Real> m_abscissa;
   mutable std::vector<std::vector<Real>> m_weights; // m_weights[n] is the n-point rule; m_weights[0] is empty.
   mutable factorization m_factorization;
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
