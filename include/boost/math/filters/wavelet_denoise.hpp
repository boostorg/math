//  Copyright Nick Thompson 2026.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_FILTERS_WAVELET_DENOISE_HPP
#define BOOST_MATH_FILTERS_WAVELET_DENOISE_HPP

#ifdef _MSC_VER
#pragma once
#endif

#include <boost/math/tools/config.hpp>

#ifdef BOOST_MATH_NO_CXX17_IF_CONSTEXPR
#error "The header <boost/math/filters/wavelet_denoise.hpp> requires C++17 or later."
#endif

#ifndef BOOST_MATH_BUILD_MODULE
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iterator>
#include <limits>
#include <memory>
#include <type_traits>
#include <utility>
#include <vector>
#ifdef BOOST_MATH_HAS_THREADS
#include <mutex>
#endif
#endif
#include <boost/math/constants/constants.hpp>
#include <boost/math/filters/daubechies.hpp>
#include <boost/math/filters/detail/wavelet_boundary.hpp>
#include <boost/math/filters/symlets.hpp>
#include <boost/math/filters/trend_filter.hpp>
#include <boost/math/policies/error_handling.hpp>
#include <boost/math/policies/policy.hpp>
#include <boost/math/special_functions/erf.hpp>

// Denoising of equally spaced samples by shrinking their coefficients in an orthogonal wavelet basis on
// the interval; the boundary construction and its rationale are in detail/wavelet_boundary.hpp.

BOOST_MATH_NAMESPACE_BEGIN namespace filters {

enum class threshold_rule { wiener, universal_soft, sure, eb_posterior_median };

namespace detail {

// One level of the transform: n_l samples -> nc = ceil(n_l/2) scaling and nd = floor(n_l/2) wavelet coefficients.
struct dwt_level
{
   std::size_t n, nc, nd;
   bool odd;
};

// The multilevel transform. After forward(), level l's details are at [plan[l].nc, plan[l].n).
template <class T>
class interval_dwt
{
public:
   explicit interval_dwt(std::vector<T> h)
      : h_(std::move(h)), K_(h_.size() / 2), left_(h_, false), right_(std::vector<T>(h_.rbegin(), h_.rend()), true)
   {
      g_.resize(h_.size());
      for (std::size_t i = 0; i < h_.size(); ++i) // g_m = (-1)^m h_{2p-1-m}
         g_[i] = (i % 2 == 0 ? T(1) : T(-1)) * h_[h_.size() - 1 - i];
      if (left_.valid() && right_.valid() && right_.valid_odd())
         width_ = (std::max)((std::max)(left_.width(false), right_.width(false)), right_.width(true));
   }
   interval_dwt(const interval_dwt&) = delete;
   interval_dwt& operator=(const interval_dwt&) = delete;

   bool valid() const { return width_ > 0; }
   std::size_t vanishing_moments() const { return K_; }
   const std::vector<T>& scaling_filter() const { return h_; }

   // The levels for n samples: while ceil(n_l/2) >= 2p and the edge windows fit; empty if n < 4p - 1.
   std::vector<dwt_level> plan(std::size_t n) const
   {
      std::vector<dwt_level> levels;
      if (!valid()) return levels;
      for (std::size_t m = n;;)
      {
         const std::size_t nc = (m + 1) / 2;
         if (nc < 2 * K_ || m < width_) break;
         levels.push_back(dwt_level{m, nc, m - nc, m % 2 == 1});
         m = nc;
      }
      return levels;
   }
   std::size_t max_levels(std::size_t n) const { return plan(n).size(); }

   // for tests: the right-edge row sets held (all-even, plus at most max_odd_cached others)
   static constexpr std::size_t max_odd_cached = 16;
   std::size_t right_cache_size() const
   {
#ifdef BOOST_MATH_HAS_THREADS
      std::lock_guard<std::mutex> lock(mutex_);
#endif
      return (even_right_chain_.rows && !even_right_chain_.rows->empty() ? 1 : 0) + odd_right_rows_.size();
   }

   // The first `levels` levels of plan(x.size()), in place; false if the boundary rows could not be generated.
   bool forward(std::vector<T>& x, std::size_t levels = std::size_t(-1)) const { return transform(x, true, levels); }
   bool inverse(std::vector<T>& x, std::size_t levels = std::size_t(-1)) const { return transform(x, false, levels); }

private:
   using level_type = wavelet_boundary_level<T>;
   using rows_type = std::vector<level_type>;

   bool transform(std::vector<T>& x, bool fwd, std::size_t levels) const
   {
      std::vector<dwt_level> p = plan(x.size());
      if (p.size() > levels) p.resize(levels);
      if (p.empty()) return true;
      std::vector<bool> odd;
      for (const dwt_level& lv : p) odd.push_back(lv.odd);
      std::shared_ptr<const rows_type> L, R;
      rows(odd, L, R);
      if (!L || !R) return false;
      std::vector<T> tmp(x.size());
      for (std::size_t s = 0; s < p.size(); ++s)
      {
         const std::size_t l = fwd ? s : p.size() - 1 - s;
         if (fwd) step((*L)[l], (*R)[l], p[l], x.data(), tmp.data());
         else unstep((*L)[l], (*R)[l], p[l], x.data(), tmp.data());
         std::copy(tmp.begin(), tmp.begin() + std::ptrdiff_t(p[l].n), x.begin());
      }
      return true;
   }

   // The left and all-even right rows grow level by level; other parity sequences live in an LRU.
   void rows(const std::vector<bool>& odd, std::shared_ptr<const rows_type>& L, std::shared_ptr<const rows_type>& R) const
   {
#ifdef BOOST_MATH_HAS_THREADS
      std::lock_guard<std::mutex> lock(mutex_);
#endif
      L = grow(left_, left_chain_, odd.size());
      if (std::find(odd.begin(), odd.end(), true) == odd.end())
      {
         R = grow(right_, even_right_chain_, odd.size());
         return;
      }
      for (auto it = odd_right_rows_.begin(); it != odd_right_rows_.end(); ++it)
         if (it->first == odd)
         {
            std::rotate(it, it + 1, odd_right_rows_.end()); // most recently used last
            R = odd_right_rows_.back().second;
            return;
         }
      rows_type r = right_.rows(odd);
      R = r.size() == odd.size() ? std::make_shared<const rows_type>(std::move(r)) : nullptr;
      if (!R) return;
      if (odd_right_rows_.size() == max_odd_cached) odd_right_rows_.erase(odd_right_rows_.begin());
      odd_right_rows_.emplace_back(odd, R);
   }

   struct chain
   {
      typename boundary_generator<T>::state state;
      std::shared_ptr<const rows_type> rows;
   };
   static std::shared_ptr<const rows_type> grow(const boundary_generator<T>& gen, chain& c, std::size_t levels)
   {
      if (!c.rows) { c.state = gen.start(); c.rows = std::make_shared<const rows_type>(); }
      if (c.rows->size() < levels && c.state.ok)
      {
         rows_type r = *c.rows;
         if (gen.extend(c.state, std::vector<bool>(levels - r.size(), false), r))
            c.rows = std::make_shared<const rows_type>(std::move(r));
      }
      return c.rows->size() >= levels ? c.rows : nullptr;
   }

   // One level: x (n samples) -> out = [coarse (nc) | detail (nd)].
   void step(const level_type& left, const level_type& right, const dwt_level& lv, const T* x, T* out) const
   {
      const std::size_t n = lv.n, nc = lv.nc, nd = lv.nd, N = K_, K = K_;
      T* a = out;
      T* d = out + nc;
      for (std::size_t i = 0; i < K; ++i)
      {
         T sa = 0, sd = 0, ra = 0;
         for (std::size_t c = 0; c < left.width; ++c)
         {
            sa += left.lowpass[i * left.width + c] * x[c];
            sd += left.highpass[i * left.width + c] * x[c];
         }
         for (std::size_t c = 0; c < right.width; ++c) ra += right.lowpass[i * right.width + c] * x[n - 1 - c];
         a[i] = sa; d[i] = sd; a[nc - 1 - i] = ra;
      }
      for (std::size_t i = 0; i < right.nhigh; ++i)
      {
         T rd = 0;
         for (std::size_t c = 0; c < right.width; ++c) rd += right.highpass[i * right.width + c] * x[n - 1 - c];
         d[nd - 1 - i] = rd;
      }
      for (std::size_t k = K; k + K < nc; ++k)
      {
         const T* xk = x + 2 * k - (N - 1);
         T sa = 0, sd = 0;
         for (std::size_t m = 0; m < 2 * N; ++m) { sa += h_[m] * xk[m]; sd += g_[m] * xk[m]; }
         a[k] = sa; d[k] = sd;
      }
   }

   void unstep(const level_type& left, const level_type& right, const dwt_level& lv, const T* in, T* x) const
   {
      const std::size_t n = lv.n, nc = lv.nc, nd = lv.nd, N = K_, K = K_;
      const T* a = in;
      const T* d = in + nc;
      std::fill(x, x + n, T(0));
      for (std::size_t i = 0; i < K; ++i)
      {
         for (std::size_t c = 0; c < left.width; ++c)
            x[c] += left.lowpass[i * left.width + c] * a[i] + left.highpass[i * left.width + c] * d[i];
         for (std::size_t c = 0; c < right.width; ++c) x[n - 1 - c] += right.lowpass[i * right.width + c] * a[nc - 1 - i];
      }
      for (std::size_t i = 0; i < right.nhigh; ++i)
         for (std::size_t c = 0; c < right.width; ++c) x[n - 1 - c] += right.highpass[i * right.width + c] * d[nd - 1 - i];
      for (std::size_t k = K; k + K < nc; ++k)
      {
         T* xk = x + 2 * k - (N - 1);
         for (std::size_t m = 0; m < 2 * N; ++m) xk[m] += h_[m] * a[k] + g_[m] * d[k];
      }
   }

   std::vector<T> h_, g_;
   std::size_t K_;
   boundary_generator<T> left_, right_;
   std::size_t width_ = 0; // set once, by the constructor
   mutable chain left_chain_, even_right_chain_;
   mutable std::vector<std::pair<std::vector<bool>, std::shared_ptr<const rows_type>>> odd_right_rows_;
#ifdef BOOST_MATH_HAS_THREADS
   mutable std::mutex mutex_;
#endif
};

} // namespace detail

namespace detail { struct wavelet_access; }

// A wavelet basis on the interval (db2..db10 or sym2..sym10, from the factories in namespace wavelet).
// Copies share the boundary rows. float is computed in double.
template <class Real>
class orthogonal_wavelet
{
public:
   orthogonal_wavelet() = default;

   // 0 for a default-constructed wavelet
   unsigned vanishing_moments() const { return impl_ ? unsigned(impl_->vanishing_moments()) : 0u; }

private:
   friend struct detail::wavelet_access;
   using compute_type = detail::trend_filter_compute_type<Real>;
   std::shared_ptr<const detail::interval_dwt<compute_type>> impl_;
};

namespace detail {

// The transform behind orthogonal_wavelet. Not part of the public interface yet.
struct wavelet_access
{
   template <class Real>
   static const std::shared_ptr<const interval_dwt<trend_filter_compute_type<Real>>>& impl(const orthogonal_wavelet<Real>& w) { return w.impl_; }
   template <class Real>
   static orthogonal_wavelet<Real> make(std::shared_ptr<const interval_dwt<trend_filter_compute_type<Real>>> impl)
   {
      orthogonal_wavelet<Real> w;
      w.impl_ = std::move(impl);
      return w;
   }
};

// h: an orthonormal scaling filter of even length >= 4 with h.size()/2 vanishing moments
template <class Real, class Container>
orthogonal_wavelet<Real> make_orthogonal_wavelet(const Container& h)
{
   using T = trend_filter_compute_type<Real>;
   std::vector<T> v;
   for (const auto& x : h) v.push_back(static_cast<T>(x));
   if (v.size() < 4 || v.size() % 2 != 0)
      return orthogonal_wavelet<Real>();
   return wavelet_access::make<Real>(std::make_shared<const interval_dwt<T>>(std::move(v)));
}

template <class Real>
const interval_dwt<trend_filter_compute_type<Real>>* wavelet_implementation(const orthogonal_wavelet<Real>& w) { return wavelet_access::impl(w).get(); }

template <class Real>
bool wavelet_valid(const orthogonal_wavelet<Real>& w) { const auto* p = wavelet_implementation(w); return p && p->valid(); }

// The number of levels the transform uses for n samples: the most possible. 0 if n is too short.
template <class Real>
std::size_t wavelet_levels(const orthogonal_wavelet<Real>& w, std::size_t n) { const auto* p = wavelet_implementation(w); return p ? p->max_levels(n) : 0; }

template <class Real>
std::vector<Real> wavelet_scaling_filter(const orthogonal_wavelet<Real>& w)
{
   std::vector<Real> h;
   if (const auto* p = wavelet_implementation(w))
      for (const auto& x : p->scaling_filter()) h.push_back(static_cast<Real>(x));
   return h;
}

} // namespace detail

namespace detail {

// The scaling filter of db2..db10 (index 0..8) or sym2..sym10 (index 9..17).
template <class T>
std::vector<T> wavelet_family_filter(unsigned index)
{
   auto v = [](const auto& a) { return std::vector<T>(a.begin(), a.end()); };
   switch (index)
   {
   case 0: return v(daubechies_scaling_filter<T, 2>());
   case 1: return v(daubechies_scaling_filter<T, 3>());
   case 2: return v(daubechies_scaling_filter<T, 4>());
   case 3: return v(daubechies_scaling_filter<T, 5>());
   case 4: return v(daubechies_scaling_filter<T, 6>());
   case 5: return v(daubechies_scaling_filter<T, 7>());
   case 6: return v(daubechies_scaling_filter<T, 8>());
   case 7: return v(daubechies_scaling_filter<T, 9>());
   case 8: return v(daubechies_scaling_filter<T, 10>());
   case 9: return v(symlet_scaling_filter<T, 2>());
   case 10: return v(symlet_scaling_filter<T, 3>());
   case 11: return v(symlet_scaling_filter<T, 4>());
   case 12: return v(symlet_scaling_filter<T, 5>());
   case 13: return v(symlet_scaling_filter<T, 6>());
   case 14: return v(symlet_scaling_filter<T, 7>());
   case 15: return v(symlet_scaling_filter<T, 8>());
   case 16: return v(symlet_scaling_filter<T, 9>());
   default: return v(symlet_scaling_filter<T, 10>());
   }
}

constexpr unsigned wavelet_family_size = 18;
// The index of db_p (sym = false) or sym_p in the family.
constexpr unsigned wavelet_index(bool sym, unsigned p) { return (sym ? 9u : 0u) + p - 2; }

template <class T, unsigned I>
const orthogonal_wavelet<T>& cached_wavelet()
{
   static const orthogonal_wavelet<T> w = make_orthogonal_wavelet<T>(wavelet_family_filter<T>(I));
   return w;
}

template <class T, std::size_t... I>
const orthogonal_wavelet<T>& cached_wavelet(unsigned index, std::index_sequence<I...>)
{
   using getter = const orthogonal_wavelet<T>& (*)();
   static constexpr getter table[] = {&cached_wavelet<T, unsigned(I)>...};
   return table[index]();
}

// The wavelet of the family (cached, built on first use), as orthogonal_wavelet<Real> sharing its rows.
template <class T>
const orthogonal_wavelet<T>& cached_wavelet(unsigned index)
{
   return cached_wavelet<T>(index, std::make_index_sequence<wavelet_family_size>());
}

template <class Real>
orthogonal_wavelet<Real> make_wavelet(unsigned p, unsigned offset, const char* function)
{
   if (p < 2 || p > 10)
   {
      policies::raise_domain_error(function, "The number of vanishing moments must be in [2, 10], but got %1%.", Real(p), policies::policy<>());
      return orthogonal_wavelet<Real>();
   }
   return wavelet_access::make<Real>(wavelet_access::impl(cached_wavelet<trend_filter_compute_type<Real>>(offset + p - 2)));
}

} // namespace detail

namespace wavelet {

// Daubechies' minimum-phase filter with p = 2..10 vanishing moments.
template <class Real>
orthogonal_wavelet<Real> daubechies(unsigned p)
{
   return detail::make_wavelet<Real>(p, 0, "boost::math::filters::wavelet::daubechies<%1%>");
}

// The least-asymmetric filter with p = 2..10 vanishing moments (see symlets.hpp for the convention).
template <class Real>
orthogonal_wavelet<Real> symlet(unsigned p)
{
   return detail::make_wavelet<Real>(p, 9, "boost::math::filters::wavelet::symlet<%1%>");
}

} // namespace wavelet

namespace detail {

// Empirical-Bayes posterior median (Johnstone & Silverman 2004, 2005; R package EbayesThresh): prior
// w Laplace(a) + (1 - w) delta_0 on the standardized coefficient, a = 1/2 (EbayesThresh's default),
// w by marginal maximum likelihood, bounded below by the weight whose threshold is the universal one.
// The cutoffs at 35 and 25 are EbayesThresh's.
namespace eb {

template <class T>
T erfc(const T& x)
{
   if constexpr (std::is_floating_point<T>::value) return std::erfc(x);
   else return BOOST_MATH_NAMESPACE::erfc(x);
}
template <class T>
T Phi(const T& z) { return erfc(-z * constants::half_root_two<T>()) / 2; }
template <class T>
T phi(const T& z) { using std::exp; return exp(-z * z / 2) * constants::one_div_root_two_pi<T>(); }

// g(x) / phi(x) - 1, g the Laplace(a) density convolved with phi
template <class T>
T beta_laplace(T x, const T& a)
{
   using std::abs;
   x = abs(x);
   const T xpa = x + a, xma = (std::min)(T(x - a), T(35));
   const T r1 = xpa < 35 ? T(Phi(T(-xpa)) / phi(xpa)) : T(1 / xpa);
   const T r2 = xma > -35 ? T(Phi(xma) / phi(xma)) : T(1 / abs(xma));
   return a / 2 * (r1 + r2) - 1;
}

// the weight whose posterior median threshold is t
template <class T>
T weight_from_threshold(const T& t, const T& a)
{
   using std::abs;
   const T tma = t - a;
   T wi = tma > -35 ? T(Phi(tma) / phi(tma)) : T(1 / abs(tma));
   wi = a * wi - beta_laplace(t, a);
   return 1 / wi;
}

// marginal maximum likelihood weight for the m standardized coefficients x, by bisection in log w
template <class T>
T weight(const T* x, std::size_t m, const T& a)
{
   using std::sqrt; using std::log;
   T wlo = weight_from_threshold(T(sqrt(2 * log(T(m)))), a), whi = 1;
   std::vector<T> b(m);
   for (std::size_t i = 0; i < m; ++i) b[i] = (std::min)(beta_laplace(x[i], a), T(1e20));
   auto score = [&](const T& w) { T s = 0; for (const T& v : b) s += v / (1 + w * v); return s; };
   if (score(T(1)) >= 0) return 1;
   if (score(wlo) <= 0) return wlo;
   for (int j = 0; j < 30; ++j)
   {
      const T wm = sqrt(wlo * whi), sm = score(wm);
      if (sm == 0) return wm;
      (sm > 0 ? wlo : whi) = wm;
   }
   return sqrt(wlo * whi);
}

template <class T>
T posterior_median(const T& x, const T& w, const T& a)
{
   using std::abs; using std::sqrt;
   const T sx = x < 0 ? T(-1) : T(1), ax = abs(x), xma = ax - a;
   T zz = 1 / a * phi(xma) * (1 / w + beta_laplace(ax, a));
   if (xma > 25) zz = T(1) / 2;
   if (zz >= 1) return T(0);
   if (zz <= 0) return sx * xma;
   const T mucor = -sqrt(T(2)) * BOOST_MATH_NAMESPACE::erfc_inv(2 * zz); // Phi^{-1}(zz)
   return sx * (std::max)(T(0), T(xma - mucor));
}

} // namespace eb

template <class T>
T soft_threshold(const T& x, const T& t) { using std::abs; return abs(x) <= t ? T(0) : (x > 0 ? T(x - t) : T(x + t)); }

// sigma from the finest details: their median absolute value over Phi^{-1}(3/4)
template <class T>
T mad_sigma(const T* d, std::size_t m)
{
   using std::abs; using std::sqrt;
   std::vector<T> a(d, d + m);
   for (T& v : a) v = abs(v);
   auto mid = a.begin() + std::ptrdiff_t(m / 2);
   std::nth_element(a.begin(), mid, a.end());
   T med = *mid;
   if (m % 2 == 0) med = (med + *std::max_element(a.begin(), mid)) / 2;
   return med / (sqrt(T(2)) * BOOST_MATH_NAMESPACE::erfc_inv(T(1) / 2));
}

// SureShrink threshold of one level (Donoho & Johnstone 1995) in units of sigma, with the hybrid rule:
// sparse levels get the level's universal threshold sqrt(2 ln m).
template <class T>
T sure_threshold(const T* d, std::size_t m, const T& sigma)
{
   using std::sqrt; using std::log; using std::pow;
   const T universal = sqrt(2 * log(T(m)));
   std::vector<T> a(m);
   T energy = 0;
   for (std::size_t i = 0; i < m; ++i) { a[i] = d[i] / sigma; a[i] *= a[i]; energy += a[i] - 1; }
   if (energy / T(m) <= pow(T(log(T(m)) / constants::ln_two<T>()), T(3) / 2) / sqrt(T(m)))
      return universal;
   std::sort(a.begin(), a.end());
   // SURE(0) = m; SURE(t = sqrt(a_i)) = m - 2 (i+1) + sum_{j<=i} a_j + (m-1-i) a_i
   T best = T(m), t2 = 0, cum = 0;
   for (std::size_t i = 0; i < m; ++i)
   {
      cum += a[i];
      const T risk = T(m) - 2 * T(i + 1) + cum + T(m - 1 - i) * a[i];
      if (risk < best) { best = risk; t2 = a[i]; }
   }
   return (std::min)(T(sqrt(t2)), universal);
}

// Shrinks the details of c (laid out as in plan) by rule; pilot: the coefficients of the Wiener pilot estimate.
// wiener (Ghael, Sayeed & Baraniuk 1997): weights theta^2 / (theta^2 + sigma^2), theta the pilot's coefficient.
template <class T>
void shrink(std::vector<T>& c, const std::vector<dwt_level>& plan, threshold_rule rule, const T& sigma, const std::vector<T>& pilot)
{
   using std::sqrt; using std::log; using std::abs;
   if (!(sigma > 0)) return;
   const T universal = sigma * sqrt(2 * log(T(c.size())));
   for (const dwt_level& lv : plan)
   {
      T* d = c.data() + lv.nc;
      const std::size_t len = lv.nd;
      T t = -1; // the soft threshold, if any
      switch (rule)
      {
      case threshold_rule::wiener:
         for (std::size_t i = 0; i < len; ++i)
         {
            const T& th = pilot[lv.nc + i];
            d[i] *= th * th / (th * th + sigma * sigma);
         }
         break;
      case threshold_rule::universal_soft:
         t = universal;
         break;
      case threshold_rule::sure:
         t = sigma * sure_threshold(d, len, sigma);
         break;
      case threshold_rule::eb_posterior_median:
      {
         const T a = T(1) / 2;
         std::vector<T> x(len);
         for (std::size_t i = 0; i < len; ++i) x[i] = d[i] / sigma;
         const T w = eb::weight(x.data(), len, a);
         for (std::size_t i = 0; i < len; ++i) d[i] = sigma * eb::posterior_median(x[i], w, a);
         break;
      }
      }
      if (t >= 0)
         for (std::size_t i = 0; i < len; ++i) d[i] = soft_threshold(d[i], t);
   }
}

enum class wavelet_status { ok, too_short, cv_too_short, invalid_wavelet, unknown_rule, bad_sigma, generator_failed, too_many_levels, non_finite };

// The Wiener pilot for W at n samples: the first of sym8, db4, sym4, db2 that differs from W and takes n samples
// (an index into the family), else W itself (wavelet_family_size).
constexpr unsigned wiener_pilots[] = {wavelet_index(true, 8), wavelet_index(false, 4), wavelet_index(true, 4), wavelet_index(false, 2)};

template <class T>
unsigned wiener_pilot(const interval_dwt<T>& W, std::size_t n)
{
   for (unsigned q : wiener_pilots)
   {
      const interval_dwt<T>& P = *wavelet_implementation(cached_wavelet<T>(q));
      if (P.scaling_filter() != W.scaling_filter() && P.max_levels(n) > 0) return q;
   }
   return wavelet_family_size;
}

// y replaced by its universal-hard estimate in P (sigma NaN: the MAD estimate); false if generation fails.
template <class T>
bool pilot_estimate(const interval_dwt<T>& P, std::vector<T>& y, T sigma)
{
   using std::isnan; using std::sqrt; using std::log; using std::abs;
   const std::vector<dwt_level> plan = P.plan(y.size());
   if (!P.forward(y)) return false;
   if (isnan(sigma)) sigma = mad_sigma(y.data() + plan[0].nc, plan[0].nd);
   const T universal = sigma * sqrt(2 * log(T(y.size())));
   for (const dwt_level& lv : plan)
      for (std::size_t i = lv.nc; i < lv.nc + lv.nd; ++i)
         if (abs(y[i]) <= universal) y[i] = 0;
   return P.inverse(y);
}

// y denoised in W by rule; sigma NaN: the MAD estimate from the finest details.
template <class T>
wavelet_status wavelet_denoise_impl(const interval_dwt<T>& W, std::vector<T>& y, threshold_rule rule, T sigma)
{
   using std::isnan;
   const std::vector<dwt_level> plan = W.plan(y.size());
   if (plan.empty()) return wavelet_status::too_short;
   std::vector<T> pilot;
   if (rule == threshold_rule::wiener)
   {
      pilot = y;
      const unsigned q = wiener_pilot(W, y.size());
      if (!pilot_estimate(q < wavelet_family_size ? *wavelet_implementation(cached_wavelet<T>(q)) : W, pilot, sigma) || !W.forward(pilot))
         return wavelet_status::generator_failed;
   }
   if (!W.forward(y)) return wavelet_status::generator_failed;
   if (isnan(sigma)) sigma = mad_sigma(y.data() + plan[0].nc, plan[0].nd);
   shrink(y, plan, rule, sigma, pilot);
   return W.inverse(y) ? wavelet_status::ok : wavelet_status::generator_failed;
}

// The first sample that is NaN or infinite, else nullptr: thresholds sort and select, which need finite values.
template <class T>
const T* non_finite(const std::vector<T>& y)
{
   using std::isfinite;
   for (const T& v : y)
      if (!(isfinite)(v)) return &v;
   return nullptr;
}

// v filled with the error handler's result for s; value is the message's %1%.
template <class RandomAccessContainer>
RandomAccessContainer wavelet_error(const RandomAccessContainer& v, const char* function, wavelet_status s, const typename RandomAccessContainer::value_type& value)
{
   static const char* const messages[] = {
      "",
      "%1% samples are too few for this wavelet (the minimum is 4p - 1).",
      "%1% samples are too few: cross-validation needs at least 14 (7 per half, the db2 minimum).",
      "The wavelet is invalid (default-constructed, or its boundary rows could not be set up; %1% vanishing moments).",
      "Unknown threshold rule %1%.",
      "sigma must be finite and >= 0, or NaN, but got %1%.",
      "The boundary rows could not be generated for %1% samples.",
      "%1% levels requested, more than the length allows for this wavelet.",
      "The samples must be finite, but one is %1%."};
   RandomAccessContainer out = v;
   std::fill(std::begin(out), std::end(out), s == wavelet_status::generator_failed
      ? policies::raise_evaluation_error(function, messages[int(s)], value, policies::policy<>())
      : policies::raise_domain_error(function, messages[int(s)], value, policies::policy<>()));
   return out;
}

// The rules the cross-validated wavelet_denoise chooses among, in order of preference on ties.
constexpr threshold_rule cv_rules[] = {threshold_rule::sure, threshold_rule::universal_soft, threshold_rule::eb_posterior_median, threshold_rule::wiener};
static_assert(std::size(cv_rules) == std::size_t(threshold_rule::eb_posterior_median) + 1, "cv_rules must list every rule");

template <class T>
struct wavelet_choice
{
   T score = (std::numeric_limits<T>::max)();
   unsigned wavelet = 0;
   threshold_rule rule = threshold_rule::sure;
};

template <class T>
void split_even_odd(const std::vector<T>& y, std::vector<T>& ye, std::vector<T>& yo)
{
   ye.resize((y.size() + 1) / 2);
   yo.resize(y.size() / 2);
   for (std::size_t i = 0; i < y.size(); ++i) (i % 2 ? yo[i / 2] : ye[i / 2]) = y[i];
}

// Even/odd cross-validation (Nason 1996) over db2..db10, sym4..sym10 and cv_rules: each half is
// denoised, the other half predicted by averaging neighbouring values (cross_validation_score).
// Wavelets that do not take both halves are skipped, and so is a Wiener candidate whose pilot fails.
template <class T>
wavelet_choice<T> wavelet_select(const std::vector<T>& y, wavelet_status& status)
{
   std::vector<T> half[2], fit[2][std::size(cv_rules)];
   split_even_odd(y, half[0], half[1]);
   // the pilot estimates of both halves, per pilot wavelet (index wavelet_family_size: W itself), computed once
   std::vector<T> estimates[wavelet_family_size + 1][2];
   int estimated[wavelet_family_size + 1] = {}; // 0: not yet, 1: done, -1: failed
   auto pilot_of = [&](unsigned q, const interval_dwt<T>& P) {
      if (q == wavelet_family_size || estimated[q] == 0)
      {
         bool ok = true;
         for (int s = 0; s < 2; ++s)
         {
            estimates[q][s] = half[s];
            ok = ok && pilot_estimate(P, estimates[q][s], std::numeric_limits<T>::quiet_NaN());
         }
         estimated[q] = ok ? 1 : -1;
      }
      return estimated[q] == 1;
   };
   wavelet_choice<T> best;
   status = wavelet_status::cv_too_short;
   for (unsigned t = 0; t < wavelet_family_size; ++t)
   {
      if (t == wavelet_index(true, 2) || t == wavelet_index(true, 3)) continue; // sym2, sym3 are db2, db3
      const interval_dwt<T>& W = *wavelet_implementation(cached_wavelet<T>(t));
      if (W.max_levels(half[0].size()) == 0 || W.max_levels(half[1].size()) == 0) continue;
      const unsigned q = wiener_pilot(W, half[1].size());
      bool wiener_ok = pilot_of(q, q < wavelet_family_size ? *wavelet_implementation(cached_wavelet<T>(q)) : W);
      for (int s = 0; s < 2; ++s)
      {
         std::vector<T> c = half[s], pilot = wiener_ok ? estimates[q][s] : std::vector<T>();
         const std::vector<dwt_level> plan = W.plan(c.size());
         if (!W.forward(c)) { status = wavelet_status::generator_failed; return best; }
         wiener_ok = wiener_ok && W.forward(pilot);
         const T sigma = mad_sigma(c.data() + plan[0].nc, plan[0].nd);
         for (std::size_t k = 0; k < std::size(cv_rules); ++k)
         {
            fit[s][k] = c;
            if (cv_rules[k] == threshold_rule::wiener && !wiener_ok) continue;
            shrink(fit[s][k], plan, cv_rules[k], sigma, pilot);
            if (!W.inverse(fit[s][k])) { status = wavelet_status::generator_failed; return best; }
         }
      }
      for (std::size_t k = 0; k < std::size(cv_rules); ++k)
      {
         if (cv_rules[k] == threshold_rule::wiener && !wiener_ok) continue;
         const T score = cross_validation_score(fit[0][k], fit[1][k], half[0], half[1]);
         if (score < best.score || status != wavelet_status::ok) { best.score = score; best.wavelet = t; best.rule = cv_rules[k]; status = wavelet_status::ok; }
      }
   }
   return best;
}

} // namespace detail

// Denoises samples of a function on an interval by thresholding the coefficients of x in the basis w.
// x.size() must be at least 4p - 1. sigma: the noise standard deviation; NaN: the median
// absolute deviation of the finest details over 0.6745.
template <class RandomAccessContainer>
RandomAccessContainer wavelet_denoise(const RandomAccessContainer& x, const orthogonal_wavelet<typename RandomAccessContainer::value_type>& w,
                                      threshold_rule rule = threshold_rule::wiener,
                                      typename RandomAccessContainer::value_type sigma = std::numeric_limits<typename RandomAccessContainer::value_type>::quiet_NaN())
{
   using Real = typename RandomAccessContainer::value_type;
   using T = detail::trend_filter_compute_type<Real>;
   using std::isnan;
   using detail::wavelet_status;
   static const char* function = "boost::math::filters::wavelet_denoise<%1%>";
   std::vector<T> y(std::begin(x), std::end(x));
   wavelet_status s = wavelet_status::ok;
   if (!detail::wavelet_valid(w))
      return detail::wavelet_error(x, function, wavelet_status::invalid_wavelet, Real(w.vanishing_moments()));
   if (int(rule) < 0 || int(rule) > int(threshold_rule::eb_posterior_median))
      return detail::wavelet_error(x, function, wavelet_status::unknown_rule, Real(int(rule)));
   if (!isnan(sigma) && !(sigma >= 0 && sigma <= (std::numeric_limits<Real>::max)()))
      return detail::wavelet_error(x, function, wavelet_status::bad_sigma, sigma);
   if (const T* bad = detail::non_finite(y))
      return detail::wavelet_error(x, function, wavelet_status::non_finite, Real(*bad));
   s = detail::wavelet_denoise_impl(*detail::wavelet_implementation(w), y, rule, T(sigma));
   return s == wavelet_status::ok ? detail::trend_filter_output(x, y) : detail::wavelet_error(x, function, s, Real(y.size()));
}

// The wavelet (db2..db10, sym2..sym10) and rule chosen by even/odd cross-validation; x.size() >= 14.
template <class RandomAccessContainer>
RandomAccessContainer wavelet_denoise(const RandomAccessContainer& x)
{
   using Real = typename RandomAccessContainer::value_type;
   using T = detail::trend_filter_compute_type<Real>;
   static const char* function = "boost::math::filters::wavelet_denoise<%1%>";
   std::vector<T> y(std::begin(x), std::end(x));
   detail::wavelet_status s;
   if (const T* bad = detail::non_finite(y))
      return detail::wavelet_error(x, function, detail::wavelet_status::non_finite, Real(*bad));
   const detail::wavelet_choice<T> c = detail::wavelet_select(y, s);
   if (s == detail::wavelet_status::ok)
      s = detail::wavelet_denoise_impl(*detail::wavelet_implementation(detail::cached_wavelet<T>(c.wavelet)), y, c.rule, std::numeric_limits<T>::quiet_NaN());
   return s == detail::wavelet_status::ok ? detail::trend_filter_output(x, y) : detail::wavelet_error(x, function, s, Real(y.size()));
}

namespace detail {

template <class RandomAccessContainer>
RandomAccessContainer wavelet_transform(const RandomAccessContainer& x, const orthogonal_wavelet<typename RandomAccessContainer::value_type>& w,
                                        std::size_t levels, bool forward, const char* function)
{
   using Real = typename RandomAccessContainer::value_type;
   using T = trend_filter_compute_type<Real>;
   std::vector<T> c(std::begin(x), std::end(x));
   const std::size_t L = wavelet_levels(w, c.size());
   if (!wavelet_valid(w)) return wavelet_error(x, function, wavelet_status::invalid_wavelet, Real(w.vanishing_moments()));
   if (L == 0) return wavelet_error(x, function, wavelet_status::too_short, Real(c.size()));
   if (levels != std::size_t(-1) && levels > L) return wavelet_error(x, function, wavelet_status::too_many_levels, Real(levels));
   const interval_dwt<T>& W = *wavelet_implementation(w);
   if (!(forward ? W.forward(c, levels) : W.inverse(c, levels)))
      return wavelet_error(x, function, wavelet_status::generator_failed, Real(c.size()));
   return trend_filter_output(x, c);
}

} // namespace detail

// The coefficients of x in the basis w, in one vector: [a_L | d_L | d_{L-1} | ... | d_1], d_1 the finest details.
// Level l takes n_{l-1} values to ceil(n_{l-1}/2) scaling and floor(n_{l-1}/2) wavelet coefficients (n_0 = n).
// levels: the number of levels L, at most (and by default) as many as n allows; x.size() must be at least 4p - 1.
template <class RandomAccessContainer>
RandomAccessContainer wavelet_transform(const RandomAccessContainer& x, const orthogonal_wavelet<typename RandomAccessContainer::value_type>& w,
                                        std::size_t levels = std::size_t(-1))
{
   return detail::wavelet_transform(x, w, levels, true, "boost::math::filters::wavelet_transform<%1%>");
}

// The inverse of wavelet_transform(x, w, levels): x from its coefficients c.
template <class RandomAccessContainer>
RandomAccessContainer inverse_wavelet_transform(const RandomAccessContainer& c, const orthogonal_wavelet<typename RandomAccessContainer::value_type>& w,
                                                std::size_t levels = std::size_t(-1))
{
   return detail::wavelet_transform(c, w, levels, false, "boost::math::filters::inverse_wavelet_transform<%1%>");
}

// Whichever of wavelet_denoise(x) and trend_filter(x) has the lower even/odd cross-validation score.
// If the trend filter fails, the wavelet result is returned.
template <class RandomAccessContainer>
RandomAccessContainer denoise(const RandomAccessContainer& x)
{
   using Real = typename RandomAccessContainer::value_type;
   using T = detail::trend_filter_compute_type<Real>;
   static const char* function = "boost::math::filters::denoise<%1%>";
   std::vector<T> y(std::begin(x), std::end(x));
   detail::wavelet_status s;
   if (const T* bad = detail::non_finite(y))
      return detail::wavelet_error(x, function, detail::wavelet_status::non_finite, Real(*bad));
   const detail::wavelet_choice<T> wc = detail::wavelet_select(y, s);
   if (s != detail::wavelet_status::ok)
      return detail::wavelet_error(x, function, s, Real(y.size()));
   // a trend filter that fails is dropped (ignore_error fills its fit with lambda, so the failure is reported by tf_ok)
   using quiet = policies::policy<policies::evaluation_error<policies::ignore_error>>;
   const T tf_sigma = detail::trend_filter_sigma(y);
   std::vector<T> tf = y;
   T tf_score;
   bool tf_ok = true;
   if (tf_sigma == 0)
   {  // trend_filter returns the data: score the halves predicting each other
      std::vector<T> ye, yo;
      detail::split_even_odd(y, ye, yo);
      tf_score = detail::cross_validation_score(ye, yo, ye, yo);
   }
   else
   {
      const std::pair<int, T> kl = detail::trend_filter_select(y, tf_sigma, &tf_score, quiet());
      if (!(tf_score < wc.score)) tf_score = (std::numeric_limits<T>::max)(); // not needed
      else tf = detail::trend_filter_fit(y, kl.first, kl.second, quiet(), &tf_ok);
   }
   if (wc.score <= tf_score || !tf_ok)
      s = detail::wavelet_denoise_impl(*detail::wavelet_implementation(detail::cached_wavelet<T>(wc.wavelet)), y, wc.rule, std::numeric_limits<T>::quiet_NaN());
   else
      y = tf;
   return s == detail::wavelet_status::ok ? detail::trend_filter_output(x, y) : detail::wavelet_error(x, function, s, Real(y.size()));
}

} BOOST_MATH_NAMESPACE_END

#endif
