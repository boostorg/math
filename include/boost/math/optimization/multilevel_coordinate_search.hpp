/*
 * Copyright Nick Thompson, 2026
 * Use, modification and distribution are subject to the
 * Boost Software License, Version 1.0. (See accompanying file
 * LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef BOOST_MATH_OPTIMIZATION_MULTILEVEL_COORDINATE_SEARCH_HPP
#define BOOST_MATH_OPTIMIZATION_MULTILEVEL_COORDINATE_SEARCH_HPP

#include <boost/math/tools/config.hpp>
#ifndef BOOST_MATH_BUILD_MODULE
#include <atomic>
#endif
#include <boost/math/optimization/detail/common.hpp>
#include <boost/math/tools/assert.hpp>
#ifndef BOOST_MATH_BUILD_MODULE
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>
#endif

BOOST_MATH_NAMESPACE_BEGIN namespace optimization {

// Huyer, W., Neumaier, A. (1999). Global optimization by multilevel coordinate search.
// Journal of Global Optimization, 14, 331-355.
// See: https://arnold-neumaier.at/ms/mcs.pdf
// Section, step and equation numbers in the comments below refer to this paper.

BOOST_MATH_EXPORT template <typename ArgumentContainer> struct multilevel_coordinate_search_parameters {
  using Real = typename ArgumentContainer::value_type;
  ArgumentContainer lower_bounds;
  ArgumentContainer upper_bounds;
  // Used as the base point x^0 of the root box; the default is the midpoint of the box.
  ArgumentContainer const *initial_guess = nullptr;
  // 0 means 1000 times the dimension.
  size_t max_function_calls = 0;
  // The number of levels s_max; 0 means 5n + 10, the value used in Section 7.3.
  size_t max_level = 0;
  // Start local searches from the base points of the boxes that reach level s_max (Section 5).
  bool local_search = true;
  // The limit on visits to Step 3 per local search (Section 7.3).
  size_t max_local_search_steps = 50;
  // smaxls in Section 7.3: the most function values in one line search of the local search.
  size_t line_search_points = 15;
};

BOOST_MATH_EXPORT template <typename ArgumentContainer>
void validate_multilevel_coordinate_search_parameters(multilevel_coordinate_search_parameters<ArgumentContainer> const &params) {
  std::ostringstream oss;
  detail::validate_bounds(params.lower_bounds, params.upper_bounds);
  for (size_t i = 0; i < params.lower_bounds.size(); ++i) {
    if (!(params.lower_bounds[i] < params.upper_bounds[i])) {
      oss << __FILE__ << ":" << __LINE__ << ":" << __func__;
      oss << ": MCS requires lower_bounds[i] < upper_bounds[i], but at i = " << i << " the bounds are "
          << params.lower_bounds[i] << " and " << params.upper_bounds[i] << ".";
      throw std::domain_error(oss.str());
    }
  }
  if (params.initial_guess) {
    detail::validate_initial_guess(*params.initial_guess, params.lower_bounds, params.upper_bounds);
  }
  if (params.max_level == 1) {
    oss << __FILE__ << ":" << __LINE__ << ":" << __func__;
    oss << ": There must be at least two levels.";
    throw std::invalid_argument(oss.str());
  }
  if (params.line_search_points < 3) {
    oss << __FILE__ << ":" << __LINE__ << ":" << __func__;
    oss << ": A line search needs at least 3 points, but line_search_points = " << params.line_search_points << ".";
    throw std::invalid_argument(oss.str());
  }
}

namespace detail {

// Minimizes g^T h + h^T G h/2 over lo <= h <= hi, where lo <= 0 <= hi and G is symmetric (row-major n x n) but may be
// indefinite. This stands in for MINQ: sweeps of exact coordinate minimization, each followed by a Newton step on the
// free variables when the Hessian restricted to them is positive definite. It finds a local minimizer, as MINQ does.
template <typename Real>
std::vector<Real> minimize_quadratic_on_box(std::vector<Real> const &g, std::vector<Real> const &G,
                                            std::vector<Real> const &lo, std::vector<Real> const &hi) {
  using std::abs;
  using std::isfinite;
  using std::sqrt;
  size_t const n = g.size();
  std::vector<Real> h(n, Real(0));
  for (size_t i = 0; i < n; ++i) {
    if (!isfinite(g[i])) {
      return h;
    }
  }
  for (auto const &Gij : G) {
    if (!isfinite(Gij)) {
      return h;
    }
  }
  // r = g + G h, the gradient at h:
  std::vector<Real> r = g;
  auto add_column = [&](size_t j, Real t) {
    for (size_t i = 0; i < n; ++i) {
      r[i] += t * G[i * n + j];
    }
  };
  // Moves h_j by t within [lo_j, hi_j], keeping r = g + G h:
  auto move = [&](size_t j, Real t) {
    Real const old = h[j];
    h[j] = (std::min)((std::max)(h[j] + t, lo[j]), hi[j]);
    add_column(j, h[j] - old);
  };
  auto value = [&]() {
    Real q = 0;
    for (size_t i = 0; i < n; ++i) {
      q += h[i] * (g[i] + r[i]) / 2;
    }
    return q;
  };
  std::vector<size_t> free_set;
  std::vector<Real> L, p;
  Real q_old = 0;
  for (size_t sweep = 0; sweep < 10 * n + 10; ++sweep) {
    for (size_t i = 0; i < n; ++i) {
      // Minimize r_i t + G_ii t^2/2 over h_i + t in [lo_i, hi_i]:
      Real const a = G[i * n + i];
      Real const tlo = lo[i] - h[i];
      Real const thi = hi[i] - h[i];
      Real t = 0;
      if (a > 0) {
        t = (std::min)((std::max)(-r[i] / a, tlo), thi);
      }
      else {
        // Not convex along coordinate i, so the minimum is at an end of the interval (or t = 0 if neither is lower):
        Real const qlo = r[i] * tlo + a * tlo * tlo / 2;
        Real const qhi = r[i] * thi + a * thi * thi / 2;
        if (qlo < 0 && qlo <= qhi) {
          t = tlo;
        }
        else if (qhi < 0) {
          t = thi;
        }
      }
      if (t != 0) {
        move(i, t);
      }
    }
    // Newton step on the free variables, by a hand-rolled Cholesky factorization:
    free_set.clear();
    for (size_t i = 0; i < n; ++i) {
      if (lo[i] < h[i] && h[i] < hi[i]) {
        free_set.push_back(i);
      }
    }
    size_t const m = free_set.size();
    bool positive_definite = m > 0;
    L.assign(m * m, Real(0));
    for (size_t j = 0; j < m && positive_definite; ++j) {
      Real d = G[free_set[j] * n + free_set[j]];
      for (size_t k = 0; k < j; ++k) {
        d -= L[j * m + k] * L[j * m + k];
      }
      if (!(d > 0)) {
        positive_definite = false;
        // w = (-G_11^{-1} G_1j, 1, 0, ...) has curvature w^T G w = d <= 0, where G_11 = L_11 L_11^T is the leading
        // block already factored. Go downhill along it as far as the box allows.
        p.assign(m, Real(0));
        p[j] = 1;
        for (size_t k = j; k-- > 0;) {
          Real s = L[j * m + k];
          for (size_t q = k + 1; q < j; ++q) {
            s += L[q * m + k] * p[q];
          }
          p[k] = -s / L[k * m + k];
        }
        Real slope = 0;
        for (size_t i = 0; i <= j; ++i) {
          slope += r[free_set[i]] * p[i];
        }
        if (slope > 0) {
          for (auto &pi : p) {
            pi = -pi;
          }
        }
        if (d < 0 || slope != 0) {
          Real t = std::numeric_limits<Real>::infinity();
          for (size_t i = 0; i <= j; ++i) {
            size_t const k = free_set[i];
            if (p[i] > 0) {
              t = (std::min)(t, (hi[k] - h[k]) / p[i]);
            }
            else if (p[i] < 0) {
              t = (std::min)(t, (lo[k] - h[k]) / p[i]);
            }
          }
          if (t > 0 && t < std::numeric_limits<Real>::infinity()) {
            for (size_t i = 0; i <= j; ++i) {
              size_t const k = free_set[i];
              move(k, t * p[i]);
            }
          }
        }
        break;
      }
      L[j * m + j] = sqrt(d);
      for (size_t i = j + 1; i < m; ++i) {
        Real s = G[free_set[i] * n + free_set[j]];
        for (size_t k = 0; k < j; ++k) {
          s -= L[i * m + k] * L[j * m + k];
        }
        L[i * m + j] = s / L[j * m + j];
      }
    }
    if (positive_definite) {
      // Solve L L^T p = -r_F:
      p.assign(m, Real(0));
      for (size_t i = 0; i < m; ++i) {
        Real s = -r[free_set[i]];
        for (size_t k = 0; k < i; ++k) {
          s -= L[i * m + k] * p[k];
        }
        p[i] = s / L[i * m + i];
      }
      for (size_t i = m; i-- > 0;) {
        Real s = p[i];
        for (size_t k = i + 1; k < m; ++k) {
          s -= L[k * m + i] * p[k];
        }
        p[i] = s / L[i * m + i];
      }
      // The largest step t <= 1 along p that stays in the box:
      Real t = 1;
      for (size_t i = 0; i < m; ++i) {
        size_t const j = free_set[i];
        if (p[i] > 0) {
          t = (std::min)(t, (hi[j] - h[j]) / p[i]);
        }
        else if (p[i] < 0) {
          t = (std::min)(t, (lo[j] - h[j]) / p[i]);
        }
      }
      if (t > 0) {
        for (size_t i = 0; i < m; ++i) {
          size_t const j = free_set[i];
          move(j, t * p[i]);
        }
      }
    }
    Real const q = value();
    if (!(q < q_old - 4 * std::numeric_limits<Real>::epsilon() * abs(q))) {
      break;
    }
    q_old = q;
  }
  return h;
}

// The state of an MCS run. multilevel_coordinate_search() below drives it. It is separate so that the boxes can be inspected between sweeps.
template <typename ArgumentContainer, class Func> class multilevel_coordinate_search_engine {
public:
  using Real = typename ArgumentContainer::value_type;
  using ResultType = std::invoke_result_t<Func, ArgumentContainer>;

  // A box as stored: instead of its bounds, its history (Section 4), so that the storage is proportional to the number
  // of function calls. That is its parent, the coordinate i it was split along, and the ith coordinates of its base
  // point x and opposite point y and of its xi1, df1, xi2, df2 below. The root is box 0.
  struct box {
    size_t parent;
    size_t i;
    Real x;
    Real y;
    Real f;
    // 0 for a box that has been split:
    size_t level;
    Real xi1, df1, xi2, df2;
  };

  // A box B[x, y] (Section 3) with base point x and opposite point y, as load() recovers it from its history.
  struct explicit_box {
    std::vector<Real> x;
    std::vector<Real> y;
    Real f;
    // 0 for a box that has been split:
    size_t level;
    // The number of times each coordinate has been split in the history of the box:
    std::vector<size_t> splits;
    // For each coordinate, two more points (the closest in the history of the box) and the differences of their
    // function values from that of the base point they were split from, for the separable model of Section 4.2:
    std::vector<Real> xi1, df1, xi2, df2;
  };

  multilevel_coordinate_search_engine(Func const &cost_function, multilevel_coordinate_search_parameters<ArgumentContainer> const &params, ResultType target_value,
             std::atomic<bool> *cancellation, std::atomic<ResultType> *current_minimum_cost,
             std::vector<std::pair<ArgumentContainer, ResultType>> *queries)
      : cost_function_(cost_function), params_(params), target_value_(target_value), cancellation_(cancellation),
        current_minimum_cost_(current_minimum_cost), queries_(queries) {
    validate_multilevel_coordinate_search_parameters(params);
    n_ = params.lower_bounds.size();
    u_.assign(params.lower_bounds.begin(), params.lower_bounds.end());
    v_.assign(params.upper_bounds.begin(), params.upper_bounds.end());
    smax_ = params.max_level ? params.max_level : 5 * n_ + 10;
    max_function_calls_ = params.max_function_calls ? params.max_function_calls : 1000 * n_;
    heaps_.resize(smax_);
  }

  // The initialization procedure of Section 3, with the initialization list {u_i, x^0_i, v_i}.
  void initialize() {
    std::vector<Real> x0(n_);
    for (size_t i = 0; i < n_; ++i) {
      x0[i] = params_.initial_guess ? Real((*params_.initial_guess)[i]) : (u_[i] + v_[i]) / 2;
    }
    init_list_.resize(n_);
    init_f_.resize(n_);
    l0_.resize(n_);
    for (size_t i = 0; i < n_; ++i) {
      auto &list = init_list_[i];
      list = {u_[i], x0[i], v_[i]};
      if (x0[i] == u_[i] || x0[i] == v_[i]) {
        list = {u_[i], (u_[i] + v_[i]) / 2, v_[i]};
      }
      l0_[i] = static_cast<size_t>(std::find(list.begin(), list.end(), x0[i]) - list.begin());
    }
    explicit_box &root = root_;
    root.x = x0;
    // The best point stays x0 if the run is over before f is ever called:
    xbest_ = x0;
    root.f = evaluate(x0);
    root.y.resize(n_);
    for (size_t i = 0; i < n_; ++i) {
      root.y[i] = (x0[i] - u_[i] > v_[i] - x0[i]) ? u_[i] : v_[i];
    }
    root.level = 1;
    root.splits.assign(n_, 0);
    root.xi1.assign(n_, std::numeric_limits<Real>::quiet_NaN());
    root.df1 = root.xi2 = root.df2 = root.xi1;
    boxes_.push_back(box{0, 0, 0, 0, root.f, root.level, 0, 0, 0, 0});
    size_t current = 0;
    for (size_t i = 0; i < n_ && !done(); ++i) {
      load(current);
      init_f_[i] = split_by_list(current, i, &current);
    }
    if (!done()) {
      push(current);
    }
    f0_ = fbest_;
    // The variability of f along each coordinate, from the quadratic interpolants through consecutive triples:
    std::vector<Real> variability(n_, Real(0));
    for (size_t i = 0; i < n_ && !done(); ++i) {
      auto const &a = init_list_[i];
      auto const &f = init_f_[i];
      Real lo = std::numeric_limits<Real>::infinity();
      Real hi = -lo;
      for (size_t l = 0; l + 2 < a.size(); ++l) {
        Real c[3];
        quadratic_through(a[l], f[l], a[l + 1], f[l + 1], a[l + 2], f[l + 2], c);
        // The range over [a_l, a_{l+2}]: the ends, and the stationary point when it lies between them.
        Real const stationary = c[2] != 0 ? -c[1] / (2 * c[2]) : a[l];
        for (Real t : {a[l], a[l + 2], std::clamp<Real>(stationary, a[l], a[l + 2])}) {
          Real const q = c[0] + t * (c[1] + t * c[2]);
          lo = (std::min)(lo, q);
          hi = (std::max)(hi, q);
        }
      }
      using std::isnan;
      variability[i] = isnan(hi - lo) ? std::numeric_limits<Real>::infinity() : hi - lo;
    }
    rank_.resize(n_);
    std::vector<size_t> order(n_);
    for (size_t i = 0; i < n_; ++i) {
      order[i] = i;
    }
    std::stable_sort(order.begin(), order.end(), [&](size_t a, size_t b) { return variability[a] > variability[b]; });
    for (size_t r = 0; r < n_; ++r) {
      rank_[order[r]] = r;
    }
  }

  // One sweep through the levels (Section 3), then the local searches of Section 5.4.
  // Returns false once the run is over.
  bool sweep() {
    if (done()) {
      return false;
    }
    bool any = false;
    for (size_t s = 1; s < smax_ && !done(); ++s) {
      if (heaps_[s].empty()) {
        continue;
      }
      any = true;
      std::pop_heap(heaps_[s].begin(), heaps_[s].end(), box_greater{this});
      size_t const b = heaps_[s].back();
      heaps_[s].pop_back();
      process(b);
    }
    if (params_.local_search && !done()) {
      shopping_basket();
    }
    candidates_.clear();
    return any && !done();
  }

  std::vector<box> const &boxes() const { return boxes_; }
  // Box b with its bounds, recovered by going back in its history: the last split along each coordinate gives the data
  // of that coordinate, and the coordinates never split keep those of the root. It takes one step per split in the
  // history of the box, so fewer than n + s_max.
  void recover(size_t b, explicit_box &out) const {
    out = root_;
    out.f = boxes_[b].f;
    out.level = boxes_[b].level;
    for (size_t k = b; k != 0; k = boxes_[k].parent) {
      box const &B = boxes_[k];
      if (out.splits[B.i]++ == 0) {
        out.x[B.i] = B.x;
        out.y[B.i] = B.y;
        out.xi1[B.i] = B.xi1;
        out.df1[B.i] = B.df1;
        out.xi2[B.i] = B.xi2;
        out.df2[B.i] = B.df2;
      }
    }
  }
  // The local minimizers found by local search (the shopping basket of Section 5.4):
  std::vector<std::pair<std::vector<Real>, Real>> const &basket() const { return basket_; }
  std::vector<Real> const &best_point() const { return xbest_; }
  Real best_value() const { return fbest_; }
  size_t function_calls() const { return function_calls_; }
  size_t max_level() const { return smax_; }

  ArgumentContainer best_argument() const { return to_argument(xbest_); }

private:
  struct box_greater {
    multilevel_coordinate_search_engine const *engine;
    bool operator()(size_t a, size_t b) const {
      Real const fa = engine->boxes_[a].f;
      Real const fb = engine->boxes_[b].f;
      return fa > fb || (fa == fb && a > b);
    }
  };

  static Real golden() { return Real(0.6180339887498948482045868343656381177203L); }

  bool done() const {
    using std::isnan;
    if (cancellation_ && *cancellation_) {
      return true;
    }
    if (!isnan(target_value_) && fbest_ <= Real(target_value_)) {
      return true;
    }
    return function_calls_ >= max_function_calls_;
  }

  ArgumentContainer to_argument(std::vector<Real> const &x) const {
    ArgumentContainer arg;
    if constexpr (detail::has_resize_v<ArgumentContainer>) {
      arg.resize(n_);
    }
    for (size_t i = 0; i < n_; ++i) {
      arg[i] = x[i];
    }
    return arg;
  }

  // Calls the cost function (NaN counts as +infinity), or returns +infinity once the run is over.
  Real evaluate(std::vector<Real> const &x) {
    using std::isnan;
    if (done()) {
      return std::numeric_limits<Real>::infinity();
    }
    auto const arg = to_argument(x);
    ResultType const value = cost_function_(arg);
    ++function_calls_;
    if (queries_) {
      queries_->push_back(std::make_pair(arg, value));
    }
    Real f = static_cast<Real>(value);
    if (isnan(f)) {
      f = std::numeric_limits<Real>::infinity();
    }
    if (f < fbest_ || function_calls_ == 1) {
      fbest_ = f;
      xbest_ = x;
      if (current_minimum_cost_) {
        *current_minimum_cost_ = static_cast<ResultType>(f);
      }
    }
    return f;
  }

  // The coefficients c of the quadratic c0 + c1 t + c2 t^2 through three points with distinct abscissae.
  static void quadratic_through(Real a, Real fa, Real b, Real fb, Real c, Real fc, Real *coef) {
    Real const dab = (fb - fa) / (b - a);
    Real const dabc = ((fc - fb) / (c - b) - dab) / (c - a);
    coef[2] = dabc;
    coef[1] = dab - dabc * (a + b);
    coef[0] = fa - a * (dab - dabc * b);
  }

  // The minimizer over [lo, hi] (in either order) of the quadratic with coefficients coef.
  static Real vertex(Real const *coef, Real lo, Real hi) {
    if (lo > hi) {
      std::swap(lo, hi);
    }
    auto q = [&](Real t) { return coef[0] + t * (coef[1] + t * coef[2]); };
    Real best = q(lo) <= q(hi) ? lo : hi;
    if (coef[2] > 0) {
      Real const t = -coef[1] / (2 * coef[2]);
      if (lo < t && t < hi && q(t) < q(best)) {
        best = t;
      }
    }
    return best;
  }

  // subint of (4), as in the authors' Matlab code: the second case is for 1000|x| >= 1.
  static Real subint(Real x, Real y) {
    using std::abs;
    Real const s = y < 0 ? Real(-1) : Real(1);
    if (1000 * abs(x) < 1) {
      return abs(y) > 1000 ? s : y;
    }
    return abs(y) > 1000 * abs(x) ? 10 * s * abs(x) : y;
  }

  void push(size_t b) {
    if (boxes_[b].level >= smax_) {
      boxes_[b].level = smax_;
      candidates_.push_back(b);
      return;
    }
    heaps_[boxes_[b].level].push_back(b);
    std::push_heap(heaps_[boxes_[b].level].begin(), heaps_[boxes_[b].level].end(), box_greater{this});
  }

  // Recovers box b into box_, the box being split.
  explicit_box const &load(size_t b) {
    recover(b, box_);
    loaded_ = b;
    return box_;
  }

  // A child of box `parent` split along coordinate i, with base point coordinate xi (value f), the other end `end` of
  // its ith interval, and the points `others` (coordinate, value) of this split, nearest first. The differences kept
  // for the other coordinates stay relative to the base points they were measured from: the model is separable.
  // The parent is box_.
  size_t add_child(size_t parent, size_t i, Real xi, Real f, Real end, size_t level,
                   std::vector<std::pair<Real, Real>> others) {
    BOOST_MATH_ASSERT(loaded_ == parent);
    box child{parent, i, xi, end, f, (std::min)(level, smax_), box_.xi1[i], box_.df1[i], box_.xi2[i], box_.df2[i]};
    // The two points closest in the history, other than the base point:
    using std::abs;
    using std::isnan;
    if (!isnan(box_.xi1[i])) {
      others.emplace_back(box_.xi1[i], box_.f + box_.df1[i]);
    }
    if (!isnan(box_.xi2[i])) {
      others.emplace_back(box_.xi2[i], box_.f + box_.df2[i]);
    }
    std::stable_sort(others.begin(), others.end(),
                     [&](auto const &a, auto const &b) { return abs(a.first - xi) < abs(b.first - xi); });
    size_t k = 0;
    for (auto const &o : others) {
      if (o.first == xi || (k == 1 && o.first == child.xi1)) {
        continue;
      }
      (k == 0 ? child.xi1 : child.xi2) = o.first;
      (k == 0 ? child.df1 : child.df2) = o.second - f;
      if (++k == 2) {
        break;
      }
    }
    boxes_.push_back(child);
    return boxes_.size() - 1;
  }

  // Splits box b along coordinate i (never split before) at the initialization list and the golden section points
  // between them (Section 3). Returns the function values at the list. When `current` is given (the initialization
  // procedure), it receives the child to split next, and that child is not pushed. Box b is box_.
  std::vector<Real> split_by_list(size_t b, size_t i, size_t *current = nullptr) {
    BOOST_MATH_ASSERT(loaded_ == b);
    auto const &a = init_list_[i];
    size_t const L = a.size();
    std::vector<Real> f(L);
    x_ = box_.x;
    for (size_t l = 0; l < L; ++l) {
      if (l == l0_[i]) {
        f[l] = box_.f;
      }
      else {
        x_[i] = a[l];
        f[l] = evaluate(x_);
      }
    }
    if (done()) {
      return f;
    }
    size_t const s = box_.level;
    boxes_[b].level = 0;
    auto others = [&](size_t l) {
      std::vector<std::pair<Real, Real>> o;
      std::vector<size_t> idx;
      for (size_t k = 0; k < L; ++k) {
        if (k != l) {
          idx.push_back(k);
        }
      }
      std::stable_sort(idx.begin(), idx.end(), [&](size_t p, size_t q) {
        using std::abs;
        return abs(a[p] - a[l]) < abs(a[q] - a[l]);
      });
      for (size_t k : idx) {
        o.emplace_back(a[k], f[k]);
      }
      return o;
    };
    size_t lbest = 0;
    for (size_t l = 1; l < L; ++l) {
      if (f[l] < f[lbest]) {
        lbest = l;
      }
    }
    // The children whose base point is the best point, left and right of it:
    size_t left = 0, right = 0;
    bool has_left = false, has_right = false;
    std::vector<size_t> children;
    if (a[0] > u_[i]) {
      children.push_back(add_child(b, i, a[0], f[0], u_[i], s + 1, others(0)));
      if (lbest == 0) {
        left = children.back();
        has_left = true;
      }
    }
    for (size_t l = 1; l < L; ++l) {
      // The part with the smaller function value gets the larger fraction:
      bool const left_larger = f[l - 1] <= f[l];
      Real const z = a[l - 1] + (left_larger ? golden() : golden() * golden()) * (a[l] - a[l - 1]);
      children.push_back(add_child(b, i, a[l - 1], f[l - 1], z, left_larger ? s + 1 : s + 2, others(l - 1)));
      if (l - 1 == lbest) {
        right = children.back();
        has_right = true;
      }
      children.push_back(add_child(b, i, a[l], f[l], z, left_larger ? s + 2 : s + 1, others(l)));
      if (l == lbest) {
        left = children.back();
        has_left = true;
      }
    }
    if (a[L - 1] < v_[i]) {
      children.push_back(add_child(b, i, a[L - 1], f[L - 1], v_[i], s + 1, others(L - 1)));
      if (lbest == L - 1) {
        right = children.back();
        has_right = true;
      }
    }
    size_t next = has_left ? left : right;
    if (has_left && has_right) {
      // The one containing the minimizer of the quadratic model through the best point and its neighbors:
      size_t const l = (std::min)((std::max)(lbest, size_t(1)), L - 2);
      Real c[3];
      quadratic_through(a[l - 1], f[l - 1], a[l], f[l], a[l + 1], f[l + 1], c);
      next = vertex(c, a[0], a[L - 1]) < a[lbest] ? left : right;
    }
    for (size_t c : children) {
      if (current && c == next) {
        continue;
      }
      push(c);
    }
    if (current) {
      *current = next;
    }
    return f;
  }

  // Splits box b along coordinate i at z (new value fz) and at the golden section point between x_i and z, giving
  // two parts, or three when z != y_i (Sections 4.1 and 4.2). Box b is box_.
  void split_at(size_t b, size_t i, Real z, Real fz, bool by_rank) {
    BOOST_MATH_ASSERT(loaded_ == b);
    using std::abs;
    explicit_box const &B = box_;
    Real const xi = B.x[i];
    Real const yi = B.y[i];
    Real const f = B.f;
    size_t const s = B.level;
    bool const near_larger = f <= fz;
    Real const fraction = near_larger ? golden() : golden() * golden();
    Real const g = xi + fraction * (z - xi);
    if (g == xi || g == z) {
      // Too narrow to split in floating point; let the box move up a level instead.
      boxes_[b].level = s + 1;
      push(b);
      return;
    }
    boxes_[b].level = 0;
    std::vector<std::pair<Real, Real>> from_x{{z, fz}};
    std::vector<std::pair<Real, Real>> from_z{{xi, f}};
    push(add_child(b, i, xi, f, g, near_larger ? s + 1 : s + 2, from_x));
    push(add_child(b, i, z, fz, g, near_larger ? s + 2 : s + 1, from_z));
    if (z != yi) {
      size_t level = s + 1;
      if (!by_rank && !(abs(yi - z) > (std::min)(abs(g - xi), abs(z - g)))) {
        level = s + 2;
      }
      push(add_child(b, i, z, fz, yi, level, from_z));
    }
  }

  // Step 2 of a sweep (Section 3) for box b: split it (Section 4), or move it up a level.
  void process(size_t b) {
    using std::abs;
    using std::isfinite;
    explicit_box const &B = load(b);
    size_t const s = B.level;
    size_t const m = *std::min_element(B.splits.begin(), B.splits.end());
    if (s > 2 * n_ * (m + 1)) {
      // Splitting by rank (Section 4.1):
      size_t i = n_;
      for (size_t j = 0; j < n_; ++j) {
        if (B.splits[j] == m && (i == n_ || rank_[j] < rank_[i])) {
          i = j;
        }
      }
      if (m == 0) {
        split_by_list(b, i);
        return;
      }
      Real const z = B.x[i] + 2 * (subint(B.x[i], B.y[i]) - B.x[i]) / 3;
      x_ = B.x;
      x_[i] = z;
      Real const fz = evaluate(x_);
      if (!done()) {
        split_at(b, i, z, fz, true);
      }
      return;
    }
    // Splitting by expected gain (Section 4.2):
    Real best_gain = std::numeric_limits<Real>::infinity();
    size_t i = n_;
    Real z = 0;
    for (size_t j = 0; j < n_; ++j) {
      Real gain;
      Real zj = 0;
      if (B.splits[j] == 0) {
        gain = *std::min_element(init_f_[j].begin(), init_f_[j].end()) - init_f_[j][l0_[j]];
      }
      else {
        Real const xj = B.x[j];
        Real const outer = subint(xj, B.y[j]);
        Real const inner = xj + (outer - xj) / 10;
        // e_j(t) = alpha (t - x_j) + beta (t - x_j)^2 through the two more points:
        Real const h1 = B.xi1[j] - xj;
        Real const h2 = B.xi2[j] - xj;
        Real const d1 = B.df1[j] / h1;
        Real const d2 = B.df2[j] / h2;
        Real const beta = (d1 - d2) / (h1 - h2);
        Real const alpha = d1 - beta * h1;
        Real const c[3] = {0, alpha, beta};
        Real const t = vertex(c, inner - xj, outer - xj);
        gain = alpha * t + beta * t * t;
        zj = xj + t;
        if (!isfinite(gain) || outer == xj) {
          gain = std::numeric_limits<Real>::infinity();
        }
      }
      if (gain < best_gain) {
        best_gain = gain;
        i = j;
        z = zj;
      }
    }
    // (5): split only if the model expects to beat the best value so far.
    if (i < n_ && B.f + best_gain < fbest_) {
      if (B.splits[i] == 0) {
        split_by_list(b, i);
        return;
      }
      x_ = B.x;
      x_[i] = z;
      Real const fz = evaluate(x_);
      if (!done()) {
        split_at(b, i, z, fz, false);
      }
      return;
    }
    boxes_[b].level = s + 1;
    push(b);
  }

  // ---------------------------------------------------------------------------------------------------------------
  // Local search (Section 5).

  std::vector<Real> along(std::vector<Real> const &x, std::vector<Real> const &p, Real alpha) const {
    std::vector<Real> y(n_);
    for (size_t i = 0; i < n_; ++i) {
      y[i] = std::clamp<Real>(x[i] + alpha * p[i], u_[i], v_[i]);
    }
    return y;
  }

  // A univariate line search along x + alpha p, alpha in [amin, amax], which stands in for gls. `samples` holds the
  // (alpha, f) pairs known so far, including alpha = 0, and receives the new ones, sorted by alpha. With only one
  // sample, the first step is `step`. It extrapolates until the best point is bracketed, then interpolates
  // parabolically (golden section steps as a safeguard), using at most `max_points` points in all.
  void line_search(std::vector<Real> const &x, std::vector<Real> const &p, Real amin, Real amax,
                   std::vector<std::pair<Real, Real>> &samples, size_t max_points, Real step) {
    using std::abs;
    using std::sqrt;
    Real const tol = sqrt(std::numeric_limits<Real>::epsilon());
    auto by_alpha = [](auto const &a, auto const &b) { return a.first < b.first; };
    std::sort(samples.begin(), samples.end(), by_alpha);
    while (samples.size() < max_points && !done()) {
      size_t const N = samples.size();
      size_t ib = 0;
      for (size_t k = 1; k < N; ++k) {
        if (samples[k].second < samples[ib].second) {
          ib = k;
        }
      }
      Real const scale = (std::max)(abs(samples.front().first), abs(samples.back().first));
      Real a;
      if (N == 1) {
        a = std::clamp<Real>(step, amin, amax);
        if (a == samples[0].first) {
          a = std::clamp<Real>(-step, amin, amax);
        }
      }
      else if (ib == 0 && samples[0].first > amin) {
        a = (std::max)(amin, samples[0].first - 2 * (samples[1].first - samples[0].first));
      }
      else if (ib == N - 1 && samples[N - 1].first < amax) {
        a = (std::min)(amax, samples[N - 1].first + 2 * (samples[N - 1].first - samples[N - 2].first));
      }
      else {
        // Bracketed, or the best point is on the bound: refine around it.
        size_t const l = ib == 0 ? 0 : (ib == N - 1 ? N - 2 : ib - 1);
        size_t const r = ib == 0 ? 1 : (ib == N - 1 ? N - 1 : ib + 1);
        Real const al = samples[l].first;
        Real const ar = samples[r].first;
        Real const ab = samples[ib].first;
        if (ar - al <= tol * (1 + scale)) {
          break;
        }
        a = std::numeric_limits<Real>::quiet_NaN();
        if (ib != 0 && ib != N - 1) {
          Real c[3];
          quadratic_through(al, samples[l].second, ab, samples[ib].second, ar, samples[r].second, c);
          if (c[2] > 0) {
            a = -c[1] / (2 * c[2]);
          }
        }
        Real const min_gap = tol * (1 + scale);
        if (!(al + min_gap < a && a < ar - min_gap && abs(a - ab) > min_gap)) {
          // Golden section step into the larger side of the best point:
          a = ab - al > ar - ab ? ab - (1 - golden()) * (ab - al) : ab + (1 - golden()) * (ar - ab);
        }
      }
      bool known = false;
      for (auto const &sample : samples) {
        if (abs(sample.first - a) <= tol * (1 + scale) / 4) {
          known = true;
        }
      }
      if (known) {
        break;
      }
      samples.emplace_back(a, evaluate(along(x, p, a)));
      std::sort(samples.begin(), samples.end(), by_alpha);
    }
  }

  // The local quadratic model of Section 5.1 around xb: f(xb + h) ~ fb + g^T h + h^T G h/2.
  struct model {
    std::vector<Real> xb;
    Real fb;
    std::vector<Real> g;
    std::vector<Real> G;
  };

  // Triple search (Section 5.1) over the coordinates `coords`, in order. For each coordinate i, `triple(i, M)` gives
  // the two other values t1, t2 of the ith coordinate (xb_i being the third) and f at xb with x_i changed to them.
  // It returns `skip` for a coordinate with no usable triple, and `stop` to end the search.
  // Without `full`, the off-diagonal entries of G are those of the previous model.
  enum class triple_status { ok, skip, stop };
  template <class Triple>
  void triple_search(model &M, std::vector<size_t> const &coords, bool full, Triple triple) {
    size_t const n = n_;
    // The triple of each coordinate: its value at the time and the two others.
    std::vector<Real> t0(n), t1(n), t2(n);
    std::vector<size_t> done_coords;
    for (size_t i : coords) {
      Real f1, f2;
      triple_status const status = triple(i, M, t1[i], f1, t2[i], f2);
      if (status == triple_status::stop) {
        return;
      }
      if (status == triple_status::skip) {
        continue;
      }
      t0[i] = M.xb[i];
      std::vector<Real> xnew = M.xb;
      Real fnew = M.fb;
      Real const h1 = t1[i] - M.xb[i];
      Real const h2 = t2[i] - M.xb[i];
      Real const d1 = (f1 - M.fb) / h1;
      Real const d2 = (f2 - M.fb) / h2;
      M.G[i * n + i] = 2 * (d1 - d2) / (h1 - h2);
      M.g[i] = d1 - M.G[i * n + i] * h1 / 2;
      if (f1 < fnew || f2 < fnew) {
        fnew = f1 < f2 ? f1 : f2;
        xnew[i] = f1 < f2 ? t1[i] : t2[i];
      }
      if (full) {
        // x^{ik}: the ith coordinate at the better of t1, t2, the kth at the one with the smaller model value.
        Real const ti = f1 < f2 ? t1[i] : t2[i];
        Real const hi = ti - M.xb[i];
        for (size_t k : done_coords) {
          // The two values of coordinate k's triple other than xb_k, which may have moved within the triple:
          Real ck[3] = {t0[k], t1[k], t2[k]};
          Real candidates[2];
          size_t nc = 0;
          for (Real c : ck) {
            if (c != M.xb[k] && nc < 2 && (nc == 0 || c != candidates[0])) {
              candidates[nc++] = c;
            }
          }
          if (nc == 0) {
            continue;
          }
          auto qk = [&](Real c) {
            Real const h = c - M.xb[k];
            return M.g[k] * h + M.G[k * n + k] * h * h / 2;
          };
          Real const tk = (nc == 1 || qk(candidates[0]) <= qk(candidates[1])) ? candidates[0] : candidates[1];
          Real const hk = tk - M.xb[k];
          std::vector<Real> x = M.xb;
          x[i] = ti;
          x[k] = tk;
          Real const fik = evaluate(x);
          if (done()) {
            return;
          }
          Real const Gik = (fik - M.fb - M.g[k] * hk - M.g[i] * hi - M.G[k * n + k] * hk * hk / 2 -
                            M.G[i * n + i] * hi * hi / 2) / (hk * hi);
          M.G[i * n + k] = M.G[k * n + i] = Gik;
          if (fik < fnew) {
            fnew = fik;
            xnew = x;
          }
        }
      }
      done_coords.push_back(i);
      if (fnew < M.fb) {
        // Re-expand the model around the new best point:
        std::vector<Real> dx(n);
        for (size_t l = 0; l < n; ++l) {
          dx[l] = xnew[l] - M.xb[l];
        }
        for (size_t k = 0; k < n; ++k) {
          for (size_t l = 0; l < n; ++l) {
            M.g[k] += M.G[k * n + l] * dx[l];
          }
        }
        M.xb = xnew;
        M.fb = fnew;
      }
    }
  }

  // The local search of Section 5.3 from x (value f). Returns the best point found and its value.
  std::pair<std::vector<Real>, Real> local_search(std::vector<Real> x, Real f) {
    using std::abs;
    size_t const n = n_;
    std::vector<Real> x0(n);
    for (size_t i = 0; i < n; ++i) {
      x0[i] = std::clamp<Real>(Real(0), u_[i], v_[i]);
    }
    std::vector<Real> d(n);
    for (size_t i = 0; i < n; ++i) {
      d[i] = (std::min)((std::min)(v_[i] - x[i], x[i] - u_[i]), Real(0.25) * (1 + abs(x[i] - x0[i])));
    }
    model M{x, f, std::vector<Real>(n, Real(0)), std::vector<Real>(n * n, Real(0))};
    std::vector<size_t> all(n);
    for (size_t i = 0; i < n; ++i) {
      all[i] = i;
    }
    // Step 1: a full triple search on the triples found by line searches along the coordinates (Section 5.2).
    bool first = true;
    bool failed = false;
    triple_search(M, all, true, [&](size_t i, model &Mi, Real &t1, Real &f1, Real &t2, Real &f2) {
      std::vector<Real> e(n, Real(0));
      e[i] = 1;
      Real const step = d[i] > 0 ? d[i] : Real(0.25) * (1 + abs(Mi.xb[i] - x0[i]));
      std::vector<std::pair<Real, Real>> samples{{Real(0), Mi.fb}};
      line_search(Mi.xb, e, u_[i] - Mi.xb[i], v_[i] - Mi.xb[i], samples, 6, step);
      if (samples.size() < 3) {
        failed = true;
        return triple_status::stop;
      }
      size_t ib = 0, i0 = 0;
      for (size_t k = 0; k < samples.size(); ++k) {
        if (samples[k].second < samples[ib].second) {
          ib = k;
        }
        if (samples[k].first == 0) {
          i0 = k;
        }
      }
      size_t const N = samples.size();
      size_t pick[2];
      if (first || ib == i0) {
        // The best point and its nearest neighbors on both sides, or the two nearest on one side:
        size_t const l = ib == 0 ? 0 : (ib == N - 1 ? N - 3 : ib - 1);
        size_t trio[3] = {l, l + 1, l + 2};
        size_t k = 0;
        for (size_t c : trio) {
          if (c != ib) {
            pick[k++] = c;
          }
        }
        if (first) {
          // Nothing has been fitted yet, so the model can move to the best point.
          Mi.xb[i] += samples[ib].first;
          Mi.fb = samples[ib].second;
          Real const shift = samples[ib].first;
          t1 = Mi.xb[i] + samples[pick[0]].first - shift;
          t2 = Mi.xb[i] + samples[pick[1]].first - shift;
          f1 = samples[pick[0]].second;
          f2 = samples[pick[1]].second;
          first = false;
          return triple_status::ok;
        }
      }
      else {
        // The best point, the start of the line search, and the best point's neighbor on the other side:
        pick[0] = ib;
        if (ib > i0) {
          pick[1] = ib + 1 < N ? ib + 1 : (i0 > 0 ? i0 - 1 : i0 + 1);
        }
        else {
          pick[1] = ib > 0 ? ib - 1 : (i0 + 1 < N ? i0 + 1 : i0 - 1);
        }
        if (pick[1] == ib || pick[1] == i0) {
          for (size_t k = 0; k < N; ++k) {
            if (k != ib && k != i0) {
              pick[1] = k;
              break;
            }
          }
        }
      }
      t1 = Mi.xb[i] + samples[pick[0]].first;
      t2 = Mi.xb[i] + samples[pick[1]].first;
      f1 = samples[pick[0]].second;
      f2 = samples[pick[1]].second;
      return triple_status::ok;
    });
    if (failed || done()) {
      return {M.xb, M.fb};
    }
    // Step 2 starts from the point Step 1 found:
    for (size_t i = 0; i < n; ++i) {
      d[i] = (std::min)((std::min)(v_[i] - M.xb[i], M.xb[i] - u_[i]), Real(0.25) * (1 + abs(M.xb[i] - x0[i])));
    }
    using std::cbrt;
    Real const delta = cbrt(std::numeric_limits<Real>::epsilon());
    // Steps 2 and 6: minimize the model over the box [-d, d] intersected with [u - x, v - x], then search along the
    // minimizer. Returns r, the ratio of the actual to the predicted gain.
    auto quadratic_step = [&](std::vector<bool> const &free) {
      std::vector<Real> lo(n), hi(n);
      for (size_t i = 0; i < n; ++i) {
        lo[i] = free[i] ? (std::max)(-d[i], u_[i] - M.xb[i]) : Real(0);
        hi[i] = free[i] ? (std::min)(d[i], v_[i] - M.xb[i]) : Real(0);
      }
      std::vector<Real> p = minimize_quadratic_on_box(M.g, M.G, lo, hi);
      Real predicted = 0;
      bool moves = false;
      for (size_t i = 0; i < n; ++i) {
        Real Gp = 0;
        for (size_t j = 0; j < n; ++j) {
          Gp += M.G[i * n + j] * p[j];
        }
        predicted += p[i] * (M.g[i] + Gp / 2);
        moves = moves || p[i] != 0;
      }
      if (!moves || !(predicted < 0)) {
        return Real(0);
      }
      Real amin = -std::numeric_limits<Real>::infinity();
      Real amax = std::numeric_limits<Real>::infinity();
      for (size_t i = 0; i < n; ++i) {
        if (p[i] > 0) {
          amax = (std::min)(amax, (v_[i] - M.xb[i]) / p[i]);
          amin = (std::max)(amin, (u_[i] - M.xb[i]) / p[i]);
        }
        else if (p[i] < 0) {
          amax = (std::min)(amax, (u_[i] - M.xb[i]) / p[i]);
          amin = (std::max)(amin, (v_[i] - M.xb[i]) / p[i]);
        }
      }
      amax = (std::max)(amax, Real(1));
      amin = (std::min)(amin, Real(0));
      Real const fold = M.fb;
      std::vector<std::pair<Real, Real>> samples{{Real(0), fold}, {Real(1), evaluate(along(M.xb, p, 1))}};
      line_search(M.xb, p, amin, amax, samples, params_.line_search_points, Real(1));
      size_t ib = 0;
      for (size_t k = 1; k < samples.size(); ++k) {
        if (samples[k].second < samples[ib].second) {
          ib = k;
        }
      }
      if (samples[ib].second < fold) {
        // Re-expand the gradient of the model around the new point:
        auto const xnew = along(M.xb, p, samples[ib].first);
        for (size_t k = 0; k < n; ++k) {
          for (size_t l = 0; l < n; ++l) {
            M.g[k] += M.G[k * n + l] * (xnew[l] - M.xb[l]);
          }
        }
        M.xb = xnew;
        M.fb = samples[ib].second;
      }
      return (fold - M.fb) / -predicted;
    };
    std::vector<bool> free(n, true);
    Real r = quadratic_step(free);
    Real f_step_start = f;
    std::vector<Real> xold = x;
    bool last_full = true;
    for (size_t step = 0; step < params_.max_local_search_steps && !done(); ++step) {
      // Step 3:
      bool improved = M.fb < f_step_start;
      Real gx = 0;
      for (size_t i = 0; i < n; ++i) {
        gx += abs(M.g[i]) * (std::max)(abs(M.xb[i]), abs(xold[i]));
      }
      bool const criterion = !improved || gx < Real(1e-18) * (f0_ - M.fb);
      std::vector<size_t> inner;
      std::vector<size_t> at_bound;
      for (size_t i = 0; i < n; ++i) {
        free[i] = u_[i] < M.xb[i] && M.xb[i] < v_[i];
        (free[i] ? inner : at_bound).push_back(i);
      }
      if (at_bound.empty() && last_full && criterion) {
        break;
      }
      f_step_start = M.fb;
      xold = M.xb;
      // Step 4: try to get away from the boundary.
      if (!at_bound.empty() && criterion) {
        for (size_t i : at_bound) {
          std::vector<Real> e(n, Real(0));
          e[i] = 1;
          std::vector<std::pair<Real, Real>> samples{{Real(0), M.fb}};
          Real const step0 = Real(0.25) * (1 + abs(M.xb[i] - x0[i]));
          line_search(M.xb, e, u_[i] - M.xb[i], v_[i] - M.xb[i], samples, params_.line_search_points, step0);
          for (auto const &sample : samples) {
            if (sample.second < M.fb) {
              M.fb = sample.second;
              M.xb[i] = std::clamp<Real>(xold[i] + sample.first, u_[i], v_[i]);
            }
          }
        }
        if (!(M.fb < f_step_start)) {
          break;
        }
        inner.clear();
        for (size_t i = 0; i < n; ++i) {
          free[i] = u_[i] < M.xb[i] && M.xb[i] < v_[i];
          if (free[i]) {
            inner.push_back(i);
          }
        }
      }
      // Step 5: a triple search at distance delta in the free coordinates.
      last_full = abs(r - 1) > Real(0.25) || criterion;
      triple_search(M, inner, last_full, [&](size_t i, model &Mi, Real &t1, Real &f1, Real &t2, Real &f2) {
        Real const h = delta * (std::max)(abs(Mi.xb[i]), Real(1));
        if (Mi.xb[i] - h < u_[i]) {
          t1 = Mi.xb[i] + h;
          t2 = Mi.xb[i] + 2 * h;
        }
        else if (Mi.xb[i] + h > v_[i]) {
          t1 = Mi.xb[i] - h;
          t2 = Mi.xb[i] - 2 * h;
        }
        else {
          t1 = Mi.xb[i] - h;
          t2 = Mi.xb[i] + h;
        }
        t1 = std::clamp<Real>(t1, u_[i], v_[i]);
        t2 = std::clamp<Real>(t2, u_[i], v_[i]);
        if (t1 == Mi.xb[i] || t2 == Mi.xb[i] || t1 == t2) {
          return triple_status::skip;
        }
        std::vector<Real> y = Mi.xb;
        y[i] = t1;
        f1 = evaluate(y);
        y[i] = t2;
        f2 = evaluate(y);
        return done() ? triple_status::stop : triple_status::ok;
      });
      if (done()) {
        break;
      }
      // Step 6:
      if (r < Real(0.25)) {
        for (auto &di : d) {
          di /= 2;
        }
      }
      else if (r > Real(0.75)) {
        for (auto &di : d) {
          di *= 2;
        }
      }
      r = quadratic_step(free);
    }
    return {M.xb, M.fb};
  }

  Real distance2(std::vector<Real> const &a, std::vector<Real> const &b) const {
    Real s = 0;
    for (size_t i = 0; i < n_; ++i) {
      s += (a[i] - b[i]) * (a[i] - b[i]);
    }
    return s;
  }

  std::vector<Real> between(std::vector<Real> const &x, std::vector<Real> const &w, Real t) const {
    std::vector<Real> y(n_);
    for (size_t i = 0; i < n_; ++i) {
      y[i] = x[i] + t * (w[i] - x[i]);
    }
    return y;
  }

  // The basket points sorted by distance to x.
  std::vector<size_t> basket_by_distance(std::vector<Real> const &x) const {
    std::vector<size_t> order(basket_.size());
    for (size_t k = 0; k < order.size(); ++k) {
      order[k] = k;
    }
    std::stable_sort(order.begin(), order.end(), [&](size_t a, size_t b) {
      return distance2(basket_[a].first, x) < distance2(basket_[b].first, x);
    });
    return order;
  }

  // Section 5.4: the candidates of this sweep, by ascending function value, start local searches unless they seem to
  // lie in the basin of attraction of a point already in the shopping basket.
  void shopping_basket() {
    std::vector<std::pair<std::vector<Real>, Real>> candidates;
    for (size_t b : candidates_) {
      candidates.emplace_back(load(b).x, boxes_[b].f);
    }
    std::stable_sort(candidates.begin(), candidates.end(),
                     [](auto const &a, auto const &b) { return a.second < b.second; });
    for (auto &candidate : candidates) {
      auto const base_point = candidate.first;
      auto &x = candidate.first;
      Real &fx = candidate.second;
      // Step 1:
      if (std::find(started_.begin(), started_.end(), base_point) != started_.end()) {
        continue;
      }
      bool discard = false;
      for (size_t k : basket_by_distance(x)) {
        auto const &w = basket_[k].first;
        Real const fw = basket_[k].second;
        if (done()) {
          return;
        }
        if (fw > fx || w == x) {
          continue;
        }
        // Step 2:
        auto const x1 = between(x, w, Real(1) / 3);
        Real const f1 = evaluate(x1);
        if (f1 > fx) {
          continue;
        }
        // Step 3:
        auto const x2 = between(x, w, Real(2) / 3);
        Real const f2 = evaluate(x2);
        if (f2 > (std::max)(f1, fw)) {
          if (f1 < fx) {
            x = x1;
            fx = f1;
          }
        }
        else if ((std::min)(f1, f2) < fw) {
          x = f1 <= f2 ? x1 : x2;
          fx = (std::min)(f1, f2);
        }
        else {
          discard = true;
          break;
        }
      }
      if (discard || done()) {
        continue;
      }
      // Step 4:
      started_.push_back(base_point);
      auto result = local_search(x, fx);
      bool is_new = true;
      for (size_t k : basket_by_distance(result.first)) {
        auto &w = basket_[k];
        if (done()) {
          return;
        }
        bool barrier = false;
        if (w.first != result.first) {
          Real const top = (std::max)(result.second, w.second);
          barrier = evaluate(between(result.first, w.first, Real(1) / 3)) > top ||
                    evaluate(between(result.first, w.first, Real(2) / 3)) > top;
        }
        if (!barrier) {
          // The same minimizer as w: keep the better of the two.
          is_new = false;
          if (result.second < w.second) {
            w = result;
          }
          break;
        }
      }
      if (is_new) {
        basket_.push_back(std::move(result));
      }
    }
  }

  Func const &cost_function_;
  multilevel_coordinate_search_parameters<ArgumentContainer> const &params_;
  ResultType target_value_;
  std::atomic<bool> *cancellation_;
  std::atomic<ResultType> *current_minimum_cost_;
  std::vector<std::pair<ArgumentContainer, ResultType>> *queries_;
  size_t n_;
  size_t smax_;
  size_t max_function_calls_;
  size_t function_calls_ = 0;
  std::vector<Real> u_, v_;
  std::vector<Real> xbest_;
  Real fbest_ = std::numeric_limits<Real>::infinity();
  // The smallest function value found in the initialization procedure:
  Real f0_ = std::numeric_limits<Real>::infinity();
  std::vector<std::vector<Real>> init_list_, init_f_;
  std::vector<size_t> l0_;
  // rank_[i] = 0 for the coordinate with the highest variability (pi_i - 1 in the paper):
  std::vector<size_t> rank_;
  std::vector<box> boxes_;
  // The root as an explicit box, and the box being split, as load() recovered it from box loaded_:
  explicit_box root_, box_;
  size_t loaded_ = 0;
  // Scratch space for the points to evaluate:
  std::vector<Real> x_;
  // The non-split boxes of each level below s_max, as min-heaps on f:
  std::vector<std::vector<size_t>> heaps_;
  // Boxes that reached level s_max in this sweep:
  std::vector<size_t> candidates_;
  std::vector<std::vector<Real>> started_;
  std::vector<std::pair<std::vector<Real>, Real>> basket_;
};

} // namespace detail

BOOST_MATH_EXPORT template <typename ArgumentContainer, class Func>
ArgumentContainer multilevel_coordinate_search(
    const Func cost_function, multilevel_coordinate_search_parameters<ArgumentContainer> const &params,
    std::invoke_result_t<Func, ArgumentContainer> target_value =
        std::numeric_limits<std::invoke_result_t<Func, ArgumentContainer>>::quiet_NaN(),
    std::atomic<bool> *cancellation = nullptr,
    std::atomic<std::invoke_result_t<Func, ArgumentContainer>> *current_minimum_cost = nullptr,
    std::vector<std::pair<ArgumentContainer, std::invoke_result_t<Func, ArgumentContainer>>> *queries = nullptr) {
  detail::multilevel_coordinate_search_engine<ArgumentContainer, Func> engine(cost_function, params, target_value, cancellation,
                                                     current_minimum_cost, queries);
  engine.initialize();
  while (engine.sweep()) {
  }
  return engine.best_argument();
}

} BOOST_MATH_NAMESPACE_END // namespace boost::math::optimization
#endif
