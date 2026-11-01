/*
 * Copyright Nick Thompson, 2026
 * Use, modification and distribution are subject to the
 * Boost Software License, Version 1.0. (See accompanying file
 * LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef BOOST_MATH_OPTIMIZATION_DYCORS_HPP
#define BOOST_MATH_OPTIMIZATION_DYCORS_HPP

#include <boost/math/tools/config.hpp>
#ifndef BOOST_MATH_BUILD_MODULE
#include <atomic>
#endif
#include <boost/math/optimization/detail/common.hpp>
#ifndef BOOST_MATH_BUILD_MODULE
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <random>
#include <sstream>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>
#endif

BOOST_MATH_NAMESPACE_BEGIN namespace optimization {

// DYCORS-LMSRBF of: Regis, R. G., Shoemaker, C. A. (2013). Combining radial basis function surrogates and dynamic
// coordinate search in high-dimensional expensive black-box optimization. Engineering Optimization 45(5), 529-555.

BOOST_MATH_EXPORT template <typename ArgumentContainer> struct dycors_parameters {
  using Real = typename ArgumentContainer::value_type;
  ArgumentContainer lower_bounds;
  ArgumentContainer upper_bounds;
  // Evaluated first, in addition to the initial design.
  ArgumentContainer const *initial_guess = nullptr;
  // 0 means 100 times the dimension.
  size_t max_function_calls = 0;
  // n_0, the points of the symmetric Latin hypercube design; 0 means 2(d + 1). At least 2d.
  size_t initial_points = 0;
  // The candidate points per iteration; 0 means min(100d, 5000). Any other number is valid.
  size_t candidates = 0;
};

BOOST_MATH_EXPORT template <typename ArgumentContainer>
void validate_dycors_parameters(dycors_parameters<ArgumentContainer> const &params) {
  std::ostringstream oss;
  detail::validate_bounds(params.lower_bounds, params.upper_bounds);
  for (size_t i = 0; i < params.lower_bounds.size(); ++i) {
    if (!(params.lower_bounds[i] < params.upper_bounds[i])) {
      oss << __FILE__ << ":" << __LINE__ << ":" << __func__;
      oss << ": DYCORS requires lower_bounds[i] < upper_bounds[i], but at i = " << i << " the bounds are "
          << params.lower_bounds[i] << " and " << params.upper_bounds[i] << ".";
      throw std::domain_error(oss.str());
    }
  }
  if (params.initial_guess) {
    detail::validate_initial_guess(*params.initial_guess, params.lower_bounds, params.upper_bounds);
  }
  // The points of a symmetric design come in pairs x, 1 - x, and the linear tail needs d independent pairs:
  if (params.initial_points != 0 && params.initial_points < 2 * params.lower_bounds.size()) {
    oss << __FILE__ << ":" << __LINE__ << ":" << __func__;
    oss << ": The linear tail of the surrogate needs at least 2d = " << 2 * params.lower_bounds.size()
        << " initial points, but initial_points = " << params.initial_points << ".";
    throw std::invalid_argument(oss.str());
  }
  size_t const initial_points = params.initial_points ? params.initial_points : 2 * (params.lower_bounds.size() + 1);
  if (params.max_function_calls != 0 && params.max_function_calls < initial_points) {
    oss << __FILE__ << ":" << __LINE__ << ":" << __func__;
    oss << ": The budget must cover the initial design, but max_function_calls = " << params.max_function_calls
        << " and initial_points = " << initial_points << ".";
    throw std::invalid_argument(oss.str());
  }
}

namespace detail {

// The state of a DYCORS run; dycors() below drives it. Everything is done in double precision on the unit cube.
template <typename ArgumentContainer, class Func, class URBG> class dycors_engine {
public:
  using Real = typename ArgumentContainer::value_type;
  using ResultType = std::invoke_result_t<Func, ArgumentContainer>;

  dycors_engine(Func const &cost_function, dycors_parameters<ArgumentContainer> const &params, URBG &gen,
                ResultType target_value, std::atomic<bool> *cancellation, std::atomic<ResultType> *current_minimum_cost,
                std::vector<std::pair<ArgumentContainer, ResultType>> *queries)
      : cost_function_(cost_function), params_(params), gen_(gen), target_value_(target_value),
        cancellation_(cancellation), current_minimum_cost_(current_minimum_cost), queries_(queries) {
    validate_dycors_parameters(params);
    d_ = params.lower_bounds.size();
    max_function_calls_ = params.max_function_calls ? params.max_function_calls : 100 * d_;
    initial_points_ = params.initial_points ? params.initial_points : 2 * (d_ + 1);
    candidates_ = params.candidates ? params.candidates : (std::min)(100 * d_, size_t(5000));
    failure_tolerance_ = (std::max)(d_, size_t(5));
    stride_ = max_function_calls_;
    Xt_.resize(d_ * stride_);
    if constexpr (detail::has_resize_v<ArgumentContainer>) {
      best_argument_.resize(d_);
    }
    for (size_t i = 0; i < d_; ++i) {
      best_argument_[i] = params.initial_guess ? (*params.initial_guess)[i]
                                               : params.lower_bounds[i] + (params.upper_bounds[i] - params.lower_bounds[i]) / 2;
    }
  }

  ArgumentContainer optimize() {
    initial_design();
    while (!done()) {
      step();
    }
    return best_argument_;
  }

private:
  // Section 3 of the paper, with lengths in units of the sides of the box:
  static constexpr double sigma_init_ = 0.2;
  static constexpr double sigma_min_ = 0.2 / 64;
  static constexpr double sigma_max_ = 1;
  static constexpr size_t success_tolerance_ = 3;
  static constexpr double weights_[4] = {0.3, 0.5, 0.8, 0.95};

  bool done() const {
    using std::isnan;
    if (cancellation_ && *cancellation_) {
      return true;
    }
    if (calls_ > 0 && !isnan(target_value_) && best_value_ <= target_value_) {
      return true;
    }
    if (calls_ > 0 && best_value_ == -std::numeric_limits<ResultType>::infinity()) {
      return true;
    }
    return calls_ >= max_function_calls_;
  }

  ArgumentContainer to_argument(double const *u) const {
    ArgumentContainer x;
    if constexpr (detail::has_resize_v<ArgumentContainer>) {
      x.resize(d_);
    }
    for (size_t i = 0; i < d_; ++i) {
      Real const lo = params_.lower_bounds[i];
      Real const hi = params_.upper_bounds[i];
      Real const xi = lo + static_cast<Real>(u[i]) * (hi - lo);
      x[i] = xi < lo ? lo : (xi > hi ? hi : xi);
    }
    return x;
  }

  // Calls the cost function at x, the next point of X_, and returns its value, with NaN counted as +infinity.
  double evaluate(ArgumentContainer const &x) {
    using std::isnan;
    ResultType const value = cost_function_(x);
    ++calls_;
    if (queries_) {
      queries_->push_back(std::make_pair(x, value));
    }
    ResultType const f = isnan(value) ? std::numeric_limits<ResultType>::infinity() : value;
    if (calls_ == 1 || f < best_value_) {
      best_value_ = f;
      best_argument_ = x;
      best_ = f_.size();
      if (current_minimum_cost_) {
        *current_minimum_cost_ = f;
      }
    }
    return static_cast<double>(f);
  }

  double uniform01() { return std::uniform_real_distribution<double>(0, 1)(gen_); }

  double cubic(double const *x, double const *y) const {
    double r2 = 0;
    for (size_t i = 0; i < d_; ++i) {
      double const t = x[i] - y[i];
      r2 += t * t;
    }
    return r2 * std::sqrt(r2);
  }

  // The symmetric Latin hypercube design of Ye, Li and Sudjianto (2000) with m points in the unit cube, appended to the
  // `first` points of X_: in each coordinate, the points take the centres of the m slices in a random order, with
  // x_{m-1-k} = 1 - x_k. Returns false if the design does not determine a linear polynomial.
  bool symmetric_latin_hypercube(size_t first, size_t m) {
    X_.resize(first * d_);
    X_.resize((first + m) * d_, 0.5);
    double *design = X_.data() + first * d_;
    std::vector<size_t> slices(m / 2);
    for (size_t i = 0; i < d_; ++i) {
      for (size_t k = 0; k < m / 2; ++k) {
        slices[k] = k;
      }
      std::shuffle(slices.begin(), slices.end(), gen_);
      for (size_t k = 0; k < m / 2; ++k) {
        size_t const s = uniform01() < 0.5 ? slices[k] : m - 1 - slices[k];
        design[k * d_ + i] = (static_cast<double>(s) + 0.5) / static_cast<double>(m);
        design[(m - 1 - k) * d_ + i] = 1 - design[k * d_ + i];
      }
    }
    return select_tail_points(first, m);
  }

  // Chooses d + 1 of the m points of X_ from `first` on, T, on which a linear polynomial is determined by its values
  // (by Gaussian elimination with partial pivoting on the rows (1, x)). Sets Pinv_ to the inverse of the matrix of
  // their rows, and PhiT_ to the matrix of phi(t_j, t_k).
  bool select_tail_points(size_t first, size_t m) {
    using std::abs;
    size_t const n = d_ + 1;
    std::vector<double> B(m * n);
    for (size_t k = 0; k < m; ++k) {
      B[k * n] = 1;
      std::copy(X_.data() + (first + k) * d_, X_.data() + (first + k + 1) * d_, B.data() + k * n + 1);
    }
    T_.clear();
    std::vector<bool> taken(m, false);
    for (size_t j = 0; j < n; ++j) {
      size_t p = m;
      for (size_t k = 0; k < m; ++k) {
        if (!taken[k] && (p == m || abs(B[k * n + j]) > abs(B[p * n + j]))) {
          p = k;
        }
      }
      if (p == m || !(abs(B[p * n + j]) > 1e-8)) {
        return false;
      }
      taken[p] = true;
      T_.push_back(first + p);
      for (size_t k = 0; k < m; ++k) {
        if (!taken[k]) {
          double const t = B[k * n + j] / B[p * n + j];
          for (size_t i = j; i < n; ++i) {
            B[k * n + i] -= t * B[p * n + i];
          }
        }
      }
    }
    // Gauss-Jordan elimination of [P_T | I], with partial pivoting:
    std::vector<double> A(n * n);
    Pinv_.assign(n * n, 0.0);
    PhiT_.resize(n * n);
    for (size_t j = 0; j < n; ++j) {
      for (size_t k = 0; k < n; ++k) {
        PhiT_[j * n + k] = cubic(X_.data() + T_[j] * d_, X_.data() + T_[k] * d_);
      }
      A[j * n] = 1;
      std::copy(X_.data() + T_[j] * d_, X_.data() + (T_[j] + 1) * d_, A.data() + j * n + 1);
      Pinv_[j * n + j] = 1;
    }
    for (size_t j = 0; j < n; ++j) {
      size_t p = j;
      for (size_t k = j + 1; k < n; ++k) {
        if (abs(A[k * n + j]) > abs(A[p * n + j])) {
          p = k;
        }
      }
      for (size_t i = 0; i < n; ++i) {
        std::swap(A[j * n + i], A[p * n + i]);
        std::swap(Pinv_[j * n + i], Pinv_[p * n + i]);
      }
      double const inv = 1 / A[j * n + j];
      for (size_t i = 0; i < n; ++i) {
        A[j * n + i] *= inv;
        Pinv_[j * n + i] *= inv;
      }
      for (size_t k = 0; k < n; ++k) {
        double const t = A[k * n + j];
        if (k != j && t != 0) {
          for (size_t i = 0; i < n; ++i) {
            A[k * n + i] -= t * A[j * n + i];
            Pinv_[k * n + i] -= t * Pinv_[j * n + i];
          }
        }
      }
    }
    return true;
  }

  // The surrogate is the cubic radial basis function interpolant with a linear tail,
  //   s(x) = sum_k lambda_k |x - x_k|^3 + a_0 + a^T x,   with sum_k lambda_k (1, x_k) = 0.
  // Its saddle point system is indefinite, but r^3 is conditionally positive definite of order 2, so it reduces to a
  // positive definite one. With l_j the Lagrange basis of linear polynomials on the points t_j of T, the kernel
  //   K(x, y) = phi(x, y) - sum_j l_j(x) phi(t_j, y) - sum_j l_j(y) phi(x, t_j) + sum_jk l_j(x) l_k(y) phi(t_j, t_k)
  // is positive definite on the other points, and s = (the linear interpolant of f on T) + sum_k alpha_k K(., x_k).
  // Its Cholesky factor grows by a row for each new point, in O(N^2) operations, and is backward stable, which
  // an LU factorization of the saddle point system bordered in the same way is not. For each point, lagrange_[k] holds
  // l(x_k), q_[k] holds phi(t_j, x_k), and r_[k] = Phi_T l(x_k) - q(x_k), so that
  // K(x_a, x_b) = phi(x_a, x_b) + l(x_a) . r(x_b) - l(x_b) . q(x_a).
  double kernel(size_t a, size_t b) const {
    size_t const n = d_ + 1;
    double k = cubic(X_.data() + a * d_, X_.data() + b * d_);
    for (size_t j = 0; j < n; ++j) {
      k += lagrange_[a * n + j] * r_[b * n + j] - lagrange_[b * n + j] * q_[a * n + j];
    }
    return k;
  }

  // Adds the point X_[point] to the surrogate. A point that repeats one already in it (numerically) is left out, as
  // the surrogate already interpolates it.
  void add_point(size_t point) {
    size_t const n = d_ + 1;
    double const *x = X_.data() + point * d_;
    for (size_t j = 0; j < n; ++j) {
      double l = Pinv_[j];
      for (size_t i = 0; i < d_; ++i) {
        l += x[i] * Pinv_[(i + 1) * n + j];
      }
      lagrange_.push_back(l);
      q_.push_back(cubic(X_.data() + T_[j] * d_, x));
    }
    for (size_t j = 0; j < n; ++j) {
      double r = -q_[point * n + j];
      for (size_t k = 0; k < n; ++k) {
        r += PhiT_[j * n + k] * lagrange_[point * n + k];
      }
      r_.push_back(r);
    }
    if (std::find(T_.begin(), T_.end(), point) != T_.end()) {
      return;
    }
    // The new row of the Cholesky factor solves L row = (K(x_b, x) for the points b already in it):
    size_t const M = nodes_.size();
    std::vector<double> row(M + 1);
    double const diagonal = kernel(point, point);
    double s = diagonal;
    for (size_t b = 0; b < M; ++b) {
      double t = kernel(nodes_[b], point);
      for (size_t k = 0; k < b; ++k) {
        t -= L_[b][k] * row[k];
      }
      row[b] = t / L_[b][b];
      s -= row[b] * row[b];
    }
    if (!(s > 1e-14 * diagonal)) {
      return;
    }
    row[M] = std::sqrt(s);
    L_.push_back(std::move(row));
    nodes_.push_back(point);
  }

  // The initial guess, if any, exactly as given, then the initial design.
  void initial_design() {
    if (params_.initial_guess && !done()) {
      for (size_t i = 0; i < d_; ++i) {
        Real const lo = params_.lower_bounds[i];
        Real const hi = params_.upper_bounds[i];
        X_.push_back(static_cast<double>(((*params_.initial_guess)[i] - lo) / (hi - lo)));
      }
      f_.push_back(evaluate(*params_.initial_guess));
    }
    size_t const first = f_.size();
    // A random design fails only with probability zero in exact arithmetic:
    while (!symmetric_latin_hypercube(first, initial_points_)) {
    }
    for (size_t k = first; k < first + initial_points_ && !done(); ++k) {
      f_.push_back(evaluate(to_argument(X_.data() + k * d_)));
    }
    n0_ = f_.size();
    if (n0_ < first + initial_points_) {
      return;
    }
    for (size_t k = 0; k < n0_; ++k) {
      add_point(k);
    }
  }

  // The coefficients lambda_k of the surrogate (one for each point evaluated), followed by a_0 and a. As in Regis and
  // Shoemaker (2007), the function values above their median are replaced by the median, so that a few large values
  // do not make the surrogate oscillate; values that are not finite are first replaced by the largest finite one.
  std::vector<double> fit() const {
    using std::isfinite;
    size_t const N = f_.size();
    size_t const n = d_ + 1;
    std::vector<double> y = f_;
    double worst = -std::numeric_limits<double>::infinity();
    for (double v : y) {
      if (isfinite(v)) {
        worst = (std::max)(worst, v);
      }
    }
    for (double &v : y) {
      v = isfinite(v) ? v : (isfinite(worst) ? worst : 0);
    }
    std::vector<double> sorted = y;
    std::nth_element(sorted.begin(), sorted.begin() + static_cast<std::ptrdiff_t>(N / 2), sorted.end());
    double median = sorted[N / 2];
    if (N % 2 == 0) {
      median = (median + *std::max_element(sorted.begin(), sorted.begin() + static_cast<std::ptrdiff_t>(N / 2))) / 2;
    }
    for (double &v : y) {
      v = (std::min)(v, median);
    }
    // K alpha = y - (the linear interpolant of y on T), by the Cholesky factor:
    size_t const M = nodes_.size();
    std::vector<double> alpha(M);
    for (size_t b = 0; b < M; ++b) {
      double t = y[nodes_[b]];
      for (size_t j = 0; j < n; ++j) {
        t -= lagrange_[nodes_[b] * n + j] * y[T_[j]];
      }
      for (size_t k = 0; k < b; ++k) {
        t -= L_[b][k] * alpha[k];
      }
      alpha[b] = t / L_[b][b];
    }
    for (size_t b = M; b-- > 0;) {
      alpha[b] /= L_[b][b];
      for (size_t k = 0; k < b; ++k) {
        alpha[k] -= L_[b][k] * alpha[b];
      }
    }
    // The coefficients of the points of T make the lambda_k orthogonal to linear polynomials, and the tail
    // interpolates what is left of y on T:
    std::vector<double> c(N + n, 0.0);
    for (size_t b = 0; b < M; ++b) {
      c[nodes_[b]] = alpha[b];
      for (size_t j = 0; j < n; ++j) {
        c[T_[j]] -= alpha[b] * lagrange_[nodes_[b] * n + j];
      }
    }
    for (size_t j = 0; j < n; ++j) {
      double g = y[T_[j]];
      for (size_t k = 0; k < N; ++k) {
        g -= c[k] * q_[k * n + j];
      }
      for (size_t i = 0; i < n; ++i) {
        c[N + i] += Pinv_[i * n + j] * g;
      }
    }
    return c;
  }

  // One iteration of Section 3 (Step 2): fit the surrogate, perturb the best point into candidates, evaluate the best
  // candidate by the weighted score, and adapt the step size.
  void step() {
    using std::abs;
    using std::isfinite;
    using std::log;
    using std::sqrt;
    std::vector<double> const coefficients = fit();
    size_t const N = f_.size();
    double const *xb = X_.data() + best_ * d_;

    // The probability of perturbing each coordinate, which decreases from min(20/d, 1) over the budget, but not below
    // 1/d (n_0 counts the initial guess):
    size_t const n0 = n0_;
    double p = (std::min)(20.0 / static_cast<double>(d_), 1.0);
    if (max_function_calls_ > n0 + 1) {
      p *= 1 - log(static_cast<double>(calls_ - n0 + 1)) / log(static_cast<double>(max_function_calls_ - n0));
    }
    p = (std::max)(p, 1 / static_cast<double>(d_));
    // Each candidate differs from the best point b in a few coordinates: candidate k has the indices and values
    // index_[e], value_[e] for offset_[k] <= e < offset_[k + 1].
    size_t const K = candidates_;
    index_.clear();
    value_.clear();
    offset_.assign(1, 0);
    std::normal_distribution<double> normal(0, 1);
    std::uniform_int_distribution<size_t> pick(0, d_ - 1);
    for (size_t k = 0; k < K; ++k) {
      size_t const begin = index_.size();
      for (size_t i = 0; i < d_; ++i) {
        if (uniform01() < p) {
          index_.push_back(i);
          value_.push_back(xb[i] + sigma_ * normal(gen_));
        }
      }
      if (index_.size() == begin) {
        size_t const i = pick(gen_);
        index_.push_back(i);
        value_.push_back(xb[i] + sigma_ * normal(gen_));
      }
      // Steps that leave the box are reflected at its faces until they are inside, as in the paper; the clipping only
      // guards against rounding.
      for (size_t e = begin; e < value_.size(); ++e) {
        double &c = value_[e];
        for (int bounces = 0; (c < 0 || c > 1) && bounces < 64; ++bounces) {
          c = c < 0 ? -c : 2 - c;
        }
        c = (std::min)((std::max)(c, 0.0), 1.0);
      }
      offset_.push_back(index_.size());
    }

    // The surrogate's values and the distances to the nearest evaluated point, then the weighted score (Step 2(d)).
    // With the points stored by coordinates, Xt_[i * stride_ + j], and D_j = |b - x_j|^2, the squared distance from a
    // candidate c to x_j is D_j + sum (c_i - b_i)(c_i + b_i - 2 x_ji) over its perturbed coordinates i.
    for (size_t j = Xt_filled_; j < N; ++j) {
      for (size_t i = 0; i < d_; ++i) {
        Xt_[i * stride_ + j] = X_[j * d_ + i];
      }
    }
    Xt_filled_ = N;
    double const *tail = coefficients.data() + N;
    double tail_b = tail[0];
    D_.assign(N, 0.0);
    for (size_t i = 0; i < d_; ++i) {
      tail_b += tail[i + 1] * xb[i];
      double const *xi = Xt_.data() + i * stride_;
      for (size_t j = 0; j < N; ++j) {
        double const t = xb[i] - xi[j];
        D_[j] += t * t;
      }
    }
    std::vector<double> s(K), delta(K);
    r2_.resize(N);
    for (size_t k = 0; k < K; ++k) {
      std::copy(D_.begin(), D_.end(), r2_.begin());
      double v = tail_b;
      for (size_t e = offset_[k]; e < offset_[k + 1]; ++e) {
        size_t const i = index_[e];
        double const dc = value_[e] - xb[i];
        double const sc = value_[e] + xb[i];
        double const *xi = Xt_.data() + i * stride_;
        for (size_t j = 0; j < N; ++j) {
          r2_[j] += dc * (sc - 2 * xi[j]);
        }
        v += tail[i + 1] * dc;
      }
      // Four independent sums and minima; rounding can make r^2 slightly negative.
      double sum[4] = {0, 0, 0, 0};
      double nearest[4] = {std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity(),
                           std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity()};
      size_t const N4 = N - N % 4;
      size_t j = 0;
      for (; j < N4; j += 4) {
        for (size_t lane = 0; lane < 4; ++lane) {
          double const r2 = (std::max)(r2_[j + lane], 0.0);
          nearest[lane] = (std::min)(nearest[lane], r2);
          sum[lane] += coefficients[j + lane] * r2 * sqrt(r2);
        }
      }
      for (; j < N; ++j) {
        double const r2 = (std::max)(r2_[j], 0.0);
        nearest[0] = (std::min)(nearest[0], r2);
        sum[0] += coefficients[j] * r2 * sqrt(r2);
      }
      s[k] = v + ((sum[0] + sum[1]) + (sum[2] + sum[3]));
      delta[k] = sqrt((std::min)((std::min)(nearest[0], nearest[1]), (std::min)(nearest[2], nearest[3])));
    }
    auto const [s_min, s_max] = std::minmax_element(s.begin(), s.end());
    auto const [d_min, d_max] = std::minmax_element(delta.begin(), delta.end());
    double const w = weights_[iteration_++ % 4];
    size_t chosen = 0;
    double best_score = std::numeric_limits<double>::infinity();
    for (size_t k = 0; k < K; ++k) {
      double const vr = *s_max > *s_min ? (s[k] - *s_min) / (*s_max - *s_min) : 1;
      double const vd = *d_max > *d_min ? (*d_max - delta[k]) / (*d_max - *d_min) : 1;
      double const score = w * vr + (1 - w) * vd;
      if (score < best_score) {
        best_score = score;
        chosen = k;
      }
    }

    std::vector<double> u(xb, xb + d_);
    for (size_t e = offset_[chosen]; e < offset_[chosen + 1]; ++e) {
      u[index_[e]] = value_[e];
    }
    ResultType const previous = best_value_;
    X_.insert(X_.end(), u.begin(), u.end());
    f_.push_back(evaluate(to_argument(u.data())));
    add_point(N);

    // Step 2(f), where a success is an improvement on the best value by more than 1e-3 of its magnitude:
    bool const success = isfinite(previous) ? best_value_ < previous - abs(previous) / 1000 : best_value_ < previous;
    if (success) {
      ++successes_;
      failures_ = 0;
    }
    else {
      successes_ = 0;
      ++failures_;
    }
    if (successes_ >= success_tolerance_) {
      sigma_ = (std::min)(2 * sigma_, sigma_max_);
      successes_ = 0;
    }
    if (failures_ >= failure_tolerance_) {
      sigma_ = (std::max)(sigma_ / 2, sigma_min_);
      failures_ = 0;
    }
  }

  Func const &cost_function_;
  dycors_parameters<ArgumentContainer> const &params_;
  URBG &gen_;
  ResultType target_value_;
  std::atomic<bool> *cancellation_;
  std::atomic<ResultType> *current_minimum_cost_;
  std::vector<std::pair<ArgumentContainer, ResultType>> *queries_;
  size_t d_;
  size_t max_function_calls_;
  size_t initial_points_;
  size_t candidates_;
  size_t failure_tolerance_;
  size_t calls_ = 0;
  ArgumentContainer best_argument_;
  ResultType best_value_ = std::numeric_limits<ResultType>::infinity();
  // The evaluated points (unit cube, row-major) and their values, and the index of the best:
  std::vector<double> X_;
  std::vector<double> f_;
  size_t best_ = 0;
  // The step size and its counters, and the iteration number, which picks the weight:
  double sigma_ = sigma_init_;
  size_t successes_ = 0;
  size_t failures_ = 0;
  size_t iteration_ = 0;
  // The calls of the initial guess and design:
  size_t n0_ = 0;
  // The surrogate: the points T and the inverse of the matrix of their rows (1, t_j), l(x_k), q(x_k) and r(x_k) for
  // each point (row-major, N x (d + 1)), and the Cholesky factor of K on the other points (nodes_), by rows.
  std::vector<size_t> T_;
  std::vector<double> Pinv_, PhiT_;
  std::vector<double> lagrange_, q_, r_;
  std::vector<size_t> nodes_;
  // Scratch for the candidates (see step()), and the points by coordinates, with room for every call:
  std::vector<size_t> index_, offset_;
  std::vector<double> value_, D_, r2_, Xt_;
  size_t stride_ = 0;
  size_t Xt_filled_ = 0;
  std::vector<std::vector<double>> L_;
};

} // namespace detail

BOOST_MATH_EXPORT template <typename ArgumentContainer, class Func, class URBG>
ArgumentContainer dycors(
    const Func cost_function, dycors_parameters<ArgumentContainer> const &params, URBG &gen,
    std::invoke_result_t<Func, ArgumentContainer> target_value =
        std::numeric_limits<std::invoke_result_t<Func, ArgumentContainer>>::quiet_NaN(),
    std::atomic<bool> *cancellation = nullptr,
    std::atomic<std::invoke_result_t<Func, ArgumentContainer>> *current_minimum_cost = nullptr,
    std::vector<std::pair<ArgumentContainer, std::invoke_result_t<Func, ArgumentContainer>>> *queries = nullptr) {
  detail::dycors_engine<ArgumentContainer, Func, URBG> engine(cost_function, params, gen, target_value, cancellation,
                                                              current_minimum_cost, queries);
  return engine.optimize();
}

} BOOST_MATH_NAMESPACE_END // namespace boost::math::optimization
#endif
