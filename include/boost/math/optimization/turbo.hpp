/*
 * Copyright Nick Thompson, 2026
 * Use, modification and distribution are subject to the
 * Boost Software License, Version 1.0. (See accompanying file
 * LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef BOOST_MATH_OPTIMIZATION_TURBO_HPP
#define BOOST_MATH_OPTIMIZATION_TURBO_HPP

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

// TuRBO-1 of: Eriksson, D., Pearce, M., Gardner, J., Turner, R. D., Poloczek, M. (2019).
// Scalable global optimization via local Bayesian optimization. NeurIPS 32. arXiv:1910.01739.
// The defaults follow the authors' code (github.com/uber-research/TuRBO, turbo_1.py and gp.py), "the reference".
// Section numbers in the comments below refer to the paper and its supplement.

BOOST_MATH_EXPORT template <typename ArgumentContainer> struct turbo_parameters {
  using Real = typename ArgumentContainer::value_type;
  ArgumentContainer lower_bounds;
  ArgumentContainer upper_bounds;
  // Replaces the first point of the first initial design.
  ArgumentContainer const *initial_guess = nullptr;
  // 0 means 100 times the dimension.
  size_t max_function_calls = 0;
  // The points of the Latin hypercube design that starts each trust region; 0 means twice the dimension.
  size_t initial_points = 0;
  // q, the number of points proposed (by independent Thompson samples) per Gaussian process fit.
  size_t batch_size = 1;
  // The candidate points of each Thompson sample; 0 means min(100n, 1000). The paper (Section D) has min(100n, 5000),
  // but the joint sample costs a Cholesky factorization of the candidates' covariance, O(candidates^3).
  size_t candidates = 0;
  // The Gaussian process is fitted to at most this many of the trust region's points, those nearest its centre.
  size_t max_training_points = 200;
  // Adam steps on the log marginal likelihood per fit, as in the reference.
  size_t hyperparameter_steps = 50;
};

BOOST_MATH_EXPORT template <typename ArgumentContainer>
void validate_turbo_parameters(turbo_parameters<ArgumentContainer> const &params) {
  std::ostringstream oss;
  detail::validate_bounds(params.lower_bounds, params.upper_bounds);
  for (size_t i = 0; i < params.lower_bounds.size(); ++i) {
    if (!(params.lower_bounds[i] < params.upper_bounds[i])) {
      oss << __FILE__ << ":" << __LINE__ << ":" << __func__;
      oss << ": TuRBO requires lower_bounds[i] < upper_bounds[i], but at i = " << i << " the bounds are "
          << params.lower_bounds[i] << " and " << params.upper_bounds[i] << ".";
      throw std::domain_error(oss.str());
    }
  }
  if (params.initial_guess) {
    detail::validate_initial_guess(*params.initial_guess, params.lower_bounds, params.upper_bounds);
  }
  if (params.batch_size == 0) {
    oss << __FILE__ << ":" << __LINE__ << ":" << __func__;
    oss << ": The batch size must be at least 1.";
    throw std::invalid_argument(oss.str());
  }
  if (params.candidates != 0 && params.candidates < params.batch_size) {
    oss << __FILE__ << ":" << __LINE__ << ":" << __func__;
    oss << ": There must be at least as many candidates as the batch size, but candidates = " << params.candidates
        << " and batch_size = " << params.batch_size << ".";
    throw std::invalid_argument(oss.str());
  }
  if (params.max_training_points == 0) {
    oss << __FILE__ << ":" << __LINE__ << ":" << __func__;
    oss << ": The Gaussian process needs at least one training point.";
    throw std::invalid_argument(oss.str());
  }
}

namespace detail {

// Hand-rolled dense linear algebra on row-major double arrays. The kernels with long inner loops are written so that
// the compiler can vectorize them without reassociating floating point sums: either as axpy updates of contiguous
// rows, or as several rows' dot products with one shared operand, which gives independent accumulators.

// s[r] = sum_{k < m} a[r][k] b[k] for four rows a[0..3].
inline void turbo_dot4(double const *const *a, double const *b, size_t m, double *s) {
  double s0 = 0, s1 = 0, s2 = 0, s3 = 0, t0 = 0, t1 = 0, t2 = 0, t3 = 0;
  size_t k = 0;
  for (; k + 1 < m; k += 2) {
    double const x = b[k];
    double const y = b[k + 1];
    s0 += a[0][k] * x;
    t0 += a[0][k + 1] * y;
    s1 += a[1][k] * x;
    t1 += a[1][k + 1] * y;
    s2 += a[2][k] * x;
    t2 += a[2][k + 1] * y;
    s3 += a[3][k] * x;
    t3 += a[3][k + 1] * y;
  }
  if (k < m) {
    double const x = b[k];
    s0 += a[0][k] * x;
    s1 += a[1][k] * x;
    s2 += a[2][k] * x;
    s3 += a[3][k] * x;
  }
  s[0] = s0 + t0;
  s[1] = s1 + t1;
  s[2] = s2 + t2;
  s[3] = s3 + t3;
}

inline double turbo_dot(double const *a, double const *b, size_t m) {
  double s0 = 0, s1 = 0, s2 = 0, s3 = 0;
  size_t k = 0;
  for (; k + 3 < m; k += 4) {
    s0 += a[k] * b[k];
    s1 += a[k + 1] * b[k + 1];
    s2 += a[k + 2] * b[k + 2];
    s3 += a[k + 3] * b[k + 3];
  }
  for (; k < m; ++k) {
    s0 += a[k] * b[k];
  }
  return (s0 + s1) + (s2 + s3);
}

// C[i][j] -= sum_{k < len} X[i][k] X[j][k] for lo <= j <= i < hi, where row i of X starts at X + i ldx and C has
// leading dimension ldc. Tiles of 64 rows j stay in cache while every block of four rows i below them goes by.
inline void turbo_subtract_gram(double *C, size_t ldc, double const *X, size_t ldx, size_t len, size_t lo, size_t hi) {
  constexpr size_t tile = 64;
  double s[4];
  double const *rows[4];
  for (size_t j0 = lo; j0 < hi; j0 += tile) {
    size_t const j1 = (std::min)(hi, j0 + tile);
    for (size_t i = j0; i < hi; i += 4) {
      size_t const nb = (std::min)(size_t(4), hi - i);
      for (size_t b = 0; b < 4; ++b) {
        rows[b] = X + (i + (b < nb ? b : 0)) * ldx;
      }
      size_t const jend = (std::min)(j1, i + nb);
      for (size_t j = j0; j < jend; ++j) {
        double const *xj = X + j * ldx;
        if (j < i) {
          turbo_dot4(rows, xj, len, s);
          for (size_t b = 0; b < nb; ++b) {
            C[(i + b) * ldc + j] -= s[b];
          }
        }
        else {
          for (size_t b = j - i; b < nb; ++b) {
            C[(i + b) * ldc + j] -= turbo_dot(rows[b], xj, len);
          }
        }
      }
    }
  }
}

// The Cholesky factorization A = L L^T in place: on input the lower triangle of the n x n matrix A (leading dimension
// n), on output that of L. The upper triangle is not referenced. Returns false if A is not numerically positive
// definite. Right-looking by blocks of 64 columns; within a block, left-looking four rows at a time.
inline bool turbo_cholesky(double *A, size_t n) {
  using std::sqrt;
  constexpr size_t block = 64;
  double s[4];
  double *rows[4];
  double const *parts[4];
  for (size_t k0 = 0; k0 < n; k0 += block) {
    size_t const k1 = (std::min)(n, k0 + block);
    // Columns k0 to k1 of L, from the diagonal block down:
    for (size_t i = k0; i < n; i += 4) {
      size_t const nb = (std::min)(size_t(4), n - i);
      for (size_t b = 0; b < 4; ++b) {
        rows[b] = A + (i + (b < nb ? b : 0)) * n;
        parts[b] = rows[b] + k0;
      }
      size_t const jend = (std::min)(k1, i + nb);
      for (size_t j = k0; j < jend; ++j) {
        double const *lj = A + j * n;
        if (j < i) {
          turbo_dot4(parts, lj + k0, j - k0, s);
          for (size_t b = 0; b < nb; ++b) {
            rows[b][j] = (rows[b][j] - s[b]) / lj[j];
          }
          continue;
        }
        // Column j meets these rows on the diagonal: row j first, then the rows below it.
        for (size_t b = j - i; b < nb; ++b) {
          double const d = rows[b][j] - turbo_dot(parts[b], lj + k0, j - k0);
          if (b == j - i) {
            if (!(d > 0)) {
              return false;
            }
            rows[b][j] = sqrt(d);
          }
          else {
            rows[b][j] = d / lj[j];
          }
        }
      }
    }
    // The trailing submatrix loses the contribution of these columns:
    turbo_subtract_gram(A, n, A + k0, n, k1 - k0, k1, n);
  }
  return true;
}

// Solves L X = B in place for the n x m matrix B (leading dimension m), where L is the lower triangle of an n x n
// matrix (leading dimension n).
inline void turbo_forward_solve(double const *L, size_t n, double *B, size_t m) {
  for (size_t i = 0; i < n; ++i) {
    double *bi = B + i * m;
    for (size_t k = 0; k < i; ++k) {
      double const lik = L[i * n + k];
      double const *bk = B + k * m;
      for (size_t c = 0; c < m; ++c) {
        bi[c] -= lik * bk[c];
      }
    }
    double const inv = 1 / L[i * n + i];
    for (size_t c = 0; c < m; ++c) {
      bi[c] *= inv;
    }
  }
}

// Solves L L^T x = b in place.
inline void turbo_cholesky_solve(double const *L, size_t n, double *b) {
  for (size_t i = 0; i < n; ++i) {
    b[i] = (b[i] - turbo_dot(L + i * n, b, i)) / L[i * n + i];
  }
  // L^T x = b, by rows of L: once x_i is known, remove its contribution L_ik x_i from b_k for k < i.
  for (size_t i = n; i-- > 0;) {
    double const *li = L + i * n;
    double const xi = b[i] / li[i];
    b[i] = xi;
    for (size_t k = 0; k < i; ++k) {
      b[k] -= li[k] * xi;
    }
  }
}

// The lower triangle of (L L^T)^{-1} = Z^T Z, where Z = L^{-1}. Z is n x n scratch.
inline void turbo_cholesky_inverse(double const *L, size_t n, double *Z, double *Kinv) {
  for (size_t i = 0; i < n; ++i) {
    double *zi = Z + i * n;
    for (size_t j = 0; j < i; ++j) {
      zi[j] = 0;
    }
    for (size_t k = 0; k < i; ++k) {
      double const lik = L[i * n + k];
      double const *zk = Z + k * n;
      for (size_t j = 0; j <= k; ++j) {
        zi[j] -= lik * zk[j];
      }
    }
    double const inv = 1 / L[i * n + i];
    for (size_t j = 0; j < i; ++j) {
      zi[j] *= inv;
    }
    zi[i] = inv;
  }
  for (size_t a = 0; a < n; ++a) {
    for (size_t b = 0; b <= a; ++b) {
      Kinv[a * n + b] = 0;
    }
  }
  for (size_t k = 0; k < n; ++k) {
    double const *zk = Z + k * n;
    for (size_t a = 0; a <= k; ++a) {
      double const zka = zk[a];
      double *row = Kinv + a * n;
      for (size_t b = 0; b <= a; ++b) {
        row[b] += zka * zk[b];
      }
    }
  }
}

// The ARD Matern-5/2 kernel (Section C) on the unit cube: s (1 + sqrt(5) r + 5 r^2/3) exp(-sqrt(5) r), where
// r^2 = sum_i (x_i - y_i)^2 / l_i^2. inv_l2[i] = 1/l_i^2.
inline double turbo_matern52(double const *x, double const *y, double const *inv_l2, size_t d) {
  using std::exp;
  using std::sqrt;
  double r2 = 0;
  for (size_t i = 0; i < d; ++i) {
    double const t = x[i] - y[i];
    r2 += t * t * inv_l2[i];
  }
  double const r = sqrt(5 * r2);
  return (1 + r + r * r / 3) * exp(-r);
}

// The state of a TuRBO-1 run; turbo() below drives it.
template <typename ArgumentContainer, class Func, class URBG> class turbo_engine {
public:
  using Real = typename ArgumentContainer::value_type;
  using ResultType = std::invoke_result_t<Func, ArgumentContainer>;

  turbo_engine(Func const &cost_function, turbo_parameters<ArgumentContainer> const &params, URBG &gen,
               ResultType target_value, std::atomic<bool> *cancellation, std::atomic<ResultType> *current_minimum_cost,
               std::vector<std::pair<ArgumentContainer, ResultType>> *queries)
      : cost_function_(cost_function), params_(params), gen_(gen), target_value_(target_value),
        cancellation_(cancellation), current_minimum_cost_(current_minimum_cost), queries_(queries) {
    validate_turbo_parameters(params);
    n_ = params.lower_bounds.size();
    q_ = params.batch_size;
    max_function_calls_ = params.max_function_calls ? params.max_function_calls : 100 * n_;
    initial_points_ = params.initial_points ? params.initial_points : 2 * n_;
    candidates_ = params.candidates ? params.candidates : (std::min)(100 * n_, size_t(1000));
    candidates_ = (std::max)(candidates_, q_);
    // tau_fail = ceil(max(4/q, n/q)), as in the reference; the paper (Section D) has ceil(n/q).
    failure_tolerance_ = ((std::max)(size_t(4), n_) + q_ - 1) / q_;
    if constexpr (detail::has_resize_v<ArgumentContainer>) {
      best_argument_.resize(n_);
    }
    for (size_t i = 0; i < n_; ++i) {
      best_argument_[i] = params.initial_guess ? (*params.initial_guess)[i]
                                               : params.lower_bounds[i] + (params.upper_bounds[i] - params.lower_bounds[i]) / 2;
    }
  }

  ArgumentContainer optimize() {
    bool first = true;
    while (!done()) {
      restart(first);
      first = false;
      while (!done() && length_ >= length_min_) {
        step();
      }
    }
    return best_argument_;
  }

private:
  // Section D and the reference:
  static constexpr double length_init_ = 0.8;
  static constexpr double length_min_ = 0.0078125; // 2^-7
  static constexpr double length_max_ = 1.6;
  static constexpr size_t success_tolerance_ = 3;

  // The hyperparameter bounds of the reference (gp.py); the paper (Section C) has 0.1 for the upper noise bound.
  static constexpr double lengthscale_bounds_[2] = {0.005, 2.0};
  static constexpr double outputscale_bounds_[2] = {0.05, 20.0};
  static constexpr double noise_bounds_[2] = {0.0005, 0.2};

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
      x.resize(n_);
    }
    for (size_t i = 0; i < n_; ++i) {
      Real const lo = params_.lower_bounds[i];
      Real const hi = params_.upper_bounds[i];
      Real const xi = lo + static_cast<Real>(u[i]) * (hi - lo);
      x[i] = xi < lo ? lo : (xi > hi ? hi : xi);
    }
    return x;
  }

  // Calls the cost function at x (with unit-cube coordinates u), and appends the point to the trust region's data.
  // NaN counts as +infinity.
  void evaluate(ArgumentContainer const &x, double const *u) {
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
      if (current_minimum_cost_) {
        *current_minimum_cost_ = f;
      }
    }
    X_.insert(X_.end(), u, u + n_);
    y_.push_back(static_cast<double>(f));
  }

  double uniform01() { return std::uniform_real_distribution<double>(0, 1)(gen_); }

  // A random Latin hypercube design of m points in the box [lo, hi] (row-major m x n), as latin_hypercube in the
  // reference: one point in each of m equal slices of every coordinate, uniformly distributed within its slice.
  std::vector<double> latin_hypercube(size_t m, double const *lo, double const *hi) {
    std::vector<double> design(m * n_);
    std::vector<size_t> slices(m);
    for (size_t i = 0; i < n_; ++i) {
      for (size_t k = 0; k < m; ++k) {
        slices[k] = k;
      }
      std::shuffle(slices.begin(), slices.end(), gen_);
      for (size_t k = 0; k < m; ++k) {
        double const t = (static_cast<double>(slices[k]) + uniform01()) / static_cast<double>(m);
        design[k * n_ + i] = (std::min)(lo[i] + t * (hi[i] - lo[i]), hi[i]);
      }
    }
    return design;
  }

  // Starts a new trust region from a new initial design.
  void restart(bool first) {
    length_ = length_init_;
    successes_ = 0;
    failures_ = 0;
    X_.clear();
    y_.clear();
    std::vector<double> lo(n_, 0.0), hi(n_, 1.0);
    std::vector<double> design = latin_hypercube(initial_points_, lo.data(), hi.data());
    for (size_t k = 0; k < initial_points_ && !done(); ++k) {
      double *u = design.data() + k * n_;
      if (first && k == 0 && params_.initial_guess) {
        auto const &g = *params_.initial_guess;
        for (size_t i = 0; i < n_; ++i) {
          Real const lo_i = params_.lower_bounds[i];
          Real const hi_i = params_.upper_bounds[i];
          u[i] = static_cast<double>((g[i] - lo_i) / (hi_i - lo_i));
        }
        evaluate(g, u);
      }
      else {
        evaluate(to_argument(u), u);
      }
    }
  }

  // One iteration of TuRBO-1: fit the Gaussian process, Thompson sample a batch in the trust region, evaluate it, and
  // resize the trust region.
  void step() {
    using std::abs;
    using std::isfinite;
    size_t const N_all = y_.size();
    size_t center = 0;
    for (size_t k = 1; k < N_all; ++k) {
      if (y_[k] < y_[center]) {
        center = k;
      }
    }
    double const *xc = X_.data() + center * n_;
    select_training_set(xc);
    bool const model = fit();

    // The trust region (Section 2): side lengths L l_i / (prod_j l_j)^(1/n), centred at the best point, clipped to the
    // unit cube. The reference first divides by the mean lengthscale, for stability.
    std::vector<double> weights(n_, 1.0);
    if (model) {
      double mean = 0;
      for (size_t i = 0; i < n_; ++i) {
        mean += lengthscale_[i];
      }
      mean /= static_cast<double>(n_);
      double log_prod = 0;
      for (size_t i = 0; i < n_; ++i) {
        weights[i] = lengthscale_[i] / mean;
        log_prod += std::log(weights[i]);
      }
      double const scale = std::exp(-log_prod / static_cast<double>(n_));
      for (auto &w : weights) {
        w *= scale;
      }
    }
    std::vector<double> lo(n_), hi(n_);
    for (size_t i = 0; i < n_; ++i) {
      lo[i] = (std::max)(xc[i] - weights[i] * length_ / 2, 0.0);
      hi[i] = (std::min)(xc[i] + weights[i] * length_ / 2, 1.0);
    }

    // The candidates (Section D): a Latin hypercube in the trust region, then each coordinate is kept at the centre
    // unless it is perturbed, with probability min(20/n, 1); a candidate with no perturbed coordinate gets one.
    size_t const M = candidates_;
    std::vector<double> C = latin_hypercube(M, lo.data(), hi.data());
    double const p_perturb = (std::min)(20.0 / static_cast<double>(n_), 1.0);
    if (p_perturb < 1) {
      std::uniform_int_distribution<size_t> pick(0, n_ - 1);
      for (size_t k = 0; k < M; ++k) {
        double *c = C.data() + k * n_;
        bool any = false;
        for (size_t i = 0; i < n_; ++i) {
          if (uniform01() <= p_perturb) {
            any = true;
          }
          else {
            c[i] = xc[i];
          }
        }
        if (!any) {
          size_t const i = pick(gen_);
          c[i] = lo[i] + uniform01() * (hi[i] - lo[i]);
        }
      }
    }

    // Thompson sampling: q joint samples of the posterior (with the noise, as the reference samples the likelihood)
    // over the candidates, and the argmin of each, never the same candidate twice.
    std::vector<size_t> chosen;
    if (model) {
      thompson_sample(C, M, chosen);
    }
    else {
      // No finite function value yet: no model, so take random candidates.
      std::vector<size_t> order(M);
      for (size_t k = 0; k < M; ++k) {
        order[k] = k;
      }
      std::shuffle(order.begin(), order.end(), gen_);
      chosen.assign(order.begin(), order.begin() + static_cast<std::ptrdiff_t>(q_));
    }

    double const local_best = y_[center];
    double batch_best = std::numeric_limits<double>::infinity();
    for (size_t k : chosen) {
      if (done()) {
        return;
      }
      double const *u = C.data() + k * n_;
      evaluate(to_argument(u), u);
      batch_best = (std::min)(batch_best, y_.back());
    }
    // A success is an improvement on the trust region's best value by more than 1e-3 of its magnitude (reference).
    bool const success = isfinite(local_best) ? batch_best < local_best - 1e-3 * abs(local_best) : batch_best < local_best;
    if (success) {
      ++successes_;
      failures_ = 0;
    }
    else {
      successes_ = 0;
      ++failures_;
    }
    if (successes_ == success_tolerance_) {
      length_ = (std::min)(2 * length_, length_max_);
      successes_ = 0;
    }
    else if (failures_ == failure_tolerance_) {
      length_ /= 2;
      failures_ = 0;
    }
  }

  // The training set: the max_training_points points of the trust region nearest its centre xc (all of them if
  // there are fewer), in the unit cube, with the function values standardized as in the reference (median and
  // standard deviation). Non-finite values are replaced by the largest finite one.
  void select_training_set(double const *xc) {
    size_t const N_all = y_.size();
    std::vector<size_t> idx(N_all);
    for (size_t k = 0; k < N_all; ++k) {
      idx[k] = k;
    }
    size_t const N = (std::min)(N_all, params_.max_training_points);
    if (N < N_all) {
      std::vector<double> dist(N_all);
      for (size_t k = 0; k < N_all; ++k) {
        double s = 0;
        for (size_t i = 0; i < n_; ++i) {
          double const t = X_[k * n_ + i] - xc[i];
          s += t * t;
        }
        dist[k] = s;
      }
      std::nth_element(idx.begin(), idx.begin() + static_cast<std::ptrdiff_t>(N), idx.end(),
                       [&](size_t a, size_t b) { return dist[a] < dist[b] || (dist[a] == dist[b] && a < b); });
      idx.resize(N);
      std::sort(idx.begin(), idx.end());
    }
    N_ = N;
    Xt_.resize(N * n_);
    yt_.resize(N);
    using std::isfinite;
    double worst = -std::numeric_limits<double>::infinity();
    for (size_t k = 0; k < N; ++k) {
      std::copy(X_.begin() + static_cast<std::ptrdiff_t>(idx[k] * n_),
                X_.begin() + static_cast<std::ptrdiff_t>((idx[k] + 1) * n_), Xt_.begin() + static_cast<std::ptrdiff_t>(k * n_));
      yt_[k] = y_[idx[k]];
      if (isfinite(yt_[k])) {
        worst = (std::max)(worst, yt_[k]);
      }
    }
    finite_ = isfinite(worst);
    if (!finite_) {
      return;
    }
    for (auto &v : yt_) {
      if (!isfinite(v)) {
        v = worst;
      }
    }
    std::vector<double> sorted = yt_;
    std::nth_element(sorted.begin(), sorted.begin() + static_cast<std::ptrdiff_t>(N / 2), sorted.end());
    double median = sorted[N / 2];
    if (N % 2 == 0) {
      median = (median + *std::max_element(sorted.begin(), sorted.begin() + static_cast<std::ptrdiff_t>(N / 2))) / 2;
    }
    double mean = 0;
    for (double v : yt_) {
      mean += v;
    }
    mean /= static_cast<double>(N);
    double var = 0;
    for (double v : yt_) {
      var += (v - mean) * (v - mean);
    }
    using std::sqrt;
    double sigma = sqrt(var / static_cast<double>(N));
    if (!(sigma >= 1e-6)) {
      sigma = 1;
    }
    for (auto &v : yt_) {
      v = (v - median) / sigma;
    }
  }

  static double sigmoid(double t) { return 1 / (1 + std::exp(-t)); }
  static double logit(double p) { return std::log(p / (1 - p)); }
  static double from_raw(double raw, double const *bounds) { return bounds[0] + (bounds[1] - bounds[0]) * sigmoid(raw); }
  static double to_raw(double value, double const *bounds) { return logit((value - bounds[0]) / (bounds[1] - bounds[0])); }

  // Sets the hyperparameters from theta = (raw lengthscales, raw outputscale, raw noise, mean), builds the kernel
  // matrix K + noise I in L_, and factors it. pair_m_ and pair_dm_ receive the kernel's radial factor and the factor
  // of its lengthscale derivatives at each pair of training points. Returns false if the factorization fails.
  bool factor(std::vector<double> const &theta) {
    using std::exp;
    using std::sqrt;
    size_t const N = N_;
    for (size_t i = 0; i < n_; ++i) {
      lengthscale_[i] = from_raw(theta[i], lengthscale_bounds_);
      inv_l2_[i] = 1 / (lengthscale_[i] * lengthscale_[i]);
    }
    outputscale_ = from_raw(theta[n_], outputscale_bounds_);
    noise_ = from_raw(theta[n_ + 1], noise_bounds_);
    mean_ = theta[n_ + 2];
    L_.resize(N * N);
    pair_m_.resize(N * (N - 1) / 2);
    pair_dm_.resize(pair_m_.size());
    double const root5 = sqrt(5.0);
    size_t p = 0;
    for (size_t a = 0; a < N; ++a) {
      double *row = L_.data() + a * N;
      for (size_t b = 0; b < a; ++b, ++p) {
        double const *D = sqdiff_.data() + p * n_;
        double r2 = 0;
        for (size_t i = 0; i < n_; ++i) {
          r2 += D[i] * inv_l2_[i];
        }
        double const r = root5 * sqrt(r2);
        double const e = exp(-r);
        pair_m_[p] = (1 + r + r * r / 3) * e;
        // d k / d l_i = s (5/3) (1 + sqrt(5) rho) exp(-sqrt(5) rho) (x_i - y_i)^2 / l_i^3, with r = sqrt(5) rho:
        pair_dm_[p] = outputscale_ * (5.0 / 3) * (1 + r) * e;
        row[b] = outputscale_ * pair_m_[p];
      }
      row[a] = outputscale_ + noise_;
    }
    return turbo_cholesky(L_.data(), N);
  }

  // The loss, -log p(y | theta)/N (as GPyTorch's ExactMarginalLogLikelihood, up to the constant log(2 pi)/2), and
  // its gradient with respect to theta. Leaves L_ and alpha_ = (K + noise I)^{-1} (y - mean). Returns false if the
  // kernel matrix cannot be factored.
  bool loss_and_gradient(std::vector<double> const &theta, double &loss, std::vector<double> &grad) {
    using std::log;
    size_t const N = N_;
    if (!factor(theta)) {
      return false;
    }
    alpha_.resize(N);
    for (size_t a = 0; a < N; ++a) {
      alpha_[a] = yt_[a] - mean_;
    }
    turbo_cholesky_solve(L_.data(), N, alpha_.data());
    loss = 0;
    for (size_t a = 0; a < N; ++a) {
      loss += (yt_[a] - mean_) * alpha_[a] / 2 + log(L_[a * N + a]);
    }
    loss /= static_cast<double>(N);
    Z_.resize(N * N);
    Kinv_.resize(N * N);
    turbo_cholesky_inverse(L_.data(), N, Z_.data(), Kinv_.data());
    // d loss / d theta = -tr(W dK/d theta)/(2N), where W = alpha alpha^T - K^{-1}:
    double trace_noise = 0;
    double trace_scale = 0;
    std::vector<double> trace_l(n_, 0.0);
    size_t p = 0;
    for (size_t a = 0; a < N; ++a) {
      double const *kinv = Kinv_.data() + a * N;
      for (size_t b = 0; b < a; ++b, ++p) {
        double const w = 2 * (alpha_[a] * alpha_[b] - kinv[b]);
        trace_scale += w * pair_m_[p];
        double const wd = w * pair_dm_[p];
        double const *D = sqdiff_.data() + p * n_;
        for (size_t i = 0; i < n_; ++i) {
          trace_l[i] += wd * D[i];
        }
      }
      double const w = alpha_[a] * alpha_[a] - kinv[a];
      trace_noise += w;
      trace_scale += w;
    }
    double sum_alpha = 0;
    for (size_t a = 0; a < N; ++a) {
      sum_alpha += alpha_[a];
    }
    double const scale = -1 / (2 * static_cast<double>(N));
    // The chain rule through value = lo + (hi - lo) sigmoid(raw):
    auto dvalue = [](double raw, double const *bounds) {
      double const sg = sigmoid(raw);
      return (bounds[1] - bounds[0]) * sg * (1 - sg);
    };
    for (size_t i = 0; i < n_; ++i) {
      double const l = lengthscale_[i];
      grad[i] = scale * trace_l[i] / (l * l * l) * dvalue(theta[i], lengthscale_bounds_);
    }
    grad[n_] = scale * trace_scale * dvalue(theta[n_], outputscale_bounds_);
    grad[n_ + 1] = scale * trace_noise * dvalue(theta[n_ + 1], noise_bounds_);
    grad[n_ + 2] = -sum_alpha / static_cast<double>(N);
    return true;
  }

  // Fits the Gaussian process (Section C) to the training set by Adam on the loss above, over the hyperparameters
  // transformed to the real line as in GPyTorch: value = lo + (hi - lo) sigmoid(raw). Starts from the reference's
  // values every time (lengthscales 0.5, outputscale 1, noise 0.005, mean 0). Leaves L_ and alpha_ for the
  // posterior. Returns false if there is no model.
  bool fit() {
    using std::isfinite;
    using std::sqrt;
    if (!finite_) {
      return false;
    }
    size_t const N = N_;
    sqdiff_.resize(N * (N - 1) / 2 * n_);
    size_t p = 0;
    for (size_t a = 0; a < N; ++a) {
      for (size_t b = 0; b < a; ++b, ++p) {
        for (size_t i = 0; i < n_; ++i) {
          double const t = Xt_[a * n_ + i] - Xt_[b * n_ + i];
          sqdiff_[p * n_ + i] = t * t;
        }
      }
    }
    lengthscale_.resize(n_);
    inv_l2_.resize(n_);
    size_t const dim = n_ + 3;
    std::vector<double> theta(dim);
    for (size_t i = 0; i < n_; ++i) {
      theta[i] = to_raw(0.5, lengthscale_bounds_);
    }
    theta[n_] = to_raw(1.0, outputscale_bounds_);
    theta[n_ + 1] = to_raw(0.005, noise_bounds_);
    theta[n_ + 2] = 0;
    std::vector<double> m1(dim, 0.0), m2(dim, 0.0), grad(dim), last_good = theta;
    // Adam with the PyTorch defaults and the reference's learning rate 0.1:
    double const lr = 0.1, beta1 = 0.9, beta2 = 0.999, eps = 1e-8;
    double beta1_t = 1, beta2_t = 1;
    double loss;
    // Set if Adam stops on hyperparameters whose kernel matrix cannot be factored, or whose loss or gradient is not
    // finite; the model then uses the last hyperparameters that were fine.
    bool bad = false;
    for (size_t step = 0; step < params_.hyperparameter_steps; ++step) {
      bool finite = loss_and_gradient(theta, loss, grad) && isfinite(loss);
      for (size_t k = 0; k < dim && finite; ++k) {
        finite = isfinite(grad[k]);
      }
      if (!finite) {
        bad = true;
        break;
      }
      last_good = theta;
      beta1_t *= beta1;
      beta2_t *= beta2;
      for (size_t k = 0; k < dim; ++k) {
        m1[k] = beta1 * m1[k] + (1 - beta1) * grad[k];
        m2[k] = beta2 * m2[k] + (1 - beta2) * grad[k] * grad[k];
        theta[k] -= lr * (m1[k] / (1 - beta1_t)) / (sqrt(m2[k] / (1 - beta2_t)) + eps);
      }
    }
    // The model at the final (or last good) hyperparameters:
    if (bad || !factor(theta)) {
      if (!factor(last_good)) {
        return false;
      }
    }
    alpha_.resize(N);
    for (size_t a = 0; a < N; ++a) {
      alpha_[a] = yt_[a] - mean_;
    }
    turbo_cholesky_solve(L_.data(), N, alpha_.data());
    return true;
  }

  // q joint posterior samples over the M candidates C (row-major M x n), and the argmin of each.
  void thompson_sample(std::vector<double> const &C, size_t M, std::vector<size_t> &chosen) {
    using std::sqrt;
    size_t const N = N_;
    // Kxc (N x M), the posterior mean, then V = L^{-1} Kxc in place:
    V_.resize(N * M);
    std::vector<double> mu(M, mean_);
    for (size_t a = 0; a < N; ++a) {
      double const *xa = Xt_.data() + a * n_;
      double *row = V_.data() + a * M;
      for (size_t c = 0; c < M; ++c) {
        row[c] = outputscale_ * turbo_matern52(xa, C.data() + c * n_, inv_l2_.data(), n_);
        mu[c] += row[c] * alpha_[a];
      }
    }
    turbo_forward_solve(L_.data(), N, V_.data(), M);
    // The transpose, M x N, so that the rows of the Gram matrix V^T V are contiguous dot products:
    Vt_.resize(M * N);
    for (size_t a = 0; a < N; ++a) {
      for (size_t c = 0; c < M; ++c) {
        Vt_[c * N + a] = V_[a * M + c];
      }
    }
    // The posterior covariance Kcc + noise I - V^T V, then its Cholesky factor, with jitter if it is needed:
    S_.resize(M * M);
    for (size_t a = 0; a < M; ++a) {
      double const *xa = C.data() + a * n_;
      double *row = S_.data() + a * M;
      for (size_t b = 0; b < a; ++b) {
        row[b] = outputscale_ * turbo_matern52(xa, C.data() + b * n_, inv_l2_.data(), n_);
      }
      row[a] = outputscale_ + noise_;
    }
    turbo_subtract_gram(S_.data(), M, Vt_.data(), N, N, 0, M);
    // The noise bounds the covariance below, so jitter should never be needed; if it is, retry from a copy.
    covariance_ = S_;
    bool factored = turbo_cholesky(S_.data(), M);
    for (double jitter = 1e-8; !factored && jitter < 1; jitter *= 10) {
      S_ = covariance_;
      for (size_t a = 0; a < M; ++a) {
        S_[a * M + a] += jitter;
      }
      factored = turbo_cholesky(S_.data(), M);
    }
    std::normal_distribution<double> normal(0, 1);
    std::vector<double> z(M);
    std::vector<bool> taken(M, false);
    for (size_t j = 0; j < q_; ++j) {
      for (auto &zi : z) {
        zi = normal(gen_);
      }
      size_t best = M;
      double best_value = std::numeric_limits<double>::infinity();
      for (size_t c = 0; c < M; ++c) {
        // Without a factor (it should not happen, as the noise bounds the covariance below), the marginals:
        double const f = factored ? mu[c] + turbo_dot(S_.data() + c * M, z.data(), c + 1)
                                  : mu[c] + sqrt(outputscale_ + noise_) * z[c];
        if (!taken[c] && (best == M || f < best_value)) {
          best = c;
          best_value = f;
        }
      }
      taken[best] = true;
      chosen.push_back(best);
    }
  }

  Func const &cost_function_;
  turbo_parameters<ArgumentContainer> const &params_;
  URBG &gen_;
  ResultType target_value_;
  std::atomic<bool> *cancellation_;
  std::atomic<ResultType> *current_minimum_cost_;
  std::vector<std::pair<ArgumentContainer, ResultType>> *queries_;
  size_t n_;
  size_t q_;
  size_t max_function_calls_;
  size_t initial_points_;
  size_t candidates_;
  size_t failure_tolerance_;
  size_t calls_ = 0;
  ArgumentContainer best_argument_;
  ResultType best_value_ = std::numeric_limits<ResultType>::infinity();
  // The trust region: its base side length, counters, and the points evaluated since it started (unit cube):
  double length_ = length_init_;
  size_t successes_ = 0;
  size_t failures_ = 0;
  std::vector<double> X_;
  std::vector<double> y_;
  // The Gaussian process: training set, hyperparameters, and scratch.
  size_t N_ = 0;
  bool finite_ = false;
  std::vector<double> Xt_, yt_;
  std::vector<double> lengthscale_, inv_l2_;
  double outputscale_ = 1, noise_ = 0.005, mean_ = 0;
  std::vector<double> sqdiff_, pair_m_, pair_dm_;
  std::vector<double> L_, Z_, Kinv_, alpha_;
  std::vector<double> V_, Vt_, S_, covariance_;
};

} // namespace detail

BOOST_MATH_EXPORT template <typename ArgumentContainer, class Func, class URBG>
ArgumentContainer turbo(
    const Func cost_function, turbo_parameters<ArgumentContainer> const &params, URBG &gen,
    std::invoke_result_t<Func, ArgumentContainer> target_value =
        std::numeric_limits<std::invoke_result_t<Func, ArgumentContainer>>::quiet_NaN(),
    std::atomic<bool> *cancellation = nullptr,
    std::atomic<std::invoke_result_t<Func, ArgumentContainer>> *current_minimum_cost = nullptr,
    std::vector<std::pair<ArgumentContainer, std::invoke_result_t<Func, ArgumentContainer>>> *queries = nullptr) {
  detail::turbo_engine<ArgumentContainer, Func, URBG> engine(cost_function, params, gen, target_value, cancellation,
                                                             current_minimum_cost, queries);
  return engine.optimize();
}

} BOOST_MATH_NAMESPACE_END // namespace boost::math::optimization
#endif
