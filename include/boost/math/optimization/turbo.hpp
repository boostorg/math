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
#include <boost/math/optimization/detail/gaussian_process.hpp>
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
    gaussian_process_parameters gp_params;
    gp_params.hyperparameter_steps = params.hyperparameter_steps;
    gp_ = gaussian_process(gp_params);
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
    bool const model = finite_ && gp_.fit(Xt_.data(), yt_.data(), N_, n_);

    // The trust region (Section 2): side lengths L l_i / (prod_j l_j)^(1/n), centred at the best point, clipped to the
    // unit cube. The reference first divides by the mean lengthscale, for stability.
    std::vector<double> weights(n_, 1.0);
    if (model) {
      double mean = 0;
      for (size_t i = 0; i < n_; ++i) {
        mean += gp_.lengthscales()[i];
      }
      mean /= static_cast<double>(n_);
      double log_prod = 0;
      for (size_t i = 0; i < n_; ++i) {
        weights[i] = gp_.lengthscales()[i] / mean;
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

  // q joint posterior samples over the M candidates C (row-major M x n), and the argmin of each.
  void thompson_sample(std::vector<double> const &C, size_t M, std::vector<size_t> &chosen) {
    using std::sqrt;
    std::vector<double> mu;
    bool const factored = gp_.posterior(C.data(), M, mu, S_, true);
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
        double const f = factored ? mu[c] + dot(S_.data() + c * M, z.data(), c + 1)
                                  : mu[c] + sqrt(gp_.outputscale() + gp_.noise()) * z[c];
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
  // The Gaussian process, its training set, and the factor of the candidates' posterior covariance.
  size_t N_ = 0;
  bool finite_ = false;
  std::vector<double> Xt_, yt_;
  gaussian_process gp_;
  std::vector<double> S_;
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
