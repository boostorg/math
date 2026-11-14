/*
 * Copyright Nick Thompson, 2026
 * Use, modification and distribution are subject to the
 * Boost Software License, Version 1.0. (See accompanying file
 * LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef BOOST_MATH_OPTIMIZATION_QLOGNEHVI_HPP
#define BOOST_MATH_OPTIMIZATION_QLOGNEHVI_HPP

#include <boost/math/tools/config.hpp>
#ifndef BOOST_MATH_BUILD_MODULE
#include <atomic>
#endif
#include <boost/math/optimization/detail/bounded_lbfgs.hpp>
#include <boost/math/optimization/detail/common.hpp>
#include <boost/math/optimization/detail/gaussian_process.hpp>
#include <boost/math/optimization/detail/hypervolume.hpp>
#ifndef BOOST_MATH_BUILD_MODULE
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <exception>
#include <limits>
#include <mutex>
#include <numeric>
#include <random>
#include <sstream>
#include <stdexcept>
#include <thread>
#include <type_traits>
#include <utility>
#include <vector>
#endif

BOOST_MATH_NAMESPACE_BEGIN namespace optimization {

// Multi-objective Bayesian optimization by the parallel noisy expected hypervolume improvement, in its logarithmic
// form:
// Daulton, S., Balandat, M., Bakshy, E. (2021). Parallel Bayesian optimization of multiple noisy objectives with
// expected hypervolume improvement. NeurIPS 34 ("the paper": the acquisition function, its cached box
// decompositions, and the sequential greedy batches);
// Daulton, S., Balandat, M., Bakshy, E. (2020). Differentiable expected hypervolume improvement for parallel
// multi-objective Bayesian optimization. NeurIPS 33;
// Ament, S., Daulton, S., Eriksson, D., Balandat, M., Bakshy, E. (2023). Unexpected improvements to expected
// improvement for Bayesian optimization. NeurIPS 36 (the logarithm and the smoothing which make it easy to maximize).
// The defaults, and the choices the papers leave open, are those of BoTorch 0.18 (qLogNoisyExpectedHypervolumeImprovement
// with prune_baseline, SingleTaskGP, optimize_acqf with sequential = True), "the reference" below. Every objective is
// minimized, and the reference point bounds the objectives of interest from above.
//
// Where this differs from the reference:
// - the base samples are pseudo-random normal numbers in antithetic pairs, and the initial design and the starting
//   points of the maximization are a Latin hypercube and uniform random points, where it has scrambled Sobol points;
// - each restart of the maximization is its own L-BFGS run (projected, with backtracking), where it runs L-BFGS-B on
//   the sum of the acquisition function over the restarts;
// - the hypervolumes are those of the standardized objectives, so the acquisition function differs by a constant,
//   the logarithm of the product of the standard deviations, and its smoothing parameters are in those units;
// - the fit of the models follows the path of L-BFGS-B only until a hyperparameter meets a bound, and then goes on
//   over the logarithms of the hyperparameters (see gaussian_process_parameters::lbfgsb_path);
// - four fits of the models in five start from the hyperparameters of the fit before, where it always starts from
//   the modes of the priors; a fit that ends with a greater loss than the one before it is done again from there;
// - an objective that is not finite is replaced, for the models, by the worst finite value seen.

namespace detail {

// Calls f(i) for i < count on up to `threads` threads, of which the calling thread is one. The first exception thrown,
// by f or by the start of a thread, is rethrown once every thread has been joined.
template <class F> void qlognehvi_parallel(size_t count, size_t threads, F &&f) {
  threads = (std::max)(size_t(1), (std::min)(threads, count));
  std::mutex mutex;
  std::exception_ptr error;
  auto guarded = [&](auto &&g) {
    try {
      g();
    } catch (...) {
      std::lock_guard<std::mutex> lock(mutex);
      error = error ? error : std::current_exception();
    }
  };
  auto work = [&](size_t j) {
    guarded([&] {
      for (size_t i = j; i < count; i += threads) {
        f(i);
      }
    });
  };
  std::vector<std::thread> pool;
  guarded([&] {
    for (size_t j = 1; j < threads; ++j) {
      pool.emplace_back(work, j);
    }
  });
  work(0);
  for (auto &thread : pool) {
    thread.join();
  }
  if (error) {
    std::rethrow_exception(error);
  }
}

// rows x columns standard normal numbers in antithetic pairs: row 2i + 1 is minus row 2i. The reference has
// scrambled Sobol points through the inverse normal distribution function instead.
template <class URBG> void antithetic_normals(std::vector<double> &z, size_t rows, size_t columns, URBG &gen) {
  std::normal_distribution<double> normal(0, 1);
  z.resize(rows * columns);
  for (size_t r = 0; r < rows; ++r) {
    for (size_t c = 0; c < columns; ++c) {
      z[r * columns + c] = r % 2 ? -z[(r - 1) * columns + c] : normal(gen);
    }
  }
}

// The acquisition function of the paper for one candidate, in the form of Ament et al.:
//   A(x) = log (1/N) sum_s HVI(f_s(x); P_s),
// where f_s, s < N, are joint samples of the Gaussian processes of the objectives (independent of one another) at the
// conditioning set B and at x; P_s is the Pareto front of f_s(B); and HVI(y; P) is the smoothed hypervolume that y
// adds to that of P, log_hypervolume_improvement in hypervolume.hpp. B is the baseline, the points evaluated so far
// that may be optimal, and then the points of the batch that are already chosen. The region that P_s does not
// dominate is cut into boxes once per sample (Section 5 of the paper: the "cached box decompositions"), so that an
// evaluation costs a sum over boxes and no decomposition. The base samples are fixed (the sample average
// approximation), and with the sample paths of gaussian_process, for objective k
//   f_{s,k}(x) = mean_k + sum_a kernel_k(X_a, x) a_{s,k,a} + sum_j kernel_k(B_j, x) w_{s,k,j} + d_k(x) z_{s,k},
// a smooth function of x, whose gradient is that of the kernels and of the conditional standard deviation d_k.
// Everything is in the units of the models: the unit cube, and whatever the values they were fitted to were in.
BOOST_MATH_TEST_EXPORT class qlognehvi_acquisition {
public:
  // Scratch for an evaluation; concurrent evaluations need one each.
  struct workspace {
    std::vector<double> y, dy, c, sd, dvar, weight, sum, work;
  };

  // The models must be fitted, to the same points, and outlive this. tau_relu and tau_max are the smoothing
  // parameters of log_hypervolume_improvement, with the defaults of the reference.
  qlognehvi_acquisition(std::vector<gaussian_process> &models, std::vector<double> reference, size_t samples,
                        double tau_relu = 1e-6, double tau_max = 1e-3)
      : models_(models), reference_(std::move(reference)), N_(samples), tau_relu_(tau_relu), tau_max_(tau_max),
        a_(models.size()), w_(models.size()), boxes_(samples) {}

  // Conditions the samples on the `count` points B (count x n): sample s of objective k has the values
  // conditioning_mean() + L_B zeta_{s,k} there, where zeta_{s,k} is the row (k N + s) of zeta (N objectives x
  // count), and z[k N + s] for the candidate. Cuts up the region that each sample's front does not dominate, on up to
  // `threads` threads. Returns false if a covariance cannot be factored; there is then no acquisition function.
  bool condition(std::vector<double> const &B, size_t count, std::vector<double> const &zeta,
                 std::vector<double> const &z, size_t threads = 1) {
    size_t const m = models_.size(), N = N_, n = models_[0].size();
    count_ = count;
    z_ = z;
    values_.assign(N, std::vector<std::vector<double>>(count, std::vector<double>(m)));
    for (size_t k = 0; k < m; ++k) {
      auto &gp = models_[k];
      if (!gp.condition(B.data(), count)) {
        return false;
      }
      a_[k].resize(N * n);
      w_[k].resize(N * count);
      double const *L = gp.conditioning_factor().data();
      for (size_t s = 0; s < N; ++s) {
        double const *row = zeta.data() + (k * N + s) * count;
        for (size_t j = 0; j < count; ++j) {
          values_[s][j][k] = gp.conditioning_mean()[j] + dot(L + j * count, row, j + 1);
        }
        gp.sample_path(row, a_[k].data() + s * n, w_[k].data() + s * count);
      }
    }
    qlognehvi_parallel(N, threads, [&](size_t s) { boxes_[s] = nondominated_boxes(values_[s], reference_); });
    return true;
  }

  // A(x), and unless gradient is null its gradient (n numbers). The cost is that of the conditional variance of each
  // objective, O(N^2 + |B|^2) for N training points, then O(samples (N + |B|)) for each, and O(objectives) for each
  // box of each sample. Allocates nothing once the workspace has been used.
  double operator()(double const *x, double *gradient, workspace &ws) const {
    using std::exp;
    using std::log;
    using std::sqrt;
    size_t const m = models_.size(), N = N_, nb = count_, n = models_[0].size(), d = models_[0].dimension();
    size_t const stride = n + nb;
    ws.y.resize(N * m);
    ws.dy.resize(N * m);
    ws.c.resize(2 * m * stride);
    ws.sd.resize(m);
    ws.dvar.resize(m * d);
    ws.weight.resize(N);
    ws.sum.resize(stride);
    for (size_t k = 0; k < m; ++k) {
      auto const &gp = models_[k];
      // The kernel between x and the training points, then B; and the factors of its gradient.
      double *kx = ws.c.data() + 2 * k * stride, *cx = kx + stride;
      double const *X = gp.training_points().data(), *B = gp.conditioning_points().data();
      for (size_t a = 0; a < stride; ++a) {
        kx[a] = gp.kernel(a < n ? X + a * d : B + (a - n) * d, x, cx + a);
      }
      ws.sd[k] = sqrt(gp.conditional_variance(x, gradient ? ws.dvar.data() + k * d : nullptr, ws.work));
      for (size_t s = 0; s < N; ++s) {
        ws.y[s * m + k] = gp.mean() + dot(a_[k].data() + s * n, kx, n) + dot(w_[k].data() + s * nb, kx + n, nb) +
                          ws.sd[k] * z_[k * N + s];
      }
    }
    double largest = -std::numeric_limits<double>::infinity();
    for (size_t s = 0; s < N; ++s) {
      ws.weight[s] = log_hypervolume_improvement(ws.y.data() + s * m, boxes_[s],
                                                 gradient ? ws.dy.data() + s * m : nullptr, tau_relu_, tau_max_);
      largest = (std::max)(largest, ws.weight[s]);
    }
    double total = 0;
    for (size_t s = 0; s < N; ++s) {
      ws.weight[s] = exp(ws.weight[s] - largest);
      total += ws.weight[s];
    }
    if (gradient) {
      std::fill(gradient, gradient + d, 0.0);
    }
    for (size_t k = 0; gradient && k < m; ++k) {
      auto const &gp = models_[k];
      // dA/df_{s,k} = g_s, so by the chain rule the kernel at point a has the weight sum_s g_s a_{s,k,a}:
      double *sum = ws.sum.data();
      std::fill(sum, sum + stride, 0.0);
      double sum_z = 0;
      for (size_t s = 0; s < N; ++s) {
        double const g = ws.weight[s] / total * ws.dy[s * m + k];
        double const *as = a_[k].data() + s * n, *bs = w_[k].data() + s * nb;
        for (size_t a = 0; a < n; ++a) {
          sum[a] += g * as[a];
        }
        for (size_t j = 0; j < nb; ++j) {
          sum[n + j] += g * bs[j];
        }
        sum_z += g * z_[k * N + s];
      }
      // d kernel(y, x)/dx_i = -c (x_i - y_i)/l_i^2, and d sqrt(v) = dv/(2 sqrt(v)):
      double const *cx = ws.c.data() + (2 * k + 1) * stride, *inv_l2 = gp.inverse_squared_lengthscales().data();
      double const *X = gp.training_points().data(), *B = gp.conditioning_points().data();
      double const *dvar = ws.dvar.data() + k * d;
      for (size_t a = 0; a < stride; ++a) {
        double const *ya = a < n ? X + a * d : B + (a - n) * d;
        double const t = cx[a] * sum[a];
        for (size_t i = 0; i < d; ++i) {
          gradient[i] -= t * (x[i] - ya[i]) * inv_l2[i];
        }
      }
      for (size_t i = 0; i < d; ++i) {
        gradient[i] += sum_z * dvar[i] / (2 * ws.sd[k]);
      }
    }
    return largest + log(total / static_cast<double>(N));
  }

private:
  std::vector<gaussian_process> &models_;
  std::vector<double> reference_;
  size_t N_, count_ = 0;
  double tau_relu_, tau_max_;
  // For each objective, the coefficients of the sample paths (N x training points, N x count), and the normal numbers
  // for the candidate; for each sample, its values at B and its boxes.
  std::vector<std::vector<double>> a_, w_;
  std::vector<double> z_;
  std::vector<std::vector<std::vector<double>>> values_;
  std::vector<box_decomposition<double>> boxes_;
};

// The points (rows of f, count x m, all objectives minimized) which are less than the reference in every objective
// and which no other point dominates: kept[j] is set for those, and left alone for the others. order is scratch.
BOOST_MATH_TEST_EXPORT inline void mark_nondominated(double const *f, size_t count, size_t m, double const *reference,
                              std::vector<size_t> &order, std::vector<char> &kept) {
  order.clear();
  for (size_t j = 0; j < count; ++j) {
    size_t k = 0;
    while (k < m && f[j * m + k] < reference[k]) {
      ++k;
    }
    if (k == m) {
      order.push_back(j);
    }
  }
  // A point comes after those which dominate it, so it need only be compared with the front so far, which is moved
  // to the beginning of order.
  std::sort(order.begin(), order.end(), [&](size_t a, size_t b) {
    return std::lexicographical_compare(f + a * m, f + (a + 1) * m, f + b * m, f + (b + 1) * m);
  });
  size_t front = 0;
  for (size_t j : order) {
    bool dominated = false;
    for (size_t i = 0; i < front && !dominated; ++i) {
      double const *p = f + order[i] * m, *q = f + j * m;
      bool worse = false, better = false;
      for (size_t k = 0; k < m; ++k) {
        worse = worse || p[k] > q[k];
        better = better || p[k] < q[k];
      }
      dominated = better && !worse;
    }
    if (!dominated) {
      order[front++] = j;
      kept[j] = 1;
    }
  }
}

} // namespace detail

BOOST_MATH_EXPORT template <typename ArgumentContainer> struct qlognehvi_parameters {
  ArgumentContainer lower_bounds;
  ArgumentContainer upper_bounds;
  // The points of the initial design, a Latin hypercube; 0 means 2 (dimension + 1).
  size_t initial_samples = 0;
  // q, the number of points proposed per fit of the models, one after the other, each given those before it.
  size_t batch_size = 1;
  // The evaluations of qlognehvi(), the initial design among them. Not used by qlognehvi_optimizer.
  size_t max_evaluations = 100;
  // The joint samples of the models over which the hypervolume improvement is averaged.
  size_t mc_samples = 64;
  // The acquisition function is maximized by L-BFGS from `restarts` of `raw_samples` random points.
  size_t restarts = 4;
  size_t raw_samples = 128;
  // An upper bound for each objective: only points below it in every objective count. If empty, it is the greatest
  // value of each objective over the initial design, plus a tenth of the objective's range there.
  std::vector<double> reference_point;
  // Members of the initial design, in place of points of the Latin hypercube, or in addition to initial_samples
  // if there are more:
  std::vector<ArgumentContainer> const *initial_population = nullptr;
  unsigned threads = std::thread::hardware_concurrency();
};

BOOST_MATH_EXPORT template <typename ArgumentContainer>
void validate_qlognehvi_parameters(qlognehvi_parameters<ArgumentContainer> const &params) {
  std::ostringstream oss;
  detail::validate_bounds(params.lower_bounds, params.upper_bounds);
  if (params.batch_size == 0 || params.max_evaluations == 0 || params.mc_samples == 0) {
    oss << __FILE__ << ":" << __LINE__ << ":" << __func__;
    oss << ": The batch size, the number of evaluations and the number of samples must be at least 1, but got "
        << params.batch_size << ", " << params.max_evaluations << " and " << params.mc_samples << ".";
    throw std::invalid_argument(oss.str());
  }
  if (params.restarts == 0 || params.raw_samples < params.restarts) {
    oss << __FILE__ << ":" << __LINE__ << ":" << __func__;
    oss << ": There must be at least one restart, and as many raw samples as restarts, but restarts = "
        << params.restarts << " and raw_samples = " << params.raw_samples << ".";
    throw std::invalid_argument(oss.str());
  }
  if (params.reference_point.size() == 1 ||
      !std::all_of(params.reference_point.begin(), params.reference_point.end(), [](double r) { return std::isfinite(r); })) {
    oss << __FILE__ << ":" << __LINE__ << ":" << __func__;
    oss << ": The reference point must have one finite number for each of at least two objectives, or be empty.";
    throw std::domain_error(oss.str());
  }
  if (params.initial_population) {
    if (params.initial_population->size() > params.max_evaluations) {
      oss << __FILE__ << ":" << __LINE__ << ":" << __func__;
      oss << ": The initial population has " << params.initial_population->size()
          << " members, but max_evaluations is " << params.max_evaluations << ".";
      throw std::invalid_argument(oss.str());
    }
    for (auto const &x : *params.initial_population) {
      detail::validate_initial_guess(x, params.lower_bounds, params.upper_bounds);
    }
  }
  if (params.threads == 0) {
    oss << __FILE__ << ":" << __LINE__ << ":" << __func__;
    oss << ": There must be at least one thread.";
    throw std::invalid_argument(oss.str());
  }
}

// qLogNEHVI in ask and tell form, for objectives the caller evaluates. ask() returns the points to evaluate, first
// the initial design and then batches of batch_size points, and tell(values) takes their objectives, in the same
// order. front() is the nondominated points so far. ask() fits the models and maximizes the acquisition function on
// up to params.threads threads; all random numbers are drawn on the calling thread, from the generator given at
// construction, of which the optimizer keeps a copy unless URBG is a reference type, so that the points do not depend
// on the number of threads.
//
// An objective which is not finite at a point (NaN or infinite) is, for the models, the greatest finite value seen
// of that objective, so that the point is avoided but does not distort them; such a point is never in front().
// While an objective has no finite value, or if a model cannot be fitted, ask() returns uniformly distributed points.
BOOST_MATH_EXPORT template <typename ArgumentContainer, typename ObjectiveContainer, class URBG = std::mt19937_64>
class qlognehvi_optimizer {
public:
  qlognehvi_optimizer(qlognehvi_parameters<ArgumentContainer> const &params, URBG gen)
      : params_(params), gen_(std::forward<URBG>(gen)), n_(params.lower_bounds.size()) {
    validate_qlognehvi_parameters(params_);
    if (params_.initial_population) {
      pending_ = *params_.initial_population;
      params_.initial_population = nullptr;
    }
    reference_ = params_.reference_point;
  }

  // The same points are returned until their objectives are told.
  std::vector<ArgumentContainer> const &ask() { return ask(params_.batch_size); }

  // A batch of `count` points rather than batch_size; the initial design is what it is.
  std::vector<ArgumentContainer> const &ask(size_t count) {
    if (asked_) {
      return pending_;
    }
    asked_ = true;
    std::uniform_real_distribution<double> unif01(0, 1);
    if (points_.empty()) {
      // The points given, then a Latin hypercube: one point in each of `rest` equal slices of every coordinate.
      size_t const given = pending_.size();
      size_t const initial = params_.initial_samples ? params_.initial_samples : 2 * (n_ + 1);
      size_t const rest = initial > given ? initial - given : 0;
      std::vector<double> design(rest * n_);
      std::vector<size_t> slices(rest);
      for (size_t i = 0; i < n_; ++i) {
        std::iota(slices.begin(), slices.end(), size_t(0));
        std::shuffle(slices.begin(), slices.end(), gen_);
        for (size_t k = 0; k < rest; ++k) {
          design[k * n_ + i] = (static_cast<double>(slices[k]) + unif01(gen_)) / static_cast<double>(rest);
        }
      }
      for (size_t k = 0; k < rest; ++k) {
        pending_.push_back(to_argument(design.data() + k * n_));
      }
      return pending_;
    }
    std::vector<double> u(count * n_);
    for (size_t i = reference_.empty() ? 0 : propose(count, u); i < count; ++i) {
      for (size_t j = 0; j < n_; ++j) {
        u[i * n_ + j] = unif01(gen_);
      }
    }
    for (size_t i = 0; i < count; ++i) {
      pending_.push_back(to_argument(u.data() + i * n_));
    }
    return pending_;
  }

  // The objectives of the points last asked for, in their order.
  void tell(std::vector<ObjectiveContainer> values) {
    using std::isfinite;
    const size_t objectives = values_.empty() ? (values.empty() ? 0 : values[0].size()) : values_[0].size();
    auto const bad = std::find_if(values.begin(), values.end(), [&](auto const &v) { return v.size() != objectives; });
    if (!asked_ || values.size() != pending_.size() || objectives < 2 || bad != values.end() ||
        (!params_.reference_point.empty() && params_.reference_point.size() != objectives)) {
      std::ostringstream oss;
      oss << __FILE__ << ":" << __LINE__ << ":" << __func__;
      oss << ": tell takes the objectives of the " << (asked_ ? pending_.size() : 0) << " points last asked for, "
          << objectives << " (at least two) for each, but got those of " << values.size() << " points";
      if (bad != values.end()) {
        oss << ", " << bad->size() << " for one of them";
      }
      if (!params_.reference_point.empty()) {
        oss << "; the reference point has " << params_.reference_point.size();
      }
      oss << ".";
      throw std::domain_error(oss.str());
    }
    for (size_t i = 0; i < values.size(); ++i) {
      for (size_t j = 0; j < n_; ++j) {
        auto const width = params_.upper_bounds[j] - params_.lower_bounds[j];
        unit_.push_back(width > 0 ? static_cast<double>((pending_[i][j] - params_.lower_bounds[j]) / width) : 0.0);
      }
      points_.push_back(std::move(pending_[i]));
      values_.push_back(std::move(values[i]));
    }
    pending_.clear();
    asked_ = false;
    if (reference_.empty()) {
      // The nadir of what has been seen, and a tenth of the range (of the magnitude, if there is no range) beyond.
      std::vector<double> reference(objectives);
      for (size_t k = 0; k < objectives; ++k) {
        double least = std::numeric_limits<double>::infinity(), greatest = -least;
        for (auto const &v : values_) {
          double const y = static_cast<double>(v[k]);
          least = isfinite(y) ? (std::min)(least, y) : least;
          greatest = isfinite(y) ? (std::max)(greatest, y) : greatest;
        }
        if (!(least <= greatest)) {
          return;
        }
        double const range = greatest - least;
        reference[k] = greatest + (range > 0 ? range : (std::max)(std::abs(greatest), 1.0)) / 10;
      }
      reference_ = std::move(reference);
    }
  }

  // The points whose objectives are all finite and which none of the others dominates, as (argument, objectives)
  // pairs in lexicographic order of the objectives, those with equal objectives once. Empty until the first tell.
  std::vector<std::pair<ArgumentContainer, ObjectiveContainer>> front() const {
    using std::isfinite;
    std::vector<std::pair<ArgumentContainer, ObjectiveContainer>> front;
    auto finite = [](ObjectiveContainer const &v) {
      return std::all_of(v.begin(), v.end(), [](auto const &y) { return isfinite(y); });
    };
    for (size_t i = 0; i < points_.size(); ++i) {
      bool dominated = !finite(values_[i]);
      for (size_t j = 0; j < points_.size() && !dominated; ++j) {
        bool better = false, worse = false;
        for (size_t k = 0; k < values_[i].size(); ++k) {
          better = better || values_[j][k] < values_[i][k];
          worse = worse || values_[j][k] > values_[i][k];
        }
        dominated = better && !worse && finite(values_[j]);
      }
      if (!dominated) {
        front.emplace_back(points_[i], values_[i]);
      }
    }
    std::stable_sort(front.begin(), front.end(), [](auto const &a, auto const &b) {
      return std::lexicographical_compare(a.second.begin(), a.second.end(), b.second.begin(), b.second.end());
    });
    auto same = [](auto const &a, auto const &b) { return a.second == b.second; };
    front.erase(std::unique(front.begin(), front.end(), same), front.end());
    return front;
  }

  // Every point told so far, in order, and their objectives.
  std::vector<ArgumentContainer> const &points() const { return points_; }
  std::vector<ObjectiveContainer> const &objectives() const { return values_; }
  // Empty until it is known.
  std::vector<double> const &reference_point() const { return reference_; }

private:
  using Real = typename ArgumentContainer::value_type;

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

  // Writes up to `count` points of the unit cube to u, each the maximizer of the acquisition function given those
  // before it, and returns how many: fewer only if there is no model.
  size_t propose(size_t count, std::vector<double> &u) {
    using std::exp;
    using std::isfinite;
    using std::sqrt;
    size_t const N = points_.size(), m = values_[0].size(), threads = params_.threads;
    // The models: one Gaussian process for each objective, on the unit cube, fitted to the values standardized by
    // their mean and standard deviation as Standardize does in the reference.
    std::vector<double> y(m * N), reference(m);
    for (size_t k = 0; k < m; ++k) {
      double worst = -std::numeric_limits<double>::infinity();
      for (size_t i = 0; i < N; ++i) {
        y[k * N + i] = static_cast<double>(values_[i][k]);
        worst = isfinite(y[k * N + i]) ? (std::max)(worst, y[k * N + i]) : worst;
      }
      if (!isfinite(worst)) {
        return 0;
      }
      double mean = 0, variance = 0;
      for (size_t i = 0; i < N; ++i) {
        y[k * N + i] = isfinite(y[k * N + i]) ? y[k * N + i] : worst;
        mean += y[k * N + i] / static_cast<double>(N);
      }
      for (size_t i = 0; i < N; ++i) {
        variance += (y[k * N + i] - mean) * (y[k * N + i] - mean);
      }
      double const sd = N > 1 && sqrt(variance / static_cast<double>(N - 1)) >= 1e-8
                            ? sqrt(variance / static_cast<double>(N - 1))
                            : 1.0;
      for (size_t i = 0; i < N; ++i) {
        y[k * N + i] = (y[k * N + i] - mean) / sd;
      }
      reference[k] = (reference_[k] - mean) / sd;
    }
    // A fit from the hyperparameters of the last one costs a fraction of a fit from the modes of the priors, but
    // which optimum a fit ends in depends on where it starts, and with warm fits alone the points proposed are worse.
    // So every fifth fit is from the modes, and so is any after which the loss is greater than after the last one.
    bool const warm = models_.size() == m && fits_ % 5 != 0;
    ++fits_;
    if (!warm) {
      models_.assign(m, detail::gaussian_process(detail::gaussian_process_parameters::botorch(n_)));
    }
    losses_.resize(m);
    std::vector<char> fitted(m), factored(m);
    detail::qlognehvi_parallel(m, threads, [&](size_t k) {
      std::vector<double> gradient(n_ + 3);
      double loss = 0;
      for (int cold = warm ? 0 : 1; cold < 2; ++cold) {
        if (cold && warm) {
          models_[k] = detail::gaussian_process(detail::gaussian_process_parameters::botorch(n_));
        }
        fitted[k] = models_[k].fit(unit_.data(), y.data() + k * N, N, n_, !cold) &&
                    models_[k].loss_and_gradient(models_[k].hyperparameters(), loss, gradient);
        if (fitted[k] && !(loss > losses_[k])) {
          break;
        }
      }
      losses_[k] = loss;
      factored[k] = fitted[k] && models_[k].condition(unit_.data(), N);
    });
    if (std::find(fitted.begin(), fitted.end(), 0) != fitted.end()) {
      return 0;
    }
    // The baseline (prune_inferior_points_multi_objective in the reference): the points which are on the front, and
    // better than the reference point, in at least one of 2048 joint samples of the models at the points evaluated.
    // If the models' covariance there cannot be factored, they are all kept.
    std::vector<char> kept(N, 1);
    if (std::find(factored.begin(), factored.end(), 0) == factored.end()) {
      size_t const pairs = 1024, workers = (std::min)(threads, pairs);
      detail::antithetic_normals(normals_, 2 * pairs, m * N, gen_);
      std::vector<std::vector<char>> marks(workers, std::vector<char>(N, 0));
      detail::qlognehvi_parallel(workers, workers, [&](size_t t) {
        // A sample, then the one with the opposite normal numbers, which is its reflection in the mean.
        std::vector<double> f(2 * N * m);
        std::vector<size_t> order;
        for (size_t s = 2 * t; s < 2 * pairs; s += 2 * workers) {
          for (size_t k = 0; k < m; ++k) {
            double const *zeta = normals_.data() + (s * m + k) * N;
            double const *L = models_[k].conditioning_factor().data();
            for (size_t j = 0; j < N; ++j) {
              double const step = detail::dot(L + j * N, zeta, j + 1);
              f[j * m + k] = models_[k].conditioning_mean()[j] + step;
              f[(N + j) * m + k] = models_[k].conditioning_mean()[j] - step;
            }
          }
          detail::mark_nondominated(f.data(), N, m, reference.data(), order, marks[t]);
          detail::mark_nondominated(f.data() + N * m, N, m, reference.data(), order, marks[t]);
        }
      });
      for (size_t j = 0; j < N; ++j) {
        kept[j] = 0;
        for (auto const &mark : marks) {
          kept[j] = kept[j] || mark[j];
        }
      }
    }
    std::vector<double> B;
    for (size_t j = 0; j < N; ++j) {
      if (kept[j]) {
        B.insert(B.end(), unit_.begin() + static_cast<std::ptrdiff_t>(j * n_),
                 unit_.begin() + static_cast<std::ptrdiff_t>((j + 1) * n_));
      }
    }
    size_t nb = B.size() / n_;
    size_t const S = params_.mc_samples, raw = params_.raw_samples, restarts = params_.restarts;
    // zeta[(k S + s) nb + j] is the base sample of objective k at point j of B in sample s.
    std::vector<double> zeta, z, column, points(raw * n_), values(raw), lower(n_, 0.0), upper(n_);
    detail::antithetic_normals(normals_, S, m * nb, gen_);
    zeta.resize(m * S * nb);
    for (size_t s = 0; s < S; ++s) {
      for (size_t k = 0; k < m; ++k) {
        std::copy_n(normals_.data() + (s * m + k) * nb, nb, zeta.data() + (k * S + s) * nb);
      }
    }
    for (size_t i = 0; i < n_; ++i) {
      upper[i] = params_.upper_bounds[i] > params_.lower_bounds[i] ? 1.0 : 0.0;
    }
    detail::qlognehvi_acquisition acquisition(models_, reference, S);
    std::vector<detail::qlognehvi_acquisition::workspace> scratch((std::min)(threads, raw));
    std::uniform_real_distribution<double> unif01(0, 1);
    for (size_t c = 0; c < count; ++c) {
      // The samples at the candidate, then (the transpose of) those at the last point of B, if it is new.
      detail::antithetic_normals(normals_, S, m, gen_);
      z.resize(m * S);
      for (size_t s = 0; s < S; ++s) {
        for (size_t k = 0; k < m; ++k) {
          z[k * S + s] = normals_[s * m + k];
        }
      }
      if (!acquisition.condition(B, nb, zeta, z, threads)) {
        return c;
      }
      // The starting points (gen_batch_initial_conditions and initialize_q_batch in the reference): of raw_samples
      // random points, `restarts` are drawn without replacement with probabilities proportional to exp(t), t the
      // acquisition function standardized over them; the best takes the last place if it was not drawn.
      for (size_t i = 0; i < raw * n_; ++i) {
        points[i] = upper[i % n_] * unif01(gen_);
      }
      detail::qlognehvi_parallel(scratch.size(), scratch.size(), [&](size_t t) {
        for (size_t i = t; i < raw; i += scratch.size()) {
          values[i] = acquisition(points.data() + i * n_, nullptr, scratch[t]);
        }
      });
      std::vector<size_t> starts(raw);
      std::iota(starts.begin(), starts.end(), size_t(0));
      double mean = 0, variance = 0;
      for (double v : values) {
        mean += v / static_cast<double>(raw);
      }
      for (double v : values) {
        variance += (v - mean) * (v - mean) / static_cast<double>(raw > 1 ? raw - 1 : 1);
      }
      if (raw > restarts && variance > 0 && isfinite(variance)) {
        std::vector<double> weights(raw);
        for (size_t i = 0; i < raw; ++i) {
          weights[i] = exp((values[i] - mean) / sqrt(variance));
        }
        size_t const best = static_cast<size_t>(std::max_element(values.begin(), values.end()) - values.begin());
        bool has_best = false;
        for (size_t r = 0; r < restarts; ++r) {
          double total = 0;
          for (size_t i = r; i < raw; ++i) {
            total += weights[starts[i]];
          }
          double target = total * unif01(gen_);
          size_t pick = r;
          while (pick + 1 < raw && target >= weights[starts[pick]]) {
            target -= weights[starts[pick++]];
          }
          std::swap(starts[r], starts[pick]);
          has_best = has_best || starts[r] == best;
        }
        if (!has_best) {
          starts[restarts - 1] = best;
        }
      }
      starts.resize(restarts);
      // L-BFGS from each, for at most 100 iterations as the reference was run; the best end point is the candidate.
      std::vector<std::vector<double>> ends(restarts);
      std::vector<double> end_values(restarts);
      size_t const workers = (std::min)(threads, restarts);
      detail::qlognehvi_parallel(workers, workers, [&](size_t t) {
        detail::bounded_lbfgs_parameters lbfgs;
        lbfgs.max_iterations = 100;
        lbfgs.maximize = true;
        for (size_t r = t; r < restarts; r += workers) {
          ends[r].assign(points.begin() + static_cast<std::ptrdiff_t>(starts[r] * n_),
                         points.begin() + static_cast<std::ptrdiff_t>((starts[r] + 1) * n_));
          auto f = [&](std::vector<double> const &x, std::vector<double> &g) {
            return acquisition(x.data(), g.data(), scratch[t]);
          };
          end_values[r] = detail::bounded_lbfgs(f, ends[r], lower, upper, lbfgs).value;
        }
      });
      size_t best = 0;
      for (size_t r = 1; r < restarts; ++r) {
        best = end_values[r] > end_values[best] ? r : best;
      }
      std::copy(ends[best].begin(), ends[best].end(), u.begin() + static_cast<std::ptrdiff_t>(c * n_));
      // The candidate joins B as a pending point (Section 5 of the paper, cache_pending in the reference): each
      // sample gets values there, so that the next candidate is credited only with what it adds to them.
      if (c + 1 < count) {
        B.insert(B.end(), ends[best].begin(), ends[best].end());
        detail::antithetic_normals(normals_, S, m, gen_);
        column.resize(m * S * (nb + 1));
        for (size_t r = 0; r < m * S; ++r) {
          std::copy_n(zeta.data() + r * nb, nb, column.data() + r * (nb + 1));
          column[r * (nb + 1) + nb] = normals_[(r % S) * m + r / S];
        }
        zeta.swap(column);
        ++nb;
      }
    }
    return count;
  }

  qlognehvi_parameters<ArgumentContainer> params_;
  URBG gen_;
  size_t n_;
  // What has been evaluated: the points, in the unit cube too (row-major), and their objectives.
  std::vector<ArgumentContainer> points_;
  std::vector<double> unit_;
  std::vector<ObjectiveContainer> values_;
  std::vector<double> reference_;
  std::vector<ArgumentContainer> pending_;
  bool asked_ = false;
  std::vector<detail::gaussian_process> models_;
  // The number of fits of the models so far, and the loss of each model after the last.
  size_t fits_ = 0;
  std::vector<double> losses_;
  std::vector<double> normals_;
};

// Minimizes every component of objectives(x) over the bounds with max_evaluations evaluations: asks a
// qlognehvi_optimizer for batches, evaluates each on params.threads threads, and returns its front.
// The result does not depend on the number of threads.
BOOST_MATH_EXPORT template <typename ArgumentContainer, class Func, class URBG>
std::vector<std::pair<ArgumentContainer, std::invoke_result_t<Func, ArgumentContainer>>>
qlognehvi(const Func objectives, qlognehvi_parameters<ArgumentContainer> const &params, URBG &gen,
          std::atomic<bool> *cancellation = nullptr,
          std::vector<std::pair<ArgumentContainer, std::invoke_result_t<Func, ArgumentContainer>>> *queries = nullptr) {
  using ObjectiveContainer = std::invoke_result_t<Func, ArgumentContainer>;
  // The initial design may not be more than the evaluations allow:
  auto capped = params;
  size_t const initial = params.initial_samples ? params.initial_samples : 2 * (params.lower_bounds.size() + 1);
  capped.initial_samples = (std::min)(initial, params.max_evaluations);
  qlognehvi_optimizer<ArgumentContainer, ObjectiveContainer, URBG &> optimizer(capped, gen);
  for (size_t evaluated = 0;;) {
    auto const &points = optimizer.ask((std::min)(params.batch_size, params.max_evaluations - evaluated));
    std::vector<ObjectiveContainer> values(points.size());
    detail::qlognehvi_parallel(points.size(), params.threads, [&](size_t i) { values[i] = objectives(points[i]); });
    for (size_t i = 0; queries && i < points.size(); ++i) {
      queries->emplace_back(points[i], values[i]);
    }
    evaluated += points.size();
    optimizer.tell(std::move(values));
    if (evaluated >= params.max_evaluations || (cancellation && *cancellation)) {
      return optimizer.front();
    }
  }
}

} BOOST_MATH_NAMESPACE_END // namespace boost::math::optimization
#endif
