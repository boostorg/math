/*
 * Copyright Nick Thompson, 2026
 * Use, modification and distribution are subject to the
 * Boost Software License, Version 1.0. (See accompanying file
 * LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef BOOST_MATH_OPTIMIZATION_DETAIL_GAUSSIAN_PROCESS_HPP
#define BOOST_MATH_OPTIMIZATION_DETAIL_GAUSSIAN_PROCESS_HPP

#include <boost/math/optimization/detail/bounded_lbfgs.hpp>
#include <boost/math/tools/config.hpp>
#include <boost/math/tools/detail/linear_algebra/cholesky.hpp>
#ifndef BOOST_MATH_BUILD_MODULE
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <vector>
#endif

BOOST_MATH_NAMESPACE_BEGIN namespace optimization::detail {

// A Gaussian process with a constant mean and an ARD Matern-5/2 or Gaussian kernel, for the Bayesian optimizers. It
// began as the surrogate of TuRBO (turbo.hpp), whose reference code
// (github.com/uber-research/TuRBO, gp.py) the fit and its defaults follow; the other choice is the model of BoTorch.

BOOST_MATH_TEST_EXPORT enum class gaussian_process_kernel { matern52, rbf };

// A kernel k(x, y) of unit scale on the unit cube, where r^2 = sum_i (x_i - y_i)^2 / l_i^2 and inv_l2[i] = 1/l_i^2:
// the ARD Matern-5/2 kernel (Section C of the TuRBO paper), (1 + sqrt(5) r + 5 r^2/3) exp(-sqrt(5) r), or the ARD
// Gaussian kernel exp(-r^2/2). If c is given, *c receives the factor of the gradient with respect to x,
// d k / d x_i = -(*c) (x_i - y_i) / l_i^2, which is (5/3) (1 + sqrt(5) r) exp(-sqrt(5) r) or k: finite at x = y.
BOOST_MATH_TEST_EXPORT inline double kernel(gaussian_process_kernel kind, double const *x, double const *y,
                                            double const *inv_l2, size_t d, double *c = nullptr) {
  using std::exp;
  using std::sqrt;
  double r2 = 0;
  for (size_t i = 0; i < d; ++i) {
    double const t = x[i] - y[i];
    r2 += t * t * inv_l2[i];
  }
  if (kind == gaussian_process_kernel::rbf) {
    double const k = exp(-r2 / 2);
    if (c) {
      *c = k;
    }
    return k;
  }
  double const r = sqrt(5 * r2);
  double const e = exp(-r);
  if (c) {
    *c = (5.0 / 3) * (1 + r) * e;
  }
  return (1 + r + r * r / 3) * e;
}

// The defaults are those of TuRBO's reference.
BOOST_MATH_TEST_EXPORT struct gaussian_process_parameters {
  gaussian_process_kernel kernel = gaussian_process_kernel::matern52;
  // The hyperparameter bounds of the reference (gp.py); the paper (Section C) has 0.1 for the upper noise bound.
  double lengthscale_bounds[2] = {0.005, 2.0};
  double outputscale_bounds[2] = {0.05, 20.0};
  double noise_bounds[2] = {0.0005, 0.2};
  // Where every fit starts; the mean starts at 0.
  double initial_lengthscale = 0.5;
  double initial_outputscale = 1.0;
  double initial_noise = 0.005;
  // Adam steps on the log marginal likelihood per fit, and their learning rate.
  size_t hyperparameter_steps = 50;
  double learning_rate = 0.1;
  // Instead of Adam on sigmoids, fit by the bounded L-BFGS on the logarithms of the hyperparameters, for at most
  // hyperparameter_steps iterations: for bounds that span orders of magnitude or are infinite.
  bool lbfgs = false;
  // With L-BFGS: fit first over the hyperparameters themselves, with the line search of L-BFGS-B, and only then over
  // their logarithms. The posterior density often has several local maxima, and which of them a descent method ends
  // in depends on its path. This path is that of SciPy's L-BFGS-B, with which BoTorch fits its models, as long as no
  // hyperparameter meets a bound; from the modes of the priors it tends to end at simpler models (longer
  // lengthscales, more noise) than the fit over the logarithms alone, whose maxima are often higher, and the
  // Bayesian optimizers do better with it.
  bool lbfgsb_path = false;
  // If false, the L-BFGS fit keeps the outputscale at its initial value.
  bool fit_outputscale = true;
  // The (loc, scale) of LogNormal priors on each lengthscale and on the noise, whose log densities are added to the
  // log marginal likelihood (a maximum a posteriori fit); none if the scale is 0.
  double lengthscale_prior[2] = {0.0, 0.0};
  double noise_prior[2] = {0.0, 0.0};

  // The model that BoTorch 0.18 fits by default (SingleTaskGP with fit_gpytorch_mll) in n dimensions, for values
  // standardized by the caller: a Gaussian kernel of unit scale with the dimension-scaled lengthscale prior of
  // Hvarfner, Hellsten and Nardi, "Vanilla Bayesian optimization performs great in high dimensions", ICML 2024,
  // which keeps the model simple unless the data say otherwise; the fit starts at the modes of the priors.
  static gaussian_process_parameters botorch(size_t n) {
    using std::exp;
    using std::log;
    using std::sqrt;
    double const inf = std::numeric_limits<double>::infinity();
    gaussian_process_parameters p;
    p.kernel = gaussian_process_kernel::rbf;
    p.lbfgs = true;
    p.lbfgsb_path = true;
    p.fit_outputscale = false;
    p.hyperparameter_steps = 1000;
    p.lengthscale_bounds[0] = 0.025;
    p.lengthscale_bounds[1] = inf;
    p.noise_bounds[0] = 1e-4;
    p.noise_bounds[1] = inf;
    p.lengthscale_prior[0] = sqrt(2.0) + log(static_cast<double>(n)) / 2;
    p.lengthscale_prior[1] = sqrt(3.0);
    p.noise_prior[0] = -4;
    p.noise_prior[1] = 1;
    p.initial_lengthscale = exp(p.lengthscale_prior[0] - 3);
    p.initial_noise = exp(-5.0);
    return p;
  }
};

// The Gaussian process, on points of the unit cube stored row-major. With L L^T = K + noise I the Cholesky
// factorization of the kernel matrix of the N training points X plus the noise, and alpha = (K + noise I)^{-1}
// (y - mean), the posterior of the function has mean mean + k(x, X) alpha and covariance
//   Sigma(a, b) = k(a, b) - k(a, X) (K + noise I)^{-1} k(X, b),
// where k includes the outputscale.
BOOST_MATH_TEST_EXPORT class gaussian_process {
public:
  explicit gaussian_process(gaussian_process_parameters const &params = {}) : params_(params) {}

  // Fits the Gaussian process (Section C) to the N training points X (N x n) with finite values y by Adam on the
  // loss below, over the hyperparameters transformed to the real line as in GPyTorch: value = lo + (hi - lo)
  // sigmoid(raw); or by L-BFGS over their logarithms within the bounds, where BoTorch has L-BFGS-B over the
  // hyperparameters themselves. Starts from the initial values of the parameters, or if warm is set and there has
  // been a fit in this dimension, from its hyperparameters, and L-BFGS then goes over their logarithms only.
  // Leaves L_ and alpha_ for the posterior, and with L-BFGS the loss at the hyperparameters it ends with for loss()
  // (NaN with Adam, which does not evaluate it there, and after a fit that fails). Returns false if there is no model. Forgets the conditioning set.
  bool fit(double const *X, double const *y, size_t N, size_t n, bool warm = false) {
    using std::exp;
    using std::isfinite;
    using std::log;
    using std::sqrt;
    m_ = 0;
    n_ = n;
    N_ = N;
    Xt_.assign(X, X + N * n);
    yt_.assign(y, y + N);
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
      theta[i] = to_raw(params_.initial_lengthscale, params_.lengthscale_bounds);
    }
    theta[n_] = to_raw(params_.initial_outputscale, params_.outputscale_bounds);
    theta[n_ + 1] = to_raw(params_.initial_noise, params_.noise_bounds);
    theta[n_ + 2] = 0;
    warm = warm && theta_.size() == dim;
    theta = warm ? theta_ : theta;
    std::vector<double> m1(dim, 0.0), m2(dim, 0.0), grad(dim), last_good = theta;
    // Adam with the PyTorch defaults (the reference's learning rate is 0.1):
    double const lr = params_.learning_rate, beta1 = 0.9, beta2 = 0.999, eps = 1e-8;
    double beta1_t = 1, beta2_t = 1;
    // Set if Adam stops on hyperparameters whose kernel matrix cannot be factored, or whose loss or gradient is not
    // finite; the model then uses the last hyperparameters that were fine.
    bool bad = false;
    // With L-BFGS, first over the hyperparameters themselves if that is asked for, then over their logarithms.
    for (int phase = params_.lbfgsb_path && !warm ? 0 : 1; phase < 2 && params_.lbfgs && !bad; ++phase) {
      double const inf = std::numeric_limits<double>::infinity();
      std::vector<double> lower(dim, -inf), upper(dim, inf);
      for (size_t k = 0; k < n_ + 2; ++k) {
        double const *bounds = k < n_    ? params_.lengthscale_bounds
                               : k == n_ ? params_.outputscale_bounds
                                         : params_.noise_bounds;
        bool const fixed = k == n_ && !params_.fit_outputscale;
        auto raw = [&](double value) { return phase == 0 ? value : log(value); };
        theta[k] = linear_ == (phase == 0) ? theta[k] : raw(linear_ ? theta[k] : exp(theta[k]));
        lower[k] = raw(fixed ? params_.initial_outputscale : bounds[0]);
        upper[k] = fixed ? lower[k] : raw(bounds[1]);
      }
      linear_ = phase == 0;
      bounded_lbfgs_parameters lbfgs;
      lbfgs.max_iterations = params_.hyperparameter_steps;
      lbfgs.strong_wolfe = linear_;
      auto objective = [&](std::vector<double> const &t, std::vector<double> &g) {
        return loss_and_gradient(t, loss_, g) ? loss_ : std::numeric_limits<double>::quiet_NaN();
      };
      loss_ = bounded_lbfgs(objective, theta, lower, upper, lbfgs).value;
      bad = !isfinite(loss_);
    }
    linear_ = false;
    for (size_t step = 0; step < params_.hyperparameter_steps && !params_.lbfgs; ++step) {
      bool finite = loss_and_gradient(theta, loss_, grad) && isfinite(loss_);
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
    // The model at the final (or last good) hyperparameters. Adam does not evaluate the loss there: loss() is NaN.
    double const nan = std::numeric_limits<double>::quiet_NaN();
    loss_ = params_.lbfgs ? loss_ : nan;
    if (bad || !factor(theta)) {
      theta = last_good;
      if (!loss_and_gradient(theta, loss_, grad)) {
        loss_ = nan;
        return false;
      }
    }
    theta_ = theta;
    alpha_.resize(N);
    for (size_t a = 0; a < N; ++a) {
      alpha_[a] = yt_[a] - mean_;
    }
    tools::detail::cholesky_solve(L_.data(), N, alpha_.data());
    return true;
  }

  // The joint posterior at the M points C (M x n): its mean mu, and in the lower triangle of S (M x M) the Cholesky
  // factor of its covariance, that of the observations (with the noise) or of the function (without), with jitter
  // on its diagonal, and more if that is needed. Returns false if the covariance cannot be factored; S is then not a
  // factor.
  bool posterior(double const *C, size_t M, std::vector<double> &mu, std::vector<double> &S, bool observation_noise,
                 double jitter = 0) {
    size_t const N = N_;
    // Kxc (N x M), the posterior mean, then V = L^{-1} Kxc in place:
    V_.resize(N * M);
    mu.assign(M, mean_);
    for (size_t a = 0; a < N; ++a) {
      double const *xa = Xt_.data() + a * n_;
      double *row = V_.data() + a * M;
      for (size_t c = 0; c < M; ++c) {
        row[c] = kernel(xa, C + c * n_);
        mu[c] += row[c] * alpha_[a];
      }
    }
    tools::detail::forward_solve(L_.data(), N, V_.data(), M);
    // The transpose, M x N, so that the rows of the Gram matrix V^T V are contiguous dot products:
    Vt_.resize(M * N);
    for (size_t a = 0; a < N; ++a) {
      for (size_t c = 0; c < M; ++c) {
        Vt_[c * N + a] = V_[a * M + c];
      }
    }
    // The posterior covariance Kcc + noise I - V^T V, then its Cholesky factor, with jitter if it is needed:
    S.resize(M * M);
    for (size_t a = 0; a < M; ++a) {
      double const *xa = C + a * n_;
      double *row = S.data() + a * M;
      for (size_t b = 0; b < a; ++b) {
        row[b] = kernel(xa, C + b * n_);
      }
      row[a] = outputscale_ + (observation_noise ? noise_ : 0.0);
    }
    tools::detail::subtract_gram(S.data(), M, Vt_.data(), N, N, 0, M);
    // With the noise, it bounds the covariance below, so jitter should never be needed; if it is, retry from a copy.
    covariance_ = S;
    for (size_t a = 0; a < M && jitter > 0; ++a) {
      S[a * M + a] += jitter;
    }
    bool factored = tools::detail::cholesky(S.data(), M);
    jitter_ = jitter;
    for (jitter = (std::max)(1e-8, 10 * jitter); !factored && jitter < 1; jitter *= 10) {
      S = covariance_;
      for (size_t a = 0; a < M; ++a) {
        S[a * M + a] += jitter;
      }
      factored = tools::detail::cholesky(S.data(), M);
      jitter_ = jitter;
    }
    return factored;
  }

  // The kernel of the fitted model, outputscale k(x, y), and if c is given the factor *c of its gradient:
  // d/dx_i = -(*c) (x_i - y_i) inverse_squared_lengthscales()[i]. So the gradient of sum_a k(x, X_a) w_a is
  // -inverse_squared_lengthscales()[i] sum_a c_a w_a (x_i - X_ai).
  double kernel(double const *x, double const *y, double *c = nullptr) const {
    double const k = outputscale_ * detail::kernel(params_.kernel, x, y, inv_l2_.data(), n_, c);
    if (c) {
      *c *= outputscale_;
    }
    return k;
  }

  // Registers the conditioning set, the m points B (m x n), m >= 0, for the conditional variance and sample paths
  // below, until the next fit: it stores V_B = L^{-1} K(X, B), the posterior mean at B, and the Cholesky factor
  // L_B L_B^T = Sigma(B, B) + jitter I. Where B has training points, Sigma(B, B) is as ill conditioned as K, so the
  // jitter is at least 1e-10 outputscale, and ten times more until the factorization succeeds, from 1e-8 as in
  // posterior(), which this calls. Returns false if it does not; there is then no conditioning set.
  bool condition(double const *B, size_t m) {
    B_.assign(B, B + m * n_);
    m_ = posterior(B, m, mu_B_, L_B_, false, 1e-10 * outputscale_) ? m : 0;
    V_B_.swap(V_);
    Vt_B_.swap(Vt_);
    jitter_B_ = jitter_;
    return m_ == m;
  }

  // Adds the point b (n numbers) to the conditioning set. What it leaves is what condition() does for the larger set,
  // up to rounding, as a row of a Cholesky factor depends only on the rows before it; but only the new column of V_B
  // and the new row of L_B are computed, in O(N^2 + N m + m^2 + (N + m) n). If there is no set to add to, or the new
  // row cannot be factored with the jitter of the others, it is condition() that is called.
  bool append(double const *b) {
    size_t const N = N_, m = m_, M = m + 1;
    B_.resize(m * n_);
    B_.insert(B_.end(), b, b + n_);
    // A column more for V_B (N x m) and L_B (m x m), from the last row, so that nothing is overwritten before it is
    // moved; and a row more for L_B and for the transpose of V_B.
    auto widen = [&](std::vector<double> &A, size_t rows, size_t rows_after) {
      A.resize(rows_after * M);
      for (size_t r = rows; r-- > 1;) {
        std::copy_backward(A.data() + r * m, A.data() + (r + 1) * m, A.data() + r * M + m);
      }
    };
    widen(V_B_, N, N);
    widen(L_B_, m, M);
    Vt_B_.resize(M * N);
    double *v = Vt_B_.data() + m * N, *row = L_B_.data() + m * M;
    mu_B_.push_back(mean_);
    for (size_t a = 0; a < N; ++a) {
      v[a] = kernel(Xt_.data() + a * n_, b);
      mu_B_[m] += v[a] * alpha_[a];
    }
    tools::detail::forward_solve(L_.data(), N, v, 1);
    for (size_t a = 0; a < N; ++a) {
      V_B_[a * M + m] = v[a];
    }
    for (size_t j = 0; j < m; ++j) {
      row[j] = kernel(b, B_.data() + j * n_);
    }
    row[m] = outputscale_;
    tools::detail::subtract_gram(L_B_.data(), M, Vt_B_.data(), N, N, 0, M, m);
    row[m] += jitter_B_;
    m_ = m > 0 && tools::detail::cholesky(L_B_.data(), M, m) ? M : 0;
    refactored_B_ = m_ != M;
    return !refactored_B_ || condition(std::vector<double>(B_).data(), M);
  }

  // The posterior variance of the function at x given its values at the conditioning set, and if grad is not null
  // its gradient with respect to x (n numbers):
  //   d(x)^2 = Sigma(x, x) - Sigma(x, B) (Sigma(B, B) + jitter I)^{-1} Sigma(B, x) = outputscale - |v|^2 - |u|^2,
  // where v = L^{-1} k(X, x) and u = L_B^{-1} c with c = Sigma(B, x) = k(B, x) - V_B^T v, and
  //   grad d^2 = -2 (dk(X, x)/dx)^T L^{-T} (v - V_B gamma) - 2 (dk(B, x)/dx)^T gamma, gamma = L_B^{-T} u.
  // With no conditioning set it is the posterior variance Sigma(x, x). It is close to 0 (the jitter) at the points
  // of B, and never less than variance_floor(): below that it is the floor, with zero gradient. The cost is
  // O(N^2 + N m + m^2 + (N + m) n); work is scratch that is resized to 2 (N + m) only if it is smaller, so a call
  // does not allocate after the first, and calls with separate scratch may run concurrently. A caller that has them
  // gives the kernel between x and the N points of X, then the m of B, and then its N + m gradient factors, as k.
  double conditional_variance(double const *x, double *grad, std::vector<double> &work,
                              double const *k = nullptr) const {
    size_t const N = N_, m = m_;
    if (work.size() < 2 * (N + m)) {
      work.resize(2 * (N + m));
    }
    double *v = work.data(), *u = v + N;
    double const *c = k ? k + N + m : u + m;
    for (size_t a = 0; a < N + m; ++a) {
      v[a] = k ? k[a] : kernel(a < N ? Xt_.data() + a * n_ : B_.data() + (a - N) * n_, x, u + m + a);
    }
    tools::detail::forward_solve(L_.data(), N, v);
    for (size_t a = 0; a < N; ++a) {
      double const *row = V_B_.data() + a * m;
      for (size_t j = 0; j < m; ++j) {
        u[j] -= row[j] * v[a];
      }
    }
    tools::detail::forward_solve(L_B_.data(), m, u);
    double const d2 = outputscale_ - tools::detail::dot(v, v, N) - tools::detail::dot(u, u, m);
    bool const floored = !(d2 > variance_floor());
    if (grad) {
      std::fill(grad, grad + n_, 0.0);
    }
    if (grad && !floored) {
      tools::detail::backward_solve(L_B_.data(), m, u);
      for (size_t a = 0; a < N; ++a) {
        v[a] -= tools::detail::dot(V_B_.data() + a * m, u, m);
      }
      tools::detail::backward_solve(L_.data(), N, v);
      // d k(y, x) / d x_i = -c (x_i - y_i) / l_i^2, for y the N points of X with weights v, then those of B with u:
      for (size_t a = 0; a < N + m; ++a) {
        double const *ya = a < N ? Xt_.data() + a * n_ : B_.data() + (a - N) * n_;
        double const w = 2 * v[a] * c[a];
        for (size_t i = 0; i < n_; ++i) {
          grad[i] += w * (x[i] - ya[i]) * inv_l2_[i];
        }
      }
    }
    return floored ? variance_floor() : d2;
  }
  double variance_floor() const { return 1e-12 * outputscale_; }

  // The coefficients of a sample path of the posterior through the conditioning set: for the m standard normal
  // numbers zeta, the values of the function at B are f(B) = conditioning_mean() + L_B zeta, and given them
  //   f(x) = mean() + sum_a kernel(X_a, x) a[a] + sum_j kernel(B_j, x) w[j] + sqrt(d(x)^2) z, z standard normal,
  // where w = L_B^{-T} zeta (m numbers; Sigma(B, B)^{-1} (f(B) - mean at B)) and a = alpha - L^{-T} V_B w (N
  // numbers). zeta and w may be the same array.
  void sample_path(double const *zeta, double *a, double *w) const {
    std::copy(zeta, zeta + m_, w);
    tools::detail::backward_solve(L_B_.data(), m_, w);
    for (size_t r = 0; r < N_; ++r) {
      a[r] = tools::detail::dot(V_B_.data() + r * m_, w, m_);
    }
    tools::detail::backward_solve(L_.data(), N_, a);
    for (size_t r = 0; r < N_; ++r) {
      a[r] = alpha_[r] - a[r];
    }
  }

  size_t dimension() const { return n_; }
  size_t size() const { return N_; }
  std::vector<double> const &lengthscales() const { return lengthscale_; }
  std::vector<double> const &inverse_squared_lengthscales() const { return inv_l2_; }
  double outputscale() const { return outputscale_; }
  double noise() const { return noise_; }
  double mean() const { return mean_; }
  double loss() const { return loss_; }
  // Whether the last append() had to factor the whole conditioning set again, which changes every row of its factor.
  bool conditioning_refactored() const { return refactored_B_; }
  std::vector<double> const &hyperparameters() const { return theta_; }
  // The training points (N x n), the lower triangle of L (N x N; the upper is not meaningful), and alpha.
  std::vector<double> const &training_points() const { return Xt_; }
  std::vector<double> const &alpha() const { return alpha_; }
  // The conditioning set: its size, points (m x n), posterior mean, V_B (N x m), and the lower triangle of L_B (m x m).
  size_t conditioning_size() const { return m_; }
  std::vector<double> const &conditioning_points() const { return B_; }
  std::vector<double> const &conditioning_mean() const { return mu_B_; }
  std::vector<double> const &conditioning_factor() const { return L_B_; }

private:
  static double sigmoid(double t) { return 1 / (1 + std::exp(-t)); }
  static double logit(double p) { return std::log(p / (1 - p)); }
  // A hyperparameter as a function of the variable of the optimizer, its inverse, and its derivative.
  double from_raw(double raw, double const *bounds) const {
    return params_.lbfgs ? (linear_ ? raw : std::exp(raw)) : bounds[0] + (bounds[1] - bounds[0]) * sigmoid(raw);
  }
  double to_raw(double value, double const *bounds) const {
    return params_.lbfgs ? (linear_ ? value : std::log(value))
                         : logit((value - bounds[0]) / (bounds[1] - bounds[0]));
  }
  double dvalue(double raw, double const *bounds) const {
    double const sg = sigmoid(raw);
    return params_.lbfgs ? (linear_ ? 1.0 : std::exp(raw)) : (bounds[1] - bounds[0]) * sg * (1 - sg);
  }
  // Adds to the loss the negative logarithm of the LogNormal(prior[0], prior[1]) density at value, over N and up to
  // a constant, log(value) + (log(value) - loc)^2 / (2 scale^2), and returns its derivative with respect to value.
  double add_prior(double value, double const *prior, double &loss) const {
    if (!(prior[1] > 0)) {
      return 0;
    }
    double const N = static_cast<double>(N_), t = std::log(value) - prior[0], s2 = prior[1] * prior[1];
    loss += (std::log(value) + t * t / (2 * s2)) / N;
    return (1 + t / s2) / (value * N);
  }

  // Sets the hyperparameters from theta = (raw lengthscales, raw outputscale, raw noise, mean), builds the kernel
  // matrix K + noise I in L_, and factors it. pair_m_ and pair_dm_ receive the kernel's radial factor and the factor
  // of its lengthscale derivatives at each pair of training points. Returns false if the factorization fails.
  bool factor(std::vector<double> const &theta) {
    using std::exp;
    using std::sqrt;
    size_t const N = N_;
    for (size_t i = 0; i < n_; ++i) {
      lengthscale_[i] = from_raw(theta[i], params_.lengthscale_bounds);
      inv_l2_[i] = 1 / (lengthscale_[i] * lengthscale_[i]);
    }
    outputscale_ = from_raw(theta[n_], params_.outputscale_bounds);
    noise_ = from_raw(theta[n_ + 1], params_.noise_bounds);
    mean_ = theta[n_ + 2];
    L_.resize(N * N);
    pair_m_.resize(N * (N - 1) / 2);
    pair_dm_.resize(pair_m_.size());
    double const root5 = sqrt(5.0);
    bool const rbf = params_.kernel == gaussian_process_kernel::rbf;
    size_t p = 0;
    // The scaled squared distances, eight pairs at a time, so that their sums over the dimensions overlap:
    for (; p < pair_m_.size(); p += 8) {
      double r2[8] = {0, 0, 0, 0, 0, 0, 0, 0};
      size_t const count = (std::min)(size_t(8), pair_m_.size() - p);
      double const *D = sqdiff_.data() + p * n_;
      for (size_t i = 0; i < n_; ++i) {
        for (size_t j = 0; j < count; ++j) {
          r2[j] += D[j * n_ + i] * inv_l2_[i];
        }
      }
      std::copy_n(r2, count, pair_m_.data() + p);
    }
    p = 0;
    for (size_t a = 0; a < N; ++a) {
      double *row = L_.data() + a * N;
      for (size_t b = 0; b < a; ++b, ++p) {
        double const r2 = pair_m_[p];
        // The Gaussian kernel exp(-rho^2/2) is these formulas with r = 0, but for the constant 5/3.
        double const r = rbf ? 0.0 : root5 * sqrt(r2);
        double const e = exp(rbf ? -r2 / 2 : -r);
        pair_m_[p] = (1 + r + r * r / 3) * e;
        // d k / d l_i = s (5/3) (1 + sqrt(5) rho) exp(-sqrt(5) rho) (x_i - y_i)^2 / l_i^3, with r = sqrt(5) rho; for
        // the Gaussian kernel, k (x_i - y_i)^2 / l_i^3:
        pair_dm_[p] = outputscale_ * (rbf ? 1.0 : 5.0 / 3) * (1 + r) * e;
        row[b] = outputscale_ * pair_m_[p];
      }
      row[a] = outputscale_ + noise_;
    }
    return tools::detail::cholesky(L_.data(), N);
  }

public:
  // The loss, -(log p(y | theta) + log p(theta))/N (as GPyTorch's ExactMarginalLogLikelihood with its priors, up to
  // a constant), and its gradient with respect to theta. Leaves L_ and alpha_ = (K + noise I)^{-1} (y - mean). Returns
  // false if the kernel matrix cannot be factored. Public for the tests, as it leaves the model at theta, which is
  // what the optimizer sees and hyperparameters() returns after a fit: for each lengthscale, the outputscale and
  // the noise, the logit of its position between its bounds, or with L-BFGS its logarithm; then the mean.
  bool loss_and_gradient(std::vector<double> const &theta, double &loss, std::vector<double> &grad) {
    using std::log;
    size_t const N = N_;
    m_ = 0;
    if (!factor(theta)) {
      return false;
    }
    alpha_.resize(N);
    for (size_t a = 0; a < N; ++a) {
      alpha_[a] = yt_[a] - mean_;
    }
    tools::detail::cholesky_solve(L_.data(), N, alpha_.data());
    loss = 0;
    for (size_t a = 0; a < N; ++a) {
      loss += (yt_[a] - mean_) * alpha_[a] / 2 + log(L_[a * N + a]);
    }
    loss /= static_cast<double>(N);
    Z_.resize(N * N);
    Kinv_.resize(N * N);
    tools::detail::cholesky_inverse(L_.data(), N, Z_.data(), Kinv_.data());
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
        pair_dm_[p] *= w;
      }
      double const w = alpha_[a] * alpha_[a] - kinv[a];
      trace_noise += w;
      trace_scale += w;
    }
    // Four pairs at a time, as each of these sums is a chain of additions through memory:
    for (p = 0; p + 3 < pair_dm_.size(); p += 4) {
      double const *D = sqdiff_.data() + p * n_, *w = pair_dm_.data() + p;
      for (size_t i = 0; i < n_; ++i) {
        trace_l[i] += (w[0] * D[i] + w[1] * D[n_ + i]) + (w[2] * D[2 * n_ + i] + w[3] * D[3 * n_ + i]);
      }
    }
    for (; p < pair_dm_.size(); ++p) {
      for (size_t i = 0; i < n_; ++i) {
        trace_l[i] += pair_dm_[p] * sqdiff_[p * n_ + i];
      }
    }
    double sum_alpha = 0;
    for (size_t a = 0; a < N; ++a) {
      sum_alpha += alpha_[a];
    }
    double const scale = -1 / (2 * static_cast<double>(N));
    // The priors, and the chain rule through the value as a function of theta:
    for (size_t i = 0; i < n_; ++i) {
      double const l = lengthscale_[i];
      double const prior = add_prior(l, params_.lengthscale_prior, loss);
      grad[i] = (scale * trace_l[i] / (l * l * l) + prior) * dvalue(theta[i], params_.lengthscale_bounds);
    }
    grad[n_] = scale * trace_scale * dvalue(theta[n_], params_.outputscale_bounds);
    double const prior = add_prior(noise_, params_.noise_prior, loss);
    grad[n_ + 1] = (scale * trace_noise + prior) * dvalue(theta[n_ + 1], params_.noise_bounds);
    grad[n_ + 2] = -sum_alpha / static_cast<double>(N);
    return true;
  }

private:
  gaussian_process_parameters params_;
  // The training set, hyperparameters, and scratch.
  size_t n_ = 0;
  size_t N_ = 0;
  std::vector<double> Xt_, yt_;
  std::vector<double> lengthscale_, inv_l2_;
  double outputscale_ = 1, noise_ = 0.005, mean_ = 0, loss_ = 0;
  bool refactored_B_ = false;
  std::vector<double> sqdiff_, pair_m_, pair_dm_;
  std::vector<double> L_, Z_, Kinv_, alpha_;
  std::vector<double> V_, Vt_, covariance_;
  // The jitter with which the last posterior was factored, and that of the conditioning set.
  double jitter_ = 0, jitter_B_ = 0;
  // The hyperparameters of the last fit as the optimizer has them, and the conditioning set.
  std::vector<double> theta_;
  size_t m_ = 0;
  std::vector<double> B_, mu_B_, V_B_, Vt_B_, L_B_;
  // Whether the variables of the L-BFGS fit are, for now, the hyperparameters rather than their logarithms.
  bool linear_ = false;
};

} BOOST_MATH_NAMESPACE_END // namespace boost::math::optimization::detail
#endif
