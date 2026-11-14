/*
 * Copyright Nick Thompson, 2026
 * Use, modification and distribution are subject to the
 * Boost Software License, Version 1.0. (See accompanying file
 * LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
 */

#ifndef BOOST_MATH_BUILD_MODULE
#include <boost/math/optimization/detail/gaussian_process.hpp>
#else
import boost.math;
#endif

#include "math_unit_test.hpp"
#ifndef BOOST_MATH_BUILD_MODULE
#include <algorithm>
#include <cmath>
#include <limits>
#include <random>
#include <vector>
#endif
using boost::math::optimization::detail::cholesky;
using boost::math::optimization::detail::cholesky_solve;
using boost::math::optimization::detail::forward_solve;
using boost::math::optimization::detail::gaussian_process;
using boost::math::optimization::detail::gaussian_process_kernel;
using boost::math::optimization::detail::gaussian_process_parameters;
using boost::math::optimization::detail::kernel;

void test_cholesky() {
  // The upper triangle is not referenced:
  double const nan = std::numeric_limits<double>::quiet_NaN();
  std::vector<double> A{4, nan, nan, 12, 37, nan, -16, -43, 98};
  CHECK_TRUE(cholesky(A.data(), 3));
  double const L[9] = {2, 0, 0, 6, 1, 0, -8, 5, 3};
  for (size_t i = 0; i < 3; ++i) {
    for (size_t j = 0; j <= i; ++j) {
      CHECK_ULP_CLOSE(L[i * 3 + j], A[i * 3 + j], 2);
    }
  }
  // A x = (1, 2, 3):
  std::vector<double> x{1, 2, 3};
  cholesky_solve(A.data(), 3, x.data());
  CHECK_ULP_CLOSE(343.0 / 12, x[0], 100);
  CHECK_ULP_CLOSE(-23.0 / 3, x[1], 100);
  CHECK_ULP_CLOSE(4.0 / 3, x[2], 100);
  // L X = B for two right hand sides, the columns of B:
  std::vector<double> B{2, 4, 7, 14, -3, 6};
  forward_solve(A.data(), 3, B.data(), 2);
  double const X[6] = {1, 2, 1, 2, 0, 4};
  for (size_t i = 0; i < 6; ++i) {
    CHECK_ULP_CLOSE(X[i], B[i], 4);
  }

  std::vector<double> indefinite{1, 0, 2, 1};
  CHECK_TRUE(!cholesky(indefinite.data(), 2));

  // Larger than the blocks: A = G G^T + n I.
  size_t const n = 150;
  std::mt19937_64 gen(42);
  std::uniform_real_distribution<double> dis(-1, 1);
  std::vector<double> G(n * n);
  for (auto &g : G) {
    g = dis(gen);
  }
  std::vector<double> S(n * n);
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < n; ++j) {
      double s = i == j ? static_cast<double>(n) : 0.0;
      for (size_t k = 0; k < n; ++k) {
        s += G[i * n + k] * G[j * n + k];
      }
      S[i * n + j] = s;
    }
  }
  std::vector<double> F = S;
  CHECK_TRUE(cholesky(F.data(), n));
  double worst = 0;
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j <= i; ++j) {
      double s = 0;
      for (size_t k = 0; k <= j; ++k) {
        s += F[i * n + k] * F[j * n + k];
      }
      worst = (std::max)(worst, std::abs(s - S[i * n + j]));
    }
  }
  CHECK_LE(worst, 1e-11);
}

void test_posterior() {
  // A smooth function of two variables, sampled at 40 random points of the unit square:
  auto f = [](double const *x) { return std::sin(3 * x[0]) + x[1] * x[1] - 1; };
  size_t const N = 40;
  std::mt19937_64 gen(1234);
  std::uniform_real_distribution<double> dis(0, 1);
  std::vector<double> X(2 * N), y(N);
  for (auto &x : X) {
    x = dis(gen);
  }
  for (size_t k = 0; k < N; ++k) {
    y[k] = f(X.data() + 2 * k);
  }
  gaussian_process gp;
  CHECK_TRUE(gp.fit(X.data(), y.data(), N, 2));
  CHECK_EQUAL(gp.dimension(), size_t(2));
  CHECK_EQUAL(gp.size(), N);
  CHECK_EQUAL(gp.lengthscales().size(), size_t(2));
  for (double l : gp.lengthscales()) {
    CHECK_TRUE(0.005 < l && l < 2);
  }
  CHECK_TRUE(0.05 < gp.outputscale() && gp.outputscale() < 20);
  CHECK_TRUE(0.0005 < gp.noise() && gp.noise() < 0.2);

  // At the training points, the posterior of the function is the data, up to the (small) fitted noise:
  std::vector<double> mu, S;
  CHECK_TRUE(gp.posterior(X.data(), N, mu, S, false));
  for (size_t k = 0; k < N; ++k) {
    CHECK_ABSOLUTE_ERROR(y[k], mu[k], 0.02);
    double variance = 0;
    for (size_t j = 0; j <= k; ++j) {
      variance += S[k * N + j] * S[k * N + j];
    }
    CHECK_LE(variance, gp.noise());
  }

  // Elsewhere it predicts the function, the more uncertain the further from the data:
  size_t const M = 25;
  std::vector<double> C(2 * M);
  for (auto &c : C) {
    c = dis(gen);
  }
  CHECK_TRUE(gp.posterior(C.data(), M, mu, S, false));
  for (size_t k = 0; k < M; ++k) {
    CHECK_ABSOLUTE_ERROR(f(C.data() + 2 * k), mu[k], 0.05);
  }
  // The observations have the same mean, and their covariance is that of the function plus the noise on the diagonal.
  std::vector<double> mu_noisy, S_noisy;
  CHECK_TRUE(gp.posterior(C.data(), M, mu_noisy, S_noisy, true));
  for (size_t a = 0; a < M; ++a) {
    CHECK_ULP_CLOSE(mu[a], mu_noisy[a], 0);
    for (size_t b = 0; b <= a; ++b) {
      double c = 0, c_noisy = 0;
      for (size_t j = 0; j <= b; ++j) {
        c += S[a * M + j] * S[b * M + j];
        c_noisy += S_noisy[a * M + j] * S_noisy[b * M + j];
      }
      CHECK_ABSOLUTE_ERROR(c + (a == b ? gp.noise() : 0.0), c_noisy, 1e-12);
    }
  }
  // The marginal of the joint posterior at one of the points is the posterior at that point alone:
  std::vector<double> mu1, S1;
  CHECK_TRUE(gp.posterior(C.data() + 2 * 7, 1, mu1, S1, true));
  double variance = 0;
  for (size_t j = 0; j <= 7; ++j) {
    variance += S_noisy[7 * M + j] * S_noisy[7 * M + j];
  }
  CHECK_ABSOLUTE_ERROR(mu_noisy[7], mu1[0], 1e-12);
  CHECK_ABSOLUTE_ERROR(variance, S1[0] * S1[0], 1e-12);
  // Far from the data, the prior: the fitted mean and outputscale.
  double const far[2] = {100, -100};
  CHECK_TRUE(gp.posterior(far, 1, mu1, S1, false));
  CHECK_ABSOLUTE_ERROR(gp.mean(), mu1[0], 1e-12);
  CHECK_ABSOLUTE_ERROR(gp.outputscale(), S1[0] * S1[0], 1e-12);
}

// Both kernels: the gradient with respect to x against central differences, also next to and at y.
void test_kernel_gradient() {
  size_t const d = 5;
  std::mt19937_64 gen(5);
  std::uniform_real_distribution<double> dis(0, 1);
  for (auto kind : {gaussian_process_kernel::matern52, gaussian_process_kernel::rbf}) {
    for (double distance : {1.0, 1e-3, 1e-7, 0.0}) {
      std::vector<double> x(d), y(d), inv_l2(d);
      for (size_t i = 0; i < d; ++i) {
        y[i] = dis(gen);
        x[i] = y[i] + distance * (dis(gen) - 0.5);
        inv_l2[i] = 1 / (0.04 + dis(gen));
      }
      double c;
      double const k = kernel(kind, x.data(), y.data(), inv_l2.data(), d, &c);
      CHECK_ULP_CLOSE(k, kernel(kind, x.data(), y.data(), inv_l2.data(), d), 0);
      CHECK_ULP_CLOSE(k, kernel(kind, y.data(), x.data(), inv_l2.data(), d), 0);
      CHECK_TRUE(std::isfinite(c));
      if (distance == 0) {
        CHECK_ULP_CLOSE(1.0, k, 0);
      }
      for (size_t i = 0; i < d; ++i) {
        double const h = 1e-6, xi = x[i];
        x[i] = xi + h;
        double const up = kernel(kind, x.data(), y.data(), inv_l2.data(), d);
        x[i] = xi - h;
        double const down = kernel(kind, x.data(), y.data(), inv_l2.data(), d);
        x[i] = xi;
        // Where the kernel is at its maximum, the difference quotient is 0 up to rounding, and so is the gradient.
        CHECK_ABSOLUTE_ERROR((up - down) / (2 * h), -c * (x[i] - y[i]) * inv_l2[i], 2e-9);
      }
    }
  }
  // r^2 = 1/4 * 4 + 1/16 * 16 = 2, and sqrt(5) r = sqrt(10) for the Matern kernel.
  double const x[2] = {0.25, 0.5};
  double const y[2] = {0.75, 0.25};
  double const inv_l2[2] = {4, 16};
  CHECK_ULP_CLOSE(std::exp(-1.0), kernel(gaussian_process_kernel::rbf, x, y, inv_l2, 2), 2);
  CHECK_ULP_CLOSE((1 + std::sqrt(10.0) + 10.0 / 3) * std::exp(-std::sqrt(10.0)),
                  kernel(gaussian_process_kernel::matern52, x, y, inv_l2, 2), 4);
}

// N training points in n dimensions and values that depend on the first two only, standardized.
void training_set(size_t N, size_t n, std::vector<double> &X, std::vector<double> &y, unsigned seed) {
  std::mt19937_64 gen(seed);
  std::uniform_real_distribution<double> dis(0, 1);
  X.resize(N * n);
  y.resize(N);
  for (auto &x : X) {
    x = dis(gen);
  }
  double mean = 0, variance = 0;
  for (size_t k = 0; k < N; ++k) {
    y[k] = std::sin(5 * X[k * n]) + 2 * (X[k * n + 1] - 0.3) * (X[k * n + 1] - 0.3) + 0.01 * (dis(gen) - 0.5);
    mean += y[k] / static_cast<double>(N);
  }
  for (size_t k = 0; k < N; ++k) {
    variance += (y[k] - mean) * (y[k] - mean) / static_cast<double>(N - 1);
  }
  for (size_t k = 0; k < N; ++k) {
    y[k] = (y[k] - mean) / std::sqrt(variance);
  }
}

// The gradient of the loss against central differences: of the marginal likelihood with Adam's sigmoids and the
// Matern kernel, and of the posterior with logarithms, the Gaussian kernel and the priors.
void test_loss_gradient() {
  size_t const N = 25, n = 4;
  std::vector<double> X, y;
  training_set(N, n, X, y, 11);
  std::mt19937_64 gen(12);
  std::uniform_real_distribution<double> dis(-0.5, 0.5);
  for (bool botorch : {false, true}) {
    auto params = botorch ? gaussian_process_parameters::botorch(n) : gaussian_process_parameters();
    params.hyperparameter_steps = 3;
    gaussian_process gp(params);
    CHECK_TRUE(gp.fit(X.data(), y.data(), N, n));
    std::vector<double> theta = gp.hyperparameters();
    CHECK_EQUAL(n + 3, theta.size());
    for (auto &t : theta) {
      t += dis(gen);
    }
    std::vector<double> grad(n + 3), scratch(n + 3);
    double loss, up, down;
    CHECK_TRUE(gp.loss_and_gradient(theta, loss, grad));
    for (size_t k = 0; k < n + 3; ++k) {
      double const h = 1e-5, t = theta[k];
      theta[k] = t + h;
      CHECK_TRUE(gp.loss_and_gradient(theta, up, scratch));
      theta[k] = t - h;
      CHECK_TRUE(gp.loss_and_gradient(theta, down, scratch));
      theta[k] = t;
      CHECK_ABSOLUTE_ERROR((up - down) / (2 * h), grad[k], 1e-7);
    }
  }
}

// The fit of BoTorch's model to a function of 2 of 10 variables: long lengthscales on the others.
void test_map_fit() {
  size_t const N = 60, n = 10;
  std::vector<double> X, y;
  training_set(N, n, X, y, 21);
  gaussian_process_parameters params = gaussian_process_parameters::botorch(n);
  CHECK_ABSOLUTE_ERROR(std::exp(std::sqrt(2.0) + std::log(10.0) / 2 - 3), params.initial_lengthscale, 1e-15);
  gaussian_process gp(params);
  CHECK_TRUE(gp.fit(X.data(), y.data(), N, n));
  CHECK_EQUAL(1.0, gp.outputscale());
  CHECK_LE(1e-4, gp.noise());
  CHECK_LE(gp.noise(), 1e-2);
  for (size_t i = 0; i < n; ++i) {
    if (i < 2) {
      CHECK_LE(0.1, gp.lengthscales()[i]);
      CHECK_LE(gp.lengthscales()[i], 1.5);
    }
    else {
      CHECK_LE(5.0, gp.lengthscales()[i]);
    }
  }
  // The fit is a stationary point of the loss in the hyperparameters that are not on a bound:
  std::vector<double> theta = gp.hyperparameters(), grad(n + 3);
  double loss;
  CHECK_TRUE(gp.loss_and_gradient(theta, loss, grad));
  for (size_t k = 0; k < n + 3; ++k) {
    if (k < n || k == n + 2 || (k == n + 1 && gp.noise() > 1.001e-4)) {
      CHECK_ABSOLUTE_ERROR(0.0, grad[k], 5e-3);
    }
  }
  CHECK_LE(0.0, grad[n + 1]);
  // and predicts the function elsewhere, from training data other than its own:
  std::vector<double> C, f, mu, S;
  training_set(N, n, C, f, 21);
  std::mt19937_64 gen(22);
  std::uniform_real_distribution<double> dis(0, 1);
  for (size_t k = 0; k < N; ++k) {
    for (size_t i = 2; i < n; ++i) {
      C[k * n + i] = dis(gen);
    }
  }
  CHECK_TRUE(gp.posterior(C.data(), N, mu, S, false));
  for (size_t k = 0; k < N; ++k) {
    CHECK_ABSOLUTE_ERROR(f[k], mu[k], 0.05);
  }
}

// A fit that starts from the hyperparameters of the fit to the first forty points ends where the fit from the
// modes of the priors does; with nothing to start from, it is that fit.
void test_warm_fit() {
  size_t const N = 60, n = 10;
  std::vector<double> X, y;
  training_set(N, n, X, y, 21);
  gaussian_process cold(gaussian_process_parameters::botorch(n)), warm = cold, first = cold;
  CHECK_TRUE(cold.fit(X.data(), y.data(), N, n));
  CHECK_TRUE(first.fit(X.data(), y.data(), N, n, true));
  CHECK_TRUE(warm.fit(X.data(), y.data(), N - 20, n));
  CHECK_LE(0.01, std::abs(warm.lengthscales()[0] - cold.lengthscales()[0]));
  CHECK_TRUE(warm.fit(X.data(), y.data(), N, n, true));
  std::vector<double> grad(n + 3);
  double loss_cold, loss_warm;
  CHECK_TRUE(cold.loss_and_gradient(cold.hyperparameters(), loss_cold, grad));
  CHECK_TRUE(warm.loss_and_gradient(warm.hyperparameters(), loss_warm, grad));
  CHECK_ABSOLUTE_ERROR(loss_cold, loss_warm, 1e-3);
  for (size_t i = 0; i < n; ++i) {
    CHECK_EQUAL(cold.lengthscales()[i], first.lengthscales()[i]);
  }
  // The two variables that the function depends on; the loss is flat in the long lengthscales of the others.
  for (size_t i = 0; i < 2; ++i) {
    CHECK_ABSOLUTE_ERROR(cold.lengthscales()[i], warm.lengthscales()[i], 0.1 * cold.lengthscales()[i]);
  }
}

// The variance given the values at a conditioning set, its gradient, and the sample paths.
void test_conditional_variance() {
  size_t const N = 30, n = 3, m = 12;
  std::vector<double> X, y;
  training_set(N, n, X, y, 31);
  std::mt19937_64 gen(32);
  std::uniform_real_distribution<double> dis(0, 1);
  std::normal_distribution<double> normal;
  for (bool botorch : {false, true}) {
    gaussian_process gp(botorch ? gaussian_process_parameters::botorch(n) : gaussian_process_parameters());
    CHECK_TRUE(gp.fit(X.data(), y.data(), N, n));
    double const s = gp.outputscale();
    std::vector<double> work, grad(n), mu, S, x(n);
    // Without a conditioning set, the posterior variance:
    CHECK_EQUAL(size_t(0), gp.conditioning_size());
    for (int trial = 0; trial < 5; ++trial) {
      for (auto &v : x) {
        v = dis(gen);
      }
      CHECK_TRUE(gp.posterior(x.data(), 1, mu, S, false));
      CHECK_ABSOLUTE_ERROR(S[0] * S[0], gp.conditional_variance(x.data(), nullptr, work), 1e-12 * s);
    }
    // The first 8 training points and 4 others:
    std::vector<double> B(X.begin(), X.begin() + m * n);
    for (size_t k = 8 * n; k < m * n; ++k) {
      B[k] = dis(gen);
    }
    CHECK_TRUE(gp.condition(B.data(), m));
    CHECK_EQUAL(m, gp.conditioning_size());
    for (size_t j = 0; j < m; ++j) {
      double const d2 = gp.conditional_variance(B.data() + j * n, grad.data(), work);
      CHECK_LE(gp.variance_floor(), d2);
      CHECK_LE(d2, 1e-8 * s);
    }
    // At random points, 1e-3 and 1e-6 from a training point that is not in B, and 1e-2 and 1e-4 from points of B:
    for (int trial = 0; trial < 12; ++trial) {
      double const *center = trial % 4 == 1 ? X.data() + (N - 1) * n : B.data() + (trial % 4 == 2 ? 2 : 9) * n;
      double const radius = trial % 4 == 0 ? 1.0 : trial % 4 == 1 ? (trial < 4 ? 1e-3 : 1e-6) : trial < 4 ? 1e-2 : 1e-4;
      for (size_t i = 0; i < n; ++i) {
        x[i] = trial % 4 == 0 ? dis(gen) : center[i] + radius * (dis(gen) - 0.5);
      }
      double const d2 = gp.conditional_variance(x.data(), grad.data(), work);
      CHECK_ULP_CLOSE(d2, gp.conditional_variance(x.data(), nullptr, work), 0);
      CHECK_LE(d2, s);
      // It is the last pivot of the Cholesky factorization of the joint posterior covariance at B and x:
      std::vector<double> joint = B;
      joint.insert(joint.end(), x.begin(), x.end());
      CHECK_TRUE(gp.posterior(joint.data(), m + 1, mu, S, false, 1e-10 * s));
      CHECK_ABSOLUTE_ERROR(S[m * (m + 1) + m] * S[m * (m + 1) + m] - 1e-10 * s, d2, 1e-11 * s);
      for (size_t i = 0; i < n; ++i) {
        double const h = 1e-5, xi = x[i];
        x[i] = xi + h;
        double const up = gp.conditional_variance(x.data(), nullptr, work);
        x[i] = xi - h;
        double const down = gp.conditional_variance(x.data(), nullptr, work);
        x[i] = xi;
        CHECK_ABSOLUTE_ERROR((up - down) / (2 * h), grad[i], 1e-7 * s);
      }
    }
    // A sample path takes its values at B, and has the coefficients of the posterior mean if those are the mean:
    std::vector<double> zeta(m), a(N), w(m);
    for (auto &z : zeta) {
      z = normal(gen);
    }
    gp.sample_path(zeta.data(), a.data(), w.data());
    for (size_t j = 0; j < m; ++j) {
      double const *b = B.data() + j * n;
      double path = gp.mean(), value = gp.conditioning_mean()[j];
      for (size_t k = 0; k < N; ++k) {
        path += gp.kernel(X.data() + k * n, b) * a[k];
      }
      for (size_t k = 0; k < m; ++k) {
        path += gp.kernel(B.data() + k * n, b) * w[k];
        value += k <= j ? gp.conditioning_factor()[j * m + k] * zeta[k] : 0.0;
      }
      CHECK_ABSOLUTE_ERROR(value, path, 1e-6);
    }
    std::fill(zeta.begin(), zeta.end(), 0.0);
    gp.sample_path(zeta.data(), a.data(), zeta.data());
    for (size_t k = 0; k < N; ++k) {
      CHECK_ULP_CLOSE(gp.alpha()[k], a[k], 0);
    }
    // Refitting forgets the conditioning set.
    CHECK_TRUE(gp.fit(X.data(), y.data(), N, n));
    CHECK_EQUAL(size_t(0), gp.conditioning_size());
  }
}

// The numbers are those of BoTorch 0.18.1 (GPyTorch 1.15.2) on the same data in double precision: SingleTaskGP(X, y,
// outcome_transform=None) fitted by fit_gpytorch_mll, then model.posterior(T). Its optimizer differs, so the
// agreement is that of two optimizers on a loss that is flat in the mean and the long lengthscales.
void test_botorch() {
  size_t const N = 24, n = 4;
  double const alpha[4] = {std::sqrt(2.0), std::sqrt(3.0), std::sqrt(5.0), std::sqrt(7.0)};
  std::vector<double> X(N * n), y(N);
  double mean = 0, variance = 0;
  for (size_t k = 0; k < N; ++k) {
    double *x = X.data() + k * n;
    for (size_t i = 0; i < n; ++i) {
      x[i] = std::fmod(static_cast<double>(k + 1) * alpha[i], 1.0);
    }
    y[k] = std::sin(4 * x[0]) + 2 * (x[1] - 0.5) * (x[1] - 0.5) + 0.3 * x[2] +
           0.05 * std::sin(37.0 * static_cast<double>(k + 1));
    mean += y[k] / static_cast<double>(N);
  }
  for (size_t k = 0; k < N; ++k) {
    variance += (y[k] - mean) * (y[k] - mean) / static_cast<double>(N - 1);
  }
  for (size_t k = 0; k < N; ++k) {
    y[k] = (y[k] - mean) / std::sqrt(variance);
  }
  gaussian_process gp(gaussian_process_parameters::botorch(n));
  CHECK_TRUE(gp.fit(X.data(), y.data(), N, n));
  double const lengthscale[4] = {0.3631258728635271, 0.60197594169502655, 2.5735689593883944, 10.710054303708992};
  for (size_t i = 0; i < n; ++i) {
    CHECK_RELATIVE_ERROR(lengthscale[i], gp.lengthscales()[i], 1e-3);
  }
  CHECK_RELATIVE_ERROR(0.0027770526074676338, gp.noise(), 1e-3);
  CHECK_ABSOLUTE_ERROR(0.1397198889046373, gp.mean(), 1e-3);
  std::vector<double> T{0.5, 0.5, 0.5, 0.5, 0.1, 0.9, 0.3, 0.7, 0.95, 0.05, 0.6, 0.2, 0.3, 0.3, 0.9, 0.9};
  T.insert(T.end(), X.begin() + 5 * n, X.begin() + 6 * n);
  double const mu_T[5] = {0.68327163410366232, 0.33612466950248387, -1.0649233598350545, 0.97457308919004448,
                          0.70552295539744792};
  double const variance_T[5] = {0.0021259511137416975, 0.0079003458058752862, 0.010580162658763936,
                                0.0043979419784985208, 0.001652996562215181};
  std::vector<double> mu, S;
  CHECK_TRUE(gp.posterior(T.data(), 5, mu, S, false));
  for (size_t k = 0; k < 5; ++k) {
    double v = 0;
    for (size_t j = 0; j <= k; ++j) {
      v += S[k * 5 + j] * S[k * 5 + j];
    }
    CHECK_ABSOLUTE_ERROR(mu_T[k], mu[k], 1e-4);
    CHECK_RELATIVE_ERROR(variance_T[k], v, 1e-3);
  }
  CHECK_ABSOLUTE_ERROR(0.00031741437943433315, S[5] * S[0], 1e-6);
}

int main() {
  test_cholesky();
  test_posterior();
  test_kernel_gradient();
  test_loss_gradient();
  test_map_fit();
  test_warm_fit();
  test_conditional_variance();
  test_botorch();
  return boost::math::test::report_errors();
}
