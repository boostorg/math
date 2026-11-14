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
#ifndef BOOST_MATH_BUILD_MODULE
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <vector>
#endif

BOOST_MATH_NAMESPACE_BEGIN namespace optimization::detail {

// A Gaussian process with a constant mean and an ARD Matern-5/2 or Gaussian kernel, for the Bayesian optimizers: the
// dense linear algebra it needs, then the model. It began as the surrogate of TuRBO (turbo.hpp), whose reference code
// (github.com/uber-research/TuRBO, gp.py) the fit and its defaults follow; the other choice is the model of BoTorch.

// Hand-rolled dense linear algebra on row-major double arrays. The kernels with long inner loops are written so that
// the compiler can vectorize them without reassociating floating point sums: either as axpy updates of contiguous
// rows, or as several rows' dot products with one shared operand, which gives independent accumulators.

// s[r] = sum_{k < m} a[r][k] b[k] for four rows a[0..3].
inline void dot4(double const *const *a, double const *b, size_t m, double *s) {
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

inline double dot(double const *a, double const *b, size_t m) {
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
inline void subtract_gram(double *C, size_t ldc, double const *X, size_t ldx, size_t len, size_t lo, size_t hi) {
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
          dot4(rows, xj, len, s);
          for (size_t b = 0; b < nb; ++b) {
            C[(i + b) * ldc + j] -= s[b];
          }
        }
        else {
          for (size_t b = j - i; b < nb; ++b) {
            C[(i + b) * ldc + j] -= dot(rows[b], xj, len);
          }
        }
      }
    }
  }
}

// The Cholesky factorization A = L L^T in place: on input the lower triangle of the n x n matrix A (leading dimension
// n), on output that of L. The upper triangle is not referenced. Returns false if A is not numerically positive
// definite. Right-looking by blocks of 64 columns; within a block, left-looking four rows at a time.
BOOST_MATH_TEST_EXPORT inline bool cholesky(double *A, size_t n) {
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
          dot4(parts, lj + k0, j - k0, s);
          for (size_t b = 0; b < nb; ++b) {
            rows[b][j] = (rows[b][j] - s[b]) / lj[j];
          }
          continue;
        }
        // Column j meets these rows on the diagonal: row j first, then the rows below it.
        for (size_t b = j - i; b < nb; ++b) {
          double const d = rows[b][j] - dot(parts[b], lj + k0, j - k0);
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
    subtract_gram(A, n, A + k0, n, k1 - k0, k1, n);
  }
  return true;
}

// Solves L X = B in place for the n x m matrix B (leading dimension m), where L is the lower triangle of an n x n
// matrix (leading dimension n).
BOOST_MATH_TEST_EXPORT inline void forward_solve(double const *L, size_t n, double *B, size_t m) {
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

// Solves L x = b in place for a vector b.
BOOST_MATH_TEST_EXPORT inline void forward_solve(double const *L, size_t n, double *b) {
  for (size_t i = 0; i < n; ++i) {
    b[i] = (b[i] - dot(L + i * n, b, i)) / L[i * n + i];
  }
}

// Solves L^T x = b in place, by rows of L: once x_i is known, remove its contribution L_ik x_i from b_k for k < i.
BOOST_MATH_TEST_EXPORT inline void backward_solve(double const *L, size_t n, double *b) {
  for (size_t i = n; i-- > 0;) {
    double const *li = L + i * n;
    double const xi = b[i] / li[i];
    b[i] = xi;
    for (size_t k = 0; k < i; ++k) {
      b[k] -= li[k] * xi;
    }
  }
}

// Solves L L^T x = b in place.
BOOST_MATH_TEST_EXPORT inline void cholesky_solve(double const *L, size_t n, double *b) {
  forward_solve(L, n, b);
  backward_solve(L, n, b);
}

// The lower triangle of (L L^T)^{-1} = Z^T Z, where Z = L^{-1}. Z is n x n scratch. Four rows of Z are taken at a
// time, for which each row of Z has zeros in the three places after its diagonal.
inline void cholesky_inverse(double const *L, size_t n, double *Z, double *Kinv) {
  for (size_t i = 0; i < n; ++i) {
    double *zi = Z + i * n;
    std::fill(zi, zi + (std::min)(n, i + 4), 0.0);
    size_t k = 0;
    for (; k + 3 < i; k += 4) {
      double const *l = L + i * n + k, *z0 = Z + k * n, *z1 = z0 + n, *z2 = z1 + n, *z3 = z2 + n;
      for (size_t j = 0; j <= k + 3; ++j) {
        zi[j] = zi[j] - l[0] * z0[j] - l[1] * z1[j] - l[2] * z2[j] - l[3] * z3[j];
      }
    }
    for (; k < i; ++k) {
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
    double *row = Kinv + a * n;
    std::fill(row, row + a + 1, 0.0);
    size_t k = a;
    for (; k + 3 < n; k += 4) {
      double const *z0 = Z + k * n, *z1 = z0 + n, *z2 = z1 + n, *z3 = z2 + n;
      double const c0 = z0[a], c1 = z1[a], c2 = z2[a], c3 = z3[a];
      for (size_t b = 0; b <= a; ++b) {
        row[b] = row[b] + c0 * z0[b] + c1 * z1[b] + c2 * z2[b] + c3 * z3[b];
      }
    }
    for (; k < n; ++k) {
      double const *zk = Z + k * n;
      double const zka = zk[a];
      for (size_t b = 0; b <= a; ++b) {
        row[b] += zka * zk[b];
      }
    }
  }
}

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
  // Leaves L_ and alpha_ for the posterior. Returns false if there is no model. Forgets the conditioning set.
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
    double loss;
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
        return loss_and_gradient(t, loss, g) ? loss : std::numeric_limits<double>::quiet_NaN();
      };
      bad = !isfinite(bounded_lbfgs(objective, theta, lower, upper, lbfgs).value);
    }
    linear_ = false;
    for (size_t step = 0; step < params_.hyperparameter_steps && !params_.lbfgs; ++step) {
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
      theta = last_good;
      if (!factor(theta)) {
        return false;
      }
    }
    theta_ = theta;
    alpha_.resize(N);
    for (size_t a = 0; a < N; ++a) {
      alpha_[a] = yt_[a] - mean_;
    }
    cholesky_solve(L_.data(), N, alpha_.data());
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
    forward_solve(L_.data(), N, V_.data(), M);
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
    subtract_gram(S.data(), M, Vt_.data(), N, N, 0, M);
    // With the noise, it bounds the covariance below, so jitter should never be needed; if it is, retry from a copy.
    covariance_ = S;
    for (size_t a = 0; a < M && jitter > 0; ++a) {
      S[a * M + a] += jitter;
    }
    bool factored = cholesky(S.data(), M);
    for (jitter = (std::max)(1e-8, 10 * jitter); !factored && jitter < 1; jitter *= 10) {
      S = covariance_;
      for (size_t a = 0; a < M; ++a) {
        S[a * M + a] += jitter;
      }
      factored = cholesky(S.data(), M);
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
    return m_ == m;
  }

  // The posterior variance of the function at x given its values at the conditioning set, and if grad is not null
  // its gradient with respect to x (n numbers):
  //   d(x)^2 = Sigma(x, x) - Sigma(x, B) (Sigma(B, B) + jitter I)^{-1} Sigma(B, x) = outputscale - |v|^2 - |u|^2,
  // where v = L^{-1} k(X, x) and u = L_B^{-1} c with c = Sigma(B, x) = k(B, x) - V_B^T v, and
  //   grad d^2 = -2 (dk(X, x)/dx)^T L^{-T} (v - V_B gamma) - 2 (dk(B, x)/dx)^T gamma, gamma = L_B^{-T} u.
  // With no conditioning set it is the posterior variance Sigma(x, x). It is close to 0 (the jitter) at the points
  // of B, and never less than variance_floor(): below that it is the floor, with zero gradient. The cost is
  // O(N^2 + N m + m^2 + (N + m) n); work is scratch that is resized to 2 (N + m) only if it is smaller, so a call
  // does not allocate after the first, and calls with separate scratch may run concurrently.
  double conditional_variance(double const *x, double *grad, std::vector<double> &work) const {
    size_t const N = N_, m = m_;
    if (work.size() < 2 * (N + m)) {
      work.resize(2 * (N + m));
    }
    // The gradient factors of the kernel are kept in cx and cb.
    double *v = work.data(), *cx = v + N, *u = cx + N, *cb = u + m;
    for (size_t a = 0; a < N; ++a) {
      v[a] = kernel(Xt_.data() + a * n_, x, cx + a);
    }
    forward_solve(L_.data(), N, v);
    for (size_t j = 0; j < m; ++j) {
      u[j] = kernel(B_.data() + j * n_, x, cb + j);
    }
    for (size_t a = 0; a < N; ++a) {
      double const *row = V_B_.data() + a * m;
      for (size_t j = 0; j < m; ++j) {
        u[j] -= row[j] * v[a];
      }
    }
    forward_solve(L_B_.data(), m, u);
    double const d2 = outputscale_ - dot(v, v, N) - dot(u, u, m);
    bool const floored = !(d2 > variance_floor());
    if (grad) {
      std::fill(grad, grad + n_, 0.0);
    }
    if (grad && !floored) {
      backward_solve(L_B_.data(), m, u);
      for (size_t a = 0; a < N; ++a) {
        v[a] -= dot(V_B_.data() + a * m, u, m);
      }
      backward_solve(L_.data(), N, v);
      // d k(y, x) / d x_i = -c (x_i - y_i) / l_i^2, for y the N points of X with weights v, then those of B with u:
      for (size_t a = 0; a < N + m; ++a) {
        double const *ya = a < N ? Xt_.data() + a * n_ : B_.data() + (a - N) * n_;
        double const w = a < N ? 2 * v[a] * cx[a] : 2 * u[a - N] * cb[a - N];
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
    backward_solve(L_B_.data(), m_, w);
    for (size_t r = 0; r < N_; ++r) {
      a[r] = dot(V_B_.data() + r * m_, w, m_);
    }
    backward_solve(L_.data(), N_, a);
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
    return cholesky(L_.data(), N);
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
    cholesky_solve(L_.data(), N, alpha_.data());
    loss = 0;
    for (size_t a = 0; a < N; ++a) {
      loss += (yt_[a] - mean_) * alpha_[a] / 2 + log(L_[a * N + a]);
    }
    loss /= static_cast<double>(N);
    Z_.resize(N * N);
    Kinv_.resize(N * N);
    cholesky_inverse(L_.data(), N, Z_.data(), Kinv_.data());
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
  double outputscale_ = 1, noise_ = 0.005, mean_ = 0;
  std::vector<double> sqdiff_, pair_m_, pair_dm_;
  std::vector<double> L_, Z_, Kinv_, alpha_;
  std::vector<double> V_, Vt_, covariance_;
  // The hyperparameters of the last fit as the optimizer has them, and the conditioning set.
  std::vector<double> theta_;
  size_t m_ = 0;
  std::vector<double> B_, mu_B_, V_B_, L_B_;
  // Whether the variables of the L-BFGS fit are, for now, the hyperparameters rather than their logarithms.
  bool linear_ = false;
};

} BOOST_MATH_NAMESPACE_END // namespace boost::math::optimization::detail
#endif
