/*
 * Copyright Nick Thompson, 2024
 * Use, modification and distribution are subject to the
 * Boost Software License, Version 1.0. (See accompanying file
 * LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef BOOST_MATH_OPTIMIZATION_CMA_ES_HPP
#define BOOST_MATH_OPTIMIZATION_CMA_ES_HPP
#ifndef BOOST_MATH_BUILD_MODULE
#include <atomic>
#include <cmath>
#include <iostream>
#include <limits>
#include <random>
#include <sstream>
#include <stdexcept>
#include <utility>
#include <vector>
#endif
#include <boost/math/optimization/detail/common.hpp>
#include <boost/math/tools/assert.hpp>

// Follows the notation in:
// https://arxiv.org/pdf/1604.00772.pdf
// This is a (hopefully) faithful reproduction of the pseudocode in the arxiv review
// by Nikolaus Hansen.
// Comments referring to equations all refer to this arxiv review.
// A slide deck by the same author is given here:
// http://www.cmap.polytechnique.fr/~nikolaus.hansen/CmaTutorialGecco2023-no-audio.pdf
// which is also a very useful reference.

#ifndef BOOST_MATH_DEBUG_CMA_ES
#define BOOST_MATH_DEBUG_CMA_ES 0
#endif

BOOST_MATH_NAMESPACE_BEGIN namespace optimization {

namespace detail {

// Computes A = V diag(eigenvalues) V^T for a symmetric n x n matrix A stored row-major,
// by Householder tridiagonalization followed by the implicit QL algorithm;
// see Wilkinson & Reinsch, Handbook for Automatic Computation, Vol. II (tred2 and tql2),
// and Golub & Van Loan, Matrix Computations, 4th ed., Section 8.3.
// V is returned row-major, with the eigenvectors as its columns.
// Returns false if some eigenvalue did not converge within max_iterations QL iterations.
template <typename Real>
bool symmetric_eigendecomposition(std::vector<Real> A, size_t n, std::vector<Real> &V, std::vector<Real> &eigenvalues, size_t max_iterations = 30) {
  using std::abs;
  using std::hypot;
  using std::sqrt;
  using std::max;
  V = std::move(A);
  std::vector<Real> &d = eigenvalues;
  d.resize(n);
  std::vector<Real> e(n, Real(0));
  // Householder reduction to tridiagonal form; d is the diagonal and e the subdiagonal:
  for (size_t j = 0; j < n; ++j) {
    d[j] = V[(n - 1)*n + j];
  }
  for (size_t i = n - 1; i > 0; --i) {
    Real scale = 0;
    Real h = 0;
    for (size_t k = 0; k < i; ++k) {
      scale += abs(d[k]);
    }
    if (scale == 0) {
      e[i] = d[i - 1];
      for (size_t j = 0; j < i; ++j) {
        d[j] = V[(i - 1)*n + j];
        V[i*n + j] = 0;
        V[j*n + i] = 0;
      }
    }
    else {
      for (size_t k = 0; k < i; ++k) {
        d[k] /= scale;
        h += d[k]*d[k];
      }
      Real f = d[i - 1];
      Real g = sqrt(h);
      if (f > 0) {
        g = -g;
      }
      e[i] = scale*g;
      h -= f*g;
      d[i - 1] = f - g;
      for (size_t j = 0; j < i; ++j) {
        e[j] = 0;
      }
      for (size_t j = 0; j < i; ++j) {
        f = d[j];
        V[j*n + i] = f;
        g = e[j] + V[j*n + j]*f;
        for (size_t k = j + 1; k < i; ++k) {
          g += V[k*n + j]*d[k];
          e[k] += V[k*n + j]*f;
        }
        e[j] = g;
      }
      f = 0;
      for (size_t j = 0; j < i; ++j) {
        e[j] /= h;
        f += e[j]*d[j];
      }
      const Real hh = f/(h + h);
      for (size_t j = 0; j < i; ++j) {
        e[j] -= hh*d[j];
      }
      for (size_t j = 0; j < i; ++j) {
        f = d[j];
        g = e[j];
        for (size_t k = j; k < i; ++k) {
          V[k*n + j] -= f*e[k] + g*d[k];
        }
        d[j] = V[(i - 1)*n + j];
        V[i*n + j] = 0;
      }
    }
    d[i] = h;
  }
  // Accumulate the transformations:
  for (size_t i = 0; i + 1 < n; ++i) {
    V[(n - 1)*n + i] = V[i*n + i];
    V[i*n + i] = 1;
    const Real h = d[i + 1];
    if (h != 0) {
      for (size_t k = 0; k <= i; ++k) {
        d[k] = V[k*n + i + 1]/h;
      }
      for (size_t j = 0; j <= i; ++j) {
        Real g = 0;
        for (size_t k = 0; k <= i; ++k) {
          g += V[k*n + i + 1]*V[k*n + j];
        }
        for (size_t k = 0; k <= i; ++k) {
          V[k*n + j] -= g*d[k];
        }
      }
    }
    for (size_t k = 0; k <= i; ++k) {
      V[k*n + i + 1] = 0;
    }
  }
  for (size_t j = 0; j < n; ++j) {
    d[j] = V[(n - 1)*n + j];
    V[(n - 1)*n + j] = 0;
  }
  V[(n - 1)*n + n - 1] = 1;
  // The QL rotations act on columns of V, so apply them to the rows of V^T, which are contiguous:
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < i; ++j) {
      std::swap(V[i*n + j], V[j*n + i]);
    }
  }
  // Implicit QL iterations on the tridiagonal matrix:
  for (size_t i = 1; i < n; ++i) {
    e[i - 1] = e[i];
  }
  e[n - 1] = 0;
  Real f = 0;
  Real tst1 = 0;
  const Real eps = std::numeric_limits<Real>::epsilon();
  for (size_t l = 0; l < n; ++l) {
    // Find a small subdiagonal element:
    tst1 = (max)(tst1, abs(d[l]) + abs(e[l]));
    size_t m = l;
    while (m < n - 1 && !(abs(e[m]) <= eps*tst1)) {
      ++m;
    }
    // If m == l, d[l] is already an eigenvalue; otherwise, iterate:
    size_t iterations = 0;
    while (m > l && abs(e[l]) > eps*tst1) {
      if (iterations++ == max_iterations) {
        break;
      }
      // Compute the implicit shift:
      Real g = d[l];
      Real p = (d[l + 1] - g)/(2*e[l]);
      Real r = hypot(p, Real(1));
      if (p < 0) {
        r = -r;
      }
      d[l] = e[l]/(p + r);
      d[l + 1] = e[l]*(p + r);
      const Real dl1 = d[l + 1];
      Real h = g - d[l];
      for (size_t i = l + 2; i < n; ++i) {
        d[i] -= h;
      }
      f += h;
      p = d[m];
      Real c = 1;
      Real c2 = c;
      Real c3 = c;
      const Real el1 = e[l + 1];
      Real s = 0;
      Real s2 = 0;
      for (size_t i = m; i-- > l;) {
        c3 = c2;
        c2 = c;
        s2 = s;
        g = c*e[i];
        h = c*p;
        r = hypot(p, e[i]);
        e[i + 1] = s*r;
        s = e[i]/r;
        c = p/r;
        p = c*d[i] - s*g;
        d[i + 1] = h + s*(c*g + s*d[i]);
        for (size_t k = 0; k < n; ++k) {
          h = V[(i + 1)*n + k];
          V[(i + 1)*n + k] = s*V[i*n + k] + c*h;
          V[i*n + k] = c*V[i*n + k] - s*h;
        }
      }
      p = -s*s2*c3*el1*e[l]/dl1;
      e[l] = s*p;
      d[l] = c*p;
    }
    if (iterations > max_iterations) {
      return false;
    }
    d[l] += f;
    e[l] = 0;
  }
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < i; ++j) {
      std::swap(V[i*n + j], V[j*n + i]);
    }
  }
  return true;
}

// y = Mx, or y = M^{T}x if transpose is true, for the n x n matrix M stored row-major.
// Both accumulate along the rows of M, so the inner loop is over contiguous memory:
template <typename Real>
void matrix_vector_product(std::vector<Real> const &M, std::vector<Real> const &x, std::vector<Real> &y, bool transpose = false) {
  const size_t n = x.size();
  for (size_t i = 0; i < n; ++i) {
    y[i] = 0;
  }
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < n; ++j) {
      if (transpose) {
        y[j] += M[i*n + j]*x[i];
      }
      else {
        y[i] += M[i*n + j]*x[j];
      }
    }
  }
}

} // namespace detail

BOOST_MATH_EXPORT template <typename ArgumentContainer> struct cma_es_parameters {
  using Real = typename ArgumentContainer::value_type;
  using DimensionlessReal = decltype(Real()/Real());
  ArgumentContainer lower_bounds;
  ArgumentContainer upper_bounds;
  size_t max_generations = 1000;
  ArgumentContainer const *initial_guess = nullptr;
  // In the reference, population size = \lambda.
  // If the population size is zero, it is set to equation (48) of the reference
  // and rounded up to the nearest multiple of threads:
  size_t population_size = 0;
  // In the reference, learning_rate = c_m:
  DimensionlessReal learning_rate = 1;
};

BOOST_MATH_EXPORT template <typename ArgumentContainer>
void validate_cma_es_parameters(cma_es_parameters<ArgumentContainer> &params) {
  using Real = typename ArgumentContainer::value_type;
  using DimensionlessReal = decltype(Real()/Real());
  using std::isfinite;
  using std::isnan;
  using std::log;
  using std::ceil;
  using std::floor;

  std::ostringstream oss;
  detail::validate_bounds(params.lower_bounds, params.upper_bounds);
  if (params.initial_guess) {
    detail::validate_initial_guess(*params.initial_guess, params.lower_bounds, params.upper_bounds);
  }
  const size_t n = params.upper_bounds.size();
  // Equation 48 of the arxiv review:
  if (params.population_size == 0) {
    //auto tmp = 4.0 + floor(3*log(n));
    // But round to the nearest multiple of the thread count:
    //auto k = static_cast<size_t>(std::ceil(tmp/params.threads));
    //params.population_size = k*params.threads;
    params.population_size = static_cast<size_t>(4 + floor(3*log(n)));
  }
  if (params.learning_rate <= DimensionlessReal(0) || !isfinite(params.learning_rate)) {
    oss << __FILE__ << ":" << __LINE__ << ":" << __func__;
    oss << ": The learning rate must be > 0, but got " << params.learning_rate << ".";
    throw std::invalid_argument(oss.str());
  }
}

BOOST_MATH_EXPORT template <typename ArgumentContainer, class Func, class URBG>
ArgumentContainer cma_es(
    const Func cost_function,
    cma_es_parameters<ArgumentContainer> &params,
    URBG &gen,
    std::invoke_result_t<Func, ArgumentContainer> target_value = std::numeric_limits<std::invoke_result_t<Func, ArgumentContainer>>::quiet_NaN(),
    std::atomic<bool> *cancellation = nullptr,
    std::atomic<std::invoke_result_t<Func, ArgumentContainer>> *current_minimum_cost = nullptr,
    std::vector<std::pair<ArgumentContainer, std::invoke_result_t<Func, ArgumentContainer>>> *queries = nullptr)
 {
  using Real = typename ArgumentContainer::value_type;
  using DimensionlessReal = decltype(Real()/Real());
  using ResultType = std::invoke_result_t<Func, ArgumentContainer>;
  using std::abs;
  using std::log;
  using std::exp;
  using std::pow;
  using std::min;
  using std::max;
  using std::sqrt;
  using std::isnan;
  using std::isfinite;
  using std::uniform_real_distribution;
  using std::normal_distribution;
  validate_cma_es_parameters(params);
  // n = dimension of problem:
  const size_t n = params.lower_bounds.size();
  std::atomic<bool> target_attained = false;
  std::atomic<ResultType> lowest_cost = std::numeric_limits<ResultType>::infinity();
  ArgumentContainer best_vector;
  // p_{c} := evolution path, equation (24) of the arxiv review:
  std::vector<DimensionlessReal> p_c(n);
  // p_{\sigma} := conjugate evolution path, equation (31) of the arxiv review:
  std::vector<DimensionlessReal> p_sigma(n);
  if constexpr (detail::has_resize_v<ArgumentContainer>) {
    best_vector.resize(n, std::numeric_limits<Real>::quiet_NaN());
  }
  for (size_t i = 0; i < n; ++i) {
    p_c[i] = DimensionlessReal(0);
    p_sigma[i] = DimensionlessReal(0);
  }
  // Table 1, \mu = floor(\lambda/2):
  size_t mu = params.population_size/2;
  std::vector<DimensionlessReal> w_prime(params.population_size, std::numeric_limits<DimensionlessReal>::quiet_NaN());
  for (size_t i = 0; i < params.population_size; ++i) {
    // Equation (49), but 0-indexed:
    w_prime[i] = log(static_cast<DimensionlessReal>(params.population_size + 1)/(2*(i+1)));
  }
  // Table 1, notes at top:
  DimensionlessReal positive_weight_sum = 0;
  DimensionlessReal sq_weight_sum = 0;
  for (size_t i = 0; i < mu; ++i) {
    BOOST_MATH_ASSERT(w_prime[i] > 0);
    positive_weight_sum += w_prime[i];
    sq_weight_sum += w_prime[i]*w_prime[i];
  }
  DimensionlessReal mu_eff = positive_weight_sum*positive_weight_sum/sq_weight_sum;
  BOOST_MATH_ASSERT(1 <= mu_eff);
  BOOST_MATH_ASSERT(mu_eff <= mu);
  DimensionlessReal negative_weight_sum = 0;
  sq_weight_sum = 0;
  for (size_t i = mu; i < params.population_size; ++i) {
    BOOST_MATH_ASSERT(w_prime[i] <= 0);
    negative_weight_sum += w_prime[i];
    sq_weight_sum += w_prime[i]*w_prime[i];
  }
  DimensionlessReal mu_eff_m = negative_weight_sum*negative_weight_sum/sq_weight_sum;
  // Equation (54):
  DimensionlessReal c_m = params.learning_rate;
  // Equation (55):
  DimensionlessReal c_sigma = (mu_eff + 2)/(n + mu_eff + 5);
  BOOST_MATH_ASSERT(c_sigma < 1);
  DimensionlessReal d_sigma = 1 + 2*(max)(DimensionlessReal(0), sqrt(DimensionlessReal((mu_eff - 1)/(n + 1))) - DimensionlessReal(1)) + c_sigma;
  // Equation (56):
  DimensionlessReal c_c = (4 + mu_eff/n)/(n + 4 + 2*mu_eff/n);
  BOOST_MATH_ASSERT(c_c <= 1);
  // Equation (57):
  DimensionlessReal c_1 = DimensionlessReal(2)/(pow(n + 1.3, 2) + mu_eff);
  // Equation (58)
  DimensionlessReal c_mu = (min)(1 - c_1, 2*(DimensionlessReal(0.25)  + mu_eff  + 1/mu_eff - 2)/((n+2)*(n+2) + mu_eff));
  BOOST_MATH_ASSERT(c_1 + c_mu <= DimensionlessReal(1));
  // Equation (50):
  DimensionlessReal alpha_mu_m = 1 + c_1/c_mu;
  // Equation (51):
  DimensionlessReal alpha_mu_eff_m = 1 + 2*mu_eff_m/(mu_eff + 2);
  // Equation (52):
  DimensionlessReal alpha_m_pos_def = (1- c_1 - c_mu)/(n*c_mu);
  // Equation (53):
  std::vector<DimensionlessReal> weights(params.population_size, std::numeric_limits<DimensionlessReal>::quiet_NaN());
  for (size_t i = 0; i < mu; ++i) {
    weights[i] = w_prime[i]/positive_weight_sum;
  }
  DimensionlessReal min_alpha = (min)(alpha_mu_m, (min)(alpha_mu_eff_m, alpha_m_pos_def));
  for (size_t i = mu; i < params.population_size; ++i) {
    weights[i] = min_alpha*w_prime[i]/abs(negative_weight_sum);
  }
  // mu:= number of parents, lambda := number of offspring.
  // C, and B below, are stored row-major:
  std::vector<DimensionlessReal> C(n*n, DimensionlessReal(0));
  for (size_t i = 0; i < n; ++i) {
    C[i*n + i] = 1;
  }
  ArgumentContainer mean_vector;
  // See the footnote in Figure 6 of the arxiv review:
  // We should consider the more robust initialization described there. . . 
  Real sigma = DimensionlessReal(0.3)*(params.upper_bounds[0] - params.lower_bounds[0]);;
  if (params.initial_guess) {
    mean_vector = *params.initial_guess;
  }
  else {
    mean_vector = detail::random_initial_population(params.lower_bounds, params.upper_bounds, 1, gen)[0];
  }
  auto initial_cost = cost_function(mean_vector);
  if (!isnan(initial_cost)) {
    best_vector = mean_vector;
    lowest_cost = initial_cost;
    if (current_minimum_cost) {
      *current_minimum_cost = initial_cost;
    }
  }
#if BOOST_MATH_DEBUG_CMA_ES
  {
    std::cout << __FILE__ << ":" << __LINE__ << ":" << __func__ << "\n";
    std::cout << "\tRunning a (" << params.population_size/2 << "/" << params.population_size/2 << "_W, " << params.population_size << ")-aCMA Evolutionary Strategy on " << params.threads << " threads.\n";
    std::cout << "\tInitial mean vector: {";
    for (size_t i = 0; i < n - 1; ++i) {
      std::cout << mean_vector[i] << ", ";
    }
    std::cout << mean_vector[n - 1] << "}.\n";
    std::cout << "\tCost: " << lowest_cost << ".\n";
    std::cout << "\tInitial step length: " << sigma << ".\n";
    std::cout << "\tVariance effective selection mass: " << mu_eff << ".\n";
    std::cout << "\tLearning rate for rank-one update of covariance matrix: " << c_1 << ".\n";
    std::cout << "\tLearning rate for rank-mu update of covariance matrix: " << c_mu << ".\n";
    std::cout << "\tDecay rate for cumulation path for step-size control: " << c_sigma << ".\n";
    std::cout << "\tLearning rate for the mean: " << c_m << ".\n";
    std::cout << "\tDamping parameter for step-size update: " << d_sigma << ".\n";
  }
#endif
  size_t generation = 0;

  std::vector<std::vector<DimensionlessReal>> ys(params.population_size);
  std::vector<ArgumentContainer> xs(params.population_size);
  std::vector<ResultType> costs(params.population_size, std::numeric_limits<ResultType>::quiet_NaN());
  std::vector<DimensionlessReal> weighted_avg_y(n);
  std::vector<DimensionlessReal> z(n);
  std::vector<DimensionlessReal> B;
  std::vector<DimensionlessReal> D;
  if constexpr (detail::has_resize_v<ArgumentContainer>) {
    for (auto & x : xs) {
      x.resize(n, std::numeric_limits<Real>::quiet_NaN());
    }
  }
  for (auto & y : ys) {
    y.resize(n);
  }
  normal_distribution<DimensionlessReal> dis(DimensionlessReal(0), DimensionlessReal(1));
  do {
    if (cancellation && *cancellation) {
      break;
    }
    // TODO: The reference contends the following in
    // Section B.2 "Strategy internal numerical effort":
    // "In practice, the re-calculation of B and D needs to be done not until about
    // max(1, floor(1/(10n(c_1+c_mu)))) generations."
    // Note that sigma can be dimensionless, in which case C carries the units.
    // This is a weird decision-we're not gonna do that!
    if (!detail::symmetric_eigendecomposition(C, n, B, D)) {
      std::ostringstream oss;
      oss << __FILE__ << ":" << __LINE__ << ":" << __func__;
      oss << ": Could not decompose the covariance matrix as BDB^{T}.";
      throw std::logic_error(oss.str());
    }
    // The eigenvalues are D^2, in the notation of the survey:
    for (auto & d : D) {
      if (d <= 0 || !isfinite(d)) {
        std::ostringstream oss;
        oss << __FILE__ << ":" << __LINE__ << ":" << __func__;
        oss << ": The covariance matrix is not positive definite. This breaks the evolution path computation downstream.\n";
        oss << "C=\n";
        for (size_t i = 0; i < n; ++i) {
          for (size_t j = 0; j < n; ++j) {
            oss << C[i*n + j] << " ";
          }
          oss << "\n";
        }
        oss << "Eigenvalues: ";
        for (auto & e : D) {
          oss << e << " ";
        }
        throw std::domain_error(oss.str());
      }
      d = sqrt(d);
    }

    for (size_t k = 0; k < params.population_size; ++k) {
      auto & y = ys[k];
      auto & x = xs[k];
      BOOST_MATH_ASSERT(static_cast<size_t>(x.size()) == n);
      BOOST_MATH_ASSERT(static_cast<size_t>(y.size()) == n);
      size_t resample_counter = 0;
      do {
        // equation (39) of Figure 6:
        // Also see equation (4):
        for (size_t i = 0; i < n; ++i) {
          z[i] = dis(gen);
        }
        for (size_t i = 0; i < n; ++i) {
          z[i] *= D[i];
        }
        detail::matrix_vector_product(B, z, y);
        for (size_t i = 0; i < n; ++i) {
          BOOST_MATH_ASSERT(!isnan(mean_vector[i]));
          BOOST_MATH_ASSERT(!isnan(y[i]));
          x[i] = mean_vector[i] + sigma*y[i]; // equation (40) of Figure 6.
          // Repair by projecting onto the bounds, and update with the repaired step;
          // see "Boundaries and constraints" in Appendix B of the arxiv review:
          if (x[i] < params.lower_bounds[i] || x[i] > params.upper_bounds[i]) {
            x[i] = (max)(params.lower_bounds[i], (min)(params.upper_bounds[i], x[i]));
            y[i] = (x[i] - mean_vector[i])/sigma;
          }
        }
        costs[k] = cost_function(x);
        if (resample_counter++ == 50) {
          std::ostringstream oss;
          oss << __FILE__ << ":" << __LINE__ << ":" << __func__;
          oss << ": 50 resamples was not sufficient to find an argument to the cost function which did not return NaN.";
          oss << " Giving up.";
          throw std::domain_error(oss.str());
        }
      } while (isnan(costs[k]));

      if (queries) {
        queries->emplace_back(std::make_pair(x, costs[k]));
      }
      if (costs[k] < lowest_cost) {
        lowest_cost = costs[k];
        best_vector = x;
        if (current_minimum_cost && costs[k] < *current_minimum_cost) {
          *current_minimum_cost = costs[k];
        }
        if (lowest_cost < target_value) {
          target_attained = true;
          break;
        }
      }
    }
    if (target_attained) {
      break;
    }
    if (cancellation && *cancellation) {
      break;
    }
    auto indices = detail::best_indices(costs);
    // Equation (41), Figure 6:
    for (size_t j = 0; j < n; ++j) {
      weighted_avg_y[j] = 0;
      for (size_t i = 0; i < mu; ++i) {
        BOOST_MATH_ASSERT(!isnan(weights[i]));
        BOOST_MATH_ASSERT(!isnan(ys[indices[i]][j]));
        weighted_avg_y[j] += weights[i]*ys[indices[i]][j];
      }
    }
    // Equation (42), Figure 6:
    for (size_t j = 0; j < n; ++j) {
      const Real previous_mean = mean_vector[j];
      mean_vector[j] = mean_vector[j] + c_m*sigma*weighted_avg_y[j];
      // For c_m <= 1, the new mean is a convex combination of the old mean and points in the bounds,
      // so it can only leave the bounds by roundoff or if c_m > 1. Then project it,
      // and give the evolution paths the step that the mean actually took:
      if (mean_vector[j] < params.lower_bounds[j] || mean_vector[j] > params.upper_bounds[j]) {
        mean_vector[j] = (max)(params.lower_bounds[j], (min)(params.upper_bounds[j], mean_vector[j]));
        weighted_avg_y[j] = (mean_vector[j] - previous_mean)/(c_m*sigma);
      }
    }
    // Equation (43), Figure 6: Start with C^{-1/2}<y>_{w}
    std::vector<DimensionlessReal> inv_D_B_transpose_y(n);
    detail::matrix_vector_product(B, weighted_avg_y, inv_D_B_transpose_y, true);
    for (size_t j = 0; j < n; ++j) {
      inv_D_B_transpose_y[j] /= D[j];
    }
    std::vector<DimensionlessReal> C_inv_sqrt_y_avg(n);
    detail::matrix_vector_product(B, inv_D_B_transpose_y, C_inv_sqrt_y_avg);
    // Equation (43), Figure 6:
    DimensionlessReal p_sigma_norm = 0;
    for (size_t j = 0; j < n; ++j) {
      p_sigma[j] = (1-c_sigma)*p_sigma[j] + sqrt(c_sigma*(2-c_sigma)*mu_eff)*C_inv_sqrt_y_avg[j];
      p_sigma_norm += p_sigma[j]*p_sigma[j];
    }
    p_sigma_norm = sqrt(p_sigma_norm);
    // A: Algorithm Summary: E[||N(0,1)||]:
    const DimensionlessReal expectation_norm_0I = sqrt(static_cast<DimensionlessReal>(n))*(DimensionlessReal(1) - DimensionlessReal(1)/(4*n) + DimensionlessReal(1)/(21*n*n));
    // Equation (44), Figure 6:
    sigma = sigma*exp(c_sigma*(p_sigma_norm/expectation_norm_0I -1)/d_sigma);
    // A: Algorithm Summary:
    DimensionlessReal h_sigma = 0;
    DimensionlessReal rhs = (DimensionlessReal(1.4) + DimensionlessReal(2)/(n+1))*expectation_norm_0I*sqrt(1 - pow(1-c_sigma, 2*(generation+1)));
    if (p_sigma_norm < rhs) {
      h_sigma = 1;
    }
    // Equation (45), Figure 6:
    for (size_t j = 0; j < n; ++j) {
      p_c[j] = (1-c_c)*p_c[j] + h_sigma*sqrt(c_c*(2-c_c)*mu_eff)*weighted_avg_y[j];
    }
    DimensionlessReal delta_h_sigma = (1-h_sigma)*c_c*(2-c_c);
    DimensionlessReal weight_sum = 0;
    for (auto & w : weights) {
      weight_sum += w;
    }
    // Equation (47), Figure 6:
    DimensionlessReal K = (1 + c_1*delta_h_sigma - c_1 - c_mu*weight_sum);
    // C is symmetric, so only its upper triangle is updated, and then copied to the lower triangle:
    auto add_outer_product = [&](DimensionlessReal alpha, std::vector<DimensionlessReal> const & v) {
      for (size_t r = 0; r < n; ++r) {
        for (size_t c = r; c < n; ++c) {
          C[r*n + c] += alpha*v[r]*v[c];
        }
      }
    };
    for (auto & c : C) {
      c *= K;
    }
    add_outer_product(c_1, p_c);
    // Incorporate positive weights of Equation (46):
    for (size_t i = 0; i < params.population_size/2; ++i) {
      add_outer_product(c_mu*weights[i], ys[indices[i]]);
    }
    std::vector<DimensionlessReal> D_inv_BTy(n);
    for (size_t i = params.population_size/2; i < params.population_size; ++i) {
      auto const & y = ys[indices[i]];
      detail::matrix_vector_product(B, y, D_inv_BTy, true);
      DimensionlessReal squared_norm = 0;
      for (size_t j = 0; j < n; ++j) {
        D_inv_BTy[j] /= D[j];
        squared_norm += D_inv_BTy[j]*D_inv_BTy[j];
      }
      // A repaired step can be exactly zero:
      if (squared_norm == 0) {
        continue;
      }
      DimensionlessReal K2 = c_mu*weights[i]/squared_norm;
      add_outer_product(K2, y);
    }
    for (size_t r = 0; r < n; ++r) {
      for (size_t c = 0; c < r; ++c) {
        C[r*n + c] = C[c*n + r];
      }
    }
  } while (generation++ < params.max_generations);

  return best_vector;
}

} BOOST_MATH_NAMESPACE_END // namespace boost::math::optimization
#endif
