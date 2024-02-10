/*
 * Copyright Nick Thompson, 2024
 * Use, modification and distribution are subject to the
 * Boost Software License, Version 1.0. (See accompanying file
 * LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef BOOST_MATH_OPTIMIZATION_CMA_ES_HPP
#define BOOST_MATH_OPTIMIZATION_CMA_ES_HPP
#include <algorithm>
#include <atomic>
#include <cmath>
#include <iostream>
#include <limits>
#include <random>
#include <sstream>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>
#include <boost/math/optimization/detail/common.hpp>
#include <boost/math/tools/assert.hpp>
#if __has_include(<Eigen/Dense>)
#include <Eigen/Dense>
#else
#error "CMA-ES requires Eigen."
#endif

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

namespace boost::math::optimization {

template <typename ArgumentContainer> struct cma_es_parameters {
  using Real = typename ArgumentContainer::value_type;
  ArgumentContainer lower_bounds;
  ArgumentContainer upper_bounds;
  size_t max_generations = 1000;
  ArgumentContainer const *initial_guess = nullptr;
  // In the reference, population size = \lambda.
  // If the population size is zero, it is set to equation (48) of the reference
  // and rounded up to the nearest multiple of threads:
  size_t population_size = 0;
  // In the reference, learning_rate = c_m:
  Real learning_rate = 1;
};

template <typename ArgumentContainer>
void validate_cma_es_parameters(cma_es_parameters<ArgumentContainer> &params) {
  using Real = typename ArgumentContainer::value_type;
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
  if (params.learning_rate <= Real(0) || !std::isfinite(params.learning_rate)) {
    oss << __FILE__ << ":" << __LINE__ << ":" << __func__;
    oss << ": The learning rate must be > 0, but got " << params.learning_rate << ".";
    throw std::invalid_argument(oss.str());
  }
}

template <typename ArgumentContainer, class Func, class URBG>
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
  Eigen::Vector<Real, Eigen::Dynamic> p_c(n);
  // p_{\sigma} := conjugate evolution path, equation (31) of the arxiv review:
  Eigen::Vector<Real, Eigen::Dynamic> p_sigma(n);
  if constexpr (detail::has_resize_v<ArgumentContainer>) {
    best_vector.resize(n, std::numeric_limits<Real>::quiet_NaN());
  }
  for (size_t i = 0; i < n; ++i) {
    p_c[i] = Real(0);
    p_sigma[i] = Real(0);
  }
  // Table 1, \mu = floor(\lambda/2):
  size_t mu = params.population_size/2;
  std::vector<Real> w_prime(params.population_size, std::numeric_limits<Real>::quiet_NaN());
  for (size_t i = 0; i < params.population_size; ++i) {
    // Equation (49), but 0-indexed:
    w_prime[i] = log(static_cast<Real>(params.population_size + 1)/(2*(i+1)));
  }
  // Table 1, notes at top:
  Real positive_weight_sum = 0;
  Real sq_weight_sum = 0;
  for (size_t i = 0; i < mu; ++i) {
    BOOST_MATH_ASSERT(w_prime[i] > 0);
    positive_weight_sum += w_prime[i];
    sq_weight_sum += w_prime[i]*w_prime[i];
  }
  Real mu_eff = positive_weight_sum*positive_weight_sum/sq_weight_sum;
  BOOST_MATH_ASSERT(1 <= mu_eff);
  BOOST_MATH_ASSERT(mu_eff <= mu);
  Real negative_weight_sum = 0;
  sq_weight_sum = 0;
  for (size_t i = mu; i < params.population_size; ++i) {
    BOOST_MATH_ASSERT(w_prime[i] <= 0);
    negative_weight_sum += w_prime[i];
    sq_weight_sum += w_prime[i]*w_prime[i];
  }
  Real mu_eff_m = negative_weight_sum*negative_weight_sum/sq_weight_sum;
  // Equation (54):
  Real c_m = params.learning_rate;
  // Equation (55):
  Real c_sigma = (mu_eff + 2)/(n + mu_eff + 5);
  BOOST_MATH_ASSERT(c_sigma < 1);
  Real d_sigma = 1 + 2*(max)(Real(0), sqrt((mu_eff - 1)/(n + 1)) - 1) + c_sigma;
  // Equation (56):
  Real c_c = (4 + mu_eff/n)/(n + 4 + 2*mu_eff/n);
  BOOST_MATH_ASSERT(c_c <= 1);
  // Equation (57):
  Real c_1 = Real(2)/(pow(n + 1.3, 2) + mu_eff);
  // Equation (58)
  Real c_mu = (min)(1 - c_1, 2*(Real(0.25)  + mu_eff  + 1/mu_eff - 2)/((n+2)*(n+2) + mu_eff));
  BOOST_MATH_ASSERT(c_1 + c_mu <= Real(1));
  // Equation (50):
  Real alpha_mu_m = 1 + c_1/c_mu;
  // Equation (51):
  Real alpha_mu_eff_m = 1 + 2*mu_eff_m/(mu_eff + 2);
  // Equation (52):
  Real alpha_m_pos_def = (1- c_1 - c_mu)/(n*c_mu);
  // Equation (53):
  std::vector<Real> weights(params.population_size, std::numeric_limits<Real>::quiet_NaN());
  for (size_t i = 0; i < mu; ++i) {
    weights[i] = w_prime[i]/positive_weight_sum;
  }
  Real min_alpha = (min)(alpha_mu_m, (min)(alpha_mu_eff_m, alpha_m_pos_def));
  for (size_t i = mu; i < params.population_size; ++i) {
    weights[i] = min_alpha*w_prime[i]/abs(negative_weight_sum);
  }
  // mu:= number of parents, lambda := number of offspring.
  Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> C = Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic>::Identity(n, n);
  ArgumentContainer mean_vector;
  Real sigma = 1;
  if (params.initial_guess) {
    mean_vector = *params.initial_guess;
  }
  else {
    // See the footnote in Table 1 of the arxiv review:
    sigma = Real(0.3)*(params.upper_bounds[0] - params.lower_bounds[0]);
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

  std::vector<Eigen::Vector<Real, Eigen::Dynamic>> ys(params.population_size);
  std::vector<ArgumentContainer> xs(params.population_size);
  std::vector<ResultType> costs(params.population_size, std::numeric_limits<ResultType>::quiet_NaN());
  Eigen::Vector<Real, Eigen::Dynamic> weighted_avg_y(n);
  Eigen::Vector<Real, Eigen::Dynamic> z(n);
  if constexpr (detail::has_resize_v<ArgumentContainer>) {
    for (auto & x : xs) {
      x.resize(n, std::numeric_limits<Real>::quiet_NaN());
    }
  }
  for (auto & y : ys) {
    y.resize(n);
  }
  normal_distribution<Real> dis(Real(0), Real(1));
  do {
    if (cancellation && *cancellation) {
      break;
    }
    // TODO: The reference contends the following in
    // Section B.2 "Strategy internal numerical effort":
    // "In practice, the re-calculation of B and D needs to be done not until about
    // max(1, floor(1/(10n(c_1+c_mu)))) generations."
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic>> eigensolver(C);
    if (eigensolver.info() != Eigen::Success) {
      std::ostringstream oss;
      oss << __FILE__ << ":" << __LINE__ << ":" << __func__;
      oss << ": Could not decompose the covariance matrix as BDB^{T}.";
      throw std::logic_error(oss.str());
    }
    Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> B = eigensolver.eigenvectors();
    // Eigen returns D^2, in the notation of the survey:
    auto D = eigensolver.eigenvalues();
    // So make it better:
    for (auto & d : D) {
      if (d <= 0 || isnan(d)) {
        std::ostringstream oss;
        oss << __FILE__ << ":" << __LINE__ << ":" << __func__;
        oss << ": The covariance matrix is not positive definite. This breaks the evolution path computation downstream.\n";
        oss << "C=\n" << C << "\n";
        oss << "Eigenvalues: " << D;
        throw std::domain_error(oss.str());
      }
      d = sqrt(d);
    }

    for (size_t k = 0; k < params.population_size; ++k) {
      auto & y = ys[k];
      auto & x = xs[k];
      BOOST_MATH_ASSERT(x.size() == n);
      BOOST_MATH_ASSERT(y.size() == n);
      size_t resample_counter = 0;
      do {
        // equation (39) of Figure 6:
        // Also see equation (4):
        for (size_t i = 0; i < n; ++i) {
          z[i] = dis(gen);
        }
        Eigen::Vector<Real, Eigen::Dynamic> Dz = D.array()*z.array();
        y = B*Dz;
        for (size_t i = 0; i < n; ++i) {
          BOOST_MATH_ASSERT(!isnan(mean_vector[i]));
          BOOST_MATH_ASSERT(!isnan(y[i]));
          x[i] = mean_vector[i] + sigma*y[i]; // equation (40) of Figure 6.
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
      mean_vector[j] = mean_vector[j] + c_m*sigma*weighted_avg_y[j];
    }
    // Equation (43), Figure 6: Start with C^{-1/2}<y>_{w}
    Eigen::Vector<Real, Eigen::Dynamic> inv_D_B_transpose_y = B.transpose()*weighted_avg_y;
    for (size_t j = 0; j < inv_D_B_transpose_y.size(); ++j) {
      inv_D_B_transpose_y[j] /= D[j];
    }
    Eigen::Vector<Real, Eigen::Dynamic> C_inv_sqrt_y_avg = B*inv_D_B_transpose_y;
    // Equation (43), Figure 6:
    Real p_sigma_norm = 0;
    for (size_t j = 0; j < n; ++j) {
      p_sigma[j] = (1-c_sigma)*p_sigma[j] + sqrt(c_sigma*(2-c_sigma)*mu_eff)*C_inv_sqrt_y_avg[j];
      p_sigma_norm += p_sigma[j]*p_sigma[j];
    }
    p_sigma_norm = sqrt(p_sigma_norm);
    // A: Algorithm Summary: E[||N(0,1)||]:
    const Real expectation_norm_0I = sqrt(static_cast<Real>(n))*(Real(1) - Real(1)/(4*n) + Real(1)/(21*n*n));
    // Equation (44), Figure 6:
    sigma = sigma*exp(c_sigma*(p_sigma_norm/expectation_norm_0I -1)/d_sigma);
    // A: Algorithm Summary:
    Real h_sigma = 0;
    Real rhs = (Real(1.4) + Real(2)/(n+1))*expectation_norm_0I*sqrt(1 - pow(1-c_sigma, 2*(generation+1)));
    if (p_sigma_norm < rhs) {
      h_sigma = 1;
    }
    // Equation (45), Figure 6:
    p_c = (1-c_c)*p_c + h_sigma*sqrt(c_c*(2-c_c)*mu_eff)*weighted_avg_y;
    Real delta_h_sigma = (1-h_sigma)*c_c*(2-c_c);
    Real weight_sum = 0;
    for (auto & w : weights) {
      weight_sum += w;
    }
    // Equation (47), Figure 6:
    Real K = (1 + c_1*delta_h_sigma - c_1 - c_mu*weight_sum);
    // Can these operations be sped up using `.selfadjointView<Eigen::Upper>`?
    // Maybe: A.selfadjointView<Eigen::Lower>().rankUpdate(p_c, c_1);?
    C = K*C + c_1*p_c*p_c.transpose();
    // Incorporate positive weights of Equation (46):
    for (size_t i = 0; i < params.population_size/2; ++i) {
      C += c_mu*weights[i]*ys[indices[i]]*ys[indices[i]].transpose();
    }
    for (size_t i = params.population_size/2; i < params.population_size; ++i) {
      Eigen::Vector<Real, Eigen::Dynamic> D_inv_BTy = B.transpose()*ys[indices[i]];
      for (size_t j = 0; j < n; ++j) {
        D_inv_BTy[j] /= D[j];
      }
      Real squared_norm = D_inv_BTy.squaredNorm();
      Real K2 = c_mu*weights[i]/squared_norm;
      C += K2*ys[indices[i]]*ys[indices[i]].transpose();
    }
  } while (generation++ < params.max_generations);

  return best_vector;
}

} // namespace boost::math::optimization
#endif
