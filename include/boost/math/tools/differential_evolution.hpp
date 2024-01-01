/*
 * Copyright Nick Thompson, 2023
 * Use, modification and distribution are subject to the
 * Boost Software License, Version 1.0. (See accompanying file
 * LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef BOOST_MATH_TOOLS_DIFFERENTIAL_EVOLUTION_HPP
#define BOOST_MATH_TOOLS_DIFFERENTIAL_EVOLUTION_HPP
#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <limits>
#if __has_include(<omp.h>)
#include <omp.h>
#endif
#include <random>
#include <sstream>
#include <stdexcept>
#include <thread>
#include <type_traits>
#include <utility>
#include <vector>

namespace boost::math::tools {

namespace detail {

template <typename T, typename = void> struct has_resize : std::false_type {};

template <typename T>
struct has_resize<T, std::void_t<decltype(std::declval<T>().resize(size_t{}))>> : std::true_type {};

template <typename T> constexpr bool has_resize_v = has_resize<T>::value;

} // namespace detail

// Storn, R., Price, K. (1997). Differential evolution-a simple and efficient heuristic for global optimization over
// continuous spaces.
// Journal of global optimization, 11, 341-359.
// See:
// https://www.cp.eng.chula.ac.th/~prabhas//teaching/ec/ec2012/storn_price_de.pdf

// We provide the parameters in a struct-there are too many of them and they are too unwieldy to pass individually:
template <typename ArgumentContainer> struct differential_evolution_parameters {
  using Real = typename ArgumentContainer::value_type;
  ArgumentContainer lower_bounds;
  ArgumentContainer upper_bounds;
  Real F = static_cast<Real>(0.65);
  double crossover_ratio = 0.5;
  // Population in each generation:
  size_t NP = 200;

  size_t max_generations = 1000;
#if defined(_OPENMP)
  size_t threads = std::thread::hardware_concurrency();
#else
  size_t threads = 1;
#endif
  ArgumentContainer const *initial_guess = nullptr;
};

template <typename ArgumentContainer>
void validate_differential_evolution_parameters(differential_evolution_parameters<ArgumentContainer> const &de_params) {
  using std::isfinite;
  using std::isnan;
  std::ostringstream oss;
  if (de_params.threads == 0) {
    oss << __FILE__ << ":" << __LINE__ << ":" << __func__;
    oss << ": Requested zero threads to perform the calculation, but at least 1 is required.";
    throw std::invalid_argument(oss.str());
  }
  if (de_params.lower_bounds.size() == 0) {
    oss << __FILE__ << ":" << __LINE__ << ":" << __func__;
    oss << ": The dimension of the problem cannot be zero.";
    throw std::domain_error(oss.str());
  }
  if (de_params.upper_bounds.size() != de_params.lower_bounds.size()) {
    oss << __FILE__ << ":" << __LINE__ << ":" << __func__;
    oss << ": There must be the same number of lower bounds as upper bounds, but given ";
    oss << de_params.upper_bounds.size() << " upper bounds, and " << de_params.lower_bounds.size() << " lower bounds.";
    throw std::domain_error(oss.str());
  }
  for (size_t i = 0; i < de_params.lower_bounds.size(); ++i) {
    auto lb = de_params.lower_bounds[i];
    auto ub = de_params.upper_bounds[i];
    if (lb > ub) {
      oss << __FILE__ << ":" << __LINE__ << ":" << __func__;
      oss << ": The upper bound must be greater than or equal to the lower bound, but the upper bound is " << ub
          << " and the lower is " << lb << ".";
      throw std::domain_error(oss.str());
    }
    if (!isfinite(lb)) {
      oss << __FILE__ << ":" << __LINE__ << ":" << __func__;
      oss << ": The lower bound must be finite, but got " << lb << ".";
      oss << " For infinite bounds, emulate with std::numeric_limits<Real>::lower() or use a standard infinite->finite "
             "transform.";
      throw std::domain_error(oss.str());
    }
    if (!isfinite(ub)) {
      oss << __FILE__ << ":" << __LINE__ << ":" << __func__;
      oss << ": The upper bound must be finite, but got " << ub << ".";
      oss << " For infinite bounds, emulate with std::numeric_limits<Real>::max() or use a standard infinite->finite "
             "transform.";
      throw std::domain_error(oss.str());
    }
  }
  if (de_params.NP < 4) {
    oss << __FILE__ << ":" << __LINE__ << ":" << __func__;
    oss << ": The population size must be at least 4, but requested population size of " << de_params.NP << ".";
    throw std::invalid_argument(oss.str());
  }
  if (de_params.threads > de_params.NP) {
    oss << __FILE__ << ":" << __LINE__ << ":" << __func__;
    oss << ": There must be more individuals in the population than threads.";
    throw std::invalid_argument(oss.str());
  }
  // From: "Differential Evolution: A Practical Approach to Global Optimization (Natural Computing Series)"
  // > The scale factor, F in (0,1+), is a positive real number that controls the rate at which the population evolves.
  // > While there is no upper limit on F, effective values are seldom greater than 1.0.
  // ...
  // Also see "Limits on F", Section 2.5.1:
  // > This discontinuity at F = 1 reduces the number of mutants by half and can result in erratic convergence...
  auto F = de_params.F;
  if (isnan(F) || F >= 1 || F <= 0) {
    oss << __FILE__ << ":" << __LINE__ << ":" << __func__;
    oss << ": F in (0, 1) is required, but got F=" << F << ".";
    throw std::domain_error(oss.str());
  }
  if (de_params.max_generations < 1) {
    oss << __FILE__ << ":" << __LINE__ << ":" << __func__;
    oss << ": There must be at least one generation.";
    throw std::invalid_argument(oss.str());
  }
  if (de_params.initial_guess) {
    auto dimension = de_params.lower_bounds.size();
    if (de_params.initial_guess->size() != dimension) {
      oss << __FILE__ << ":" << __LINE__ << ":" << __func__;
      oss << ": The initial guess must have the same dimensions as the problem,";
      oss << ", but the problem size is " << dimension << " and the initial guess has "
          << de_params.initial_guess->size() << " elements.";
      throw std::domain_error(oss.str());
    }
    auto const &guess = *de_params.initial_guess;
    for (size_t i = 0; i < dimension; ++i) {
      auto lb = de_params.lower_bounds[i];
      auto ub = de_params.upper_bounds[i];
      if (guess[i] < lb || guess[i] > ub) {
        oss << __FILE__ << ":" << __LINE__ << ":" << __func__;
        oss << ": At index " << i << " the initial guess " << guess[i] << " is not in the bounds [" << lb << ", " << ub
            << "].";
        throw std::domain_error(oss.str());
      }
    }
  }
#if !defined(_OPENMP)
  if (de_params.threads != 1) {
    oss << __FILE__ << ":" << __LINE__ << ":" << __func__;
    oss << ": If OpenMP is not available, then there algorithm must run on a single thread, but requested "
        << de_params.threads << " threads.";
    throw std::invalid_argument(oss.str());
  }
#endif
}

template <typename ArgumentContainer, class Func, class URBG>
ArgumentContainer differential_evolution(
    const Func cost_function, differential_evolution_parameters<ArgumentContainer> const &de_params, URBG &g,
    std::invoke_result_t<Func, ArgumentContainer> target_value =
        std::numeric_limits<std::invoke_result_t<Func, ArgumentContainer>>::quiet_NaN(),
    std::atomic<bool> *cancellation = nullptr,
    std::vector<std::pair<ArgumentContainer, std::invoke_result_t<Func, ArgumentContainer>>> *queries = nullptr,
    std::atomic<std::invoke_result_t<Func, ArgumentContainer>> *current_minimum_cost = nullptr) {
  using Real = typename ArgumentContainer::value_type;
  validate_differential_evolution_parameters(de_params);
  constexpr bool has_resize = detail::has_resize_v<ArgumentContainer>;

  using ResultType = std::invoke_result_t<Func, ArgumentContainer>;
  using std::clamp;
  using std::isnan;
  using std::round;
  auto const NP = de_params.NP;
  std::vector<ArgumentContainer> population(NP);
  auto const dimension = de_params.lower_bounds.size();
  for (size_t i = 0; i < population.size(); ++i) {
    if constexpr (has_resize) {
      population[i].resize(dimension);
    } else {
      // Argument type must be known at compile-time; like std::array:
      if (population[i].size() != dimension) {
        std::ostringstream oss;
        oss << __FILE__ << ":" << __LINE__ << ":" << __func__;
        oss << ": For containers which do not have resize, the default size must be the same as the dimension, ";
        oss << "but the default container size is " << population[i].size() << " and the dimension of the problem is "
            << dimension << ".";
        oss << " The function argument container type is " << typeid(ArgumentContainer).name() << ".\n";
        throw std::runtime_error(oss.str());
      }
    }
  }
  // Why don't we provide an option to initialize with (say) a Gaussian distribution?
  // > If the optimum's location is fairly well known,
  // > a Gaussian distribution may prove somewhat faster, although it
  // > may also increase the probability that the population will converge prematurely.
  // > In general, uniform distributions are preferred, since they best reflect
  // > the lack of knowledge about the optimum's location.
  //  - Differential Evolution: A Practical Approach to Global Optimization
  // That said, scipy uses Latin Hypercube sampling and says self-avoiding sequences are preferable.
  // So this is something that could be investigated and potentially improved.
  using std::uniform_real_distribution;
  uniform_real_distribution<Real> dis(Real(0), Real(1));
  for (size_t i = 0; i < population.size(); ++i) {
    for (size_t j = 0; j < dimension; ++j) {
      auto const &lb = de_params.lower_bounds[j];
      auto const &ub = de_params.upper_bounds[j];
      population[i][j] = lb + dis(g) * (ub - lb);
    }
  }
  if (de_params.initial_guess) {
    population[0] = *de_params.initial_guess;
  }

  std::atomic<bool> target_attained = false;
  std::vector<ResultType> cost(NP, std::numeric_limits<ResultType>::quiet_NaN());
#if defined(_OPENMP)
#pragma omp parallel for num_threads(de_params.threads)
#endif
  for (size_t i = 0; i < cost.size(); ++i) {
    cost[i] = cost_function(population[i]);
    if (!isnan(target_value) && cost[i] <= target_value) {
      target_attained = true;
    }
    if (current_minimum_cost && cost[i] < *current_minimum_cost) {
      *current_minimum_cost = cost[i];
    }
    if (queries) {
#if defined(_OPENMP) // get rid of -Wunknown-pragmas when OpenMP is not available:
#pragma omp critical
#endif
      queries->push_back(std::make_pair(population[i], cost[i]));
    }
  }

  // This probably won't show up on any performance metrics, but just convert everything to integer ops:
  const auto crossover_int =
      static_cast<decltype(g())>(round(static_cast<double>((URBG::max)() - (URBG::min)()) * de_params.crossover_ratio));
  std::vector<URBG> generators(de_params.threads);
  for (size_t i = 0; i < de_params.threads; ++i) {
    generators[i].seed(g());
  }
  for (size_t generation = 0; generation < de_params.max_generations; ++generation) {
    if (cancellation && *cancellation) {
      break;
    }
    if (target_attained) {
      break;
    }
#if defined(_OPENMP)
#pragma omp parallel for num_threads(de_params.threads)
#endif
    for (size_t i = 0; i < NP; ++i) {
#if defined(_OPENMP)
      size_t thread_idx = omp_get_thread_num();
#else
      size_t thread_idx = 0;
#endif
      auto &gen = generators[thread_idx];
      size_t r1, r2, r3;
      do {
        r1 = gen() % NP;
      } while (r1 == i);
      do {
        r2 = gen() % NP;
      } while (r2 == i || r2 == r1);
      do {
        r3 = gen() % NP;
      } while (r3 == i || r3 == r2 || r3 == r1);
      // Hopefully the compiler optimizes this so that it's not allocated on every iteration:
      ArgumentContainer trial_vector;
      if constexpr (has_resize) {
        trial_vector.resize(dimension);
      }
      for (size_t j = 0; j < dimension; ++j) {
        // See equation (4) of the reference:
        auto guaranteed_changed_idx = gen() % NP;
        if (gen() < crossover_int || j == guaranteed_changed_idx) {
          auto tmp = population[r1][j] + de_params.F * (population[r2][j] - population[r3][j]);
          auto const &lb = de_params.lower_bounds[j];
          auto const &ub = de_params.upper_bounds[j];
          // Some others recommend regenerating the indices rather than clamping;
          // I dunno seems like it could get stuck regenerating . . .
          trial_vector[j] = clamp(tmp, lb, ub);
        } else {
          trial_vector[j] = population[i][j];
        }
      }
      auto trial_cost = cost_function(trial_vector);
      if (queries) {
#if defined(_OPENMP)
#pragma omp critical
#endif
        queries->push_back(std::make_pair(trial_vector, trial_cost));
      }

      if (isnan(trial_cost)) {
        continue;
      }
      if (trial_cost < cost[i] || isnan(cost[i])) {
        std::swap(population[i], trial_vector);
        cost[i] = trial_cost;
        if (!isnan(target_value) && cost[i] <= target_value) {
          target_attained = true;
          // In a single-threaded context, I'd put a break statement here,
          // but OpenMP does not allow break statements in for loops.
          // We'll just have to wait until the end of this generation.
        }
        if (current_minimum_cost && cost[i] < *current_minimum_cost) {
          *current_minimum_cost = cost[i];
        }
      }
    }
  }

  auto it = std::min_element(cost.begin(), cost.end());
  return population[std::distance(cost.begin(), it)];
}

} // namespace boost::math::tools
#endif
