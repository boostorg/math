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
#include <stdexcept>
#include <sstream>
#include <thread>
#include <type_traits>
#include <utility>
#include <vector>

namespace boost::math::tools {

namespace detail {

template <typename T, typename = void>
struct has_resize : std::false_type {};

template <typename T>
struct has_resize<T, std::void_t<decltype(std::declval<T>().resize(size_t{}))>> : std::true_type {};

template <typename T>
constexpr bool has_resize_v = has_resize<T>::value;

} //namespace detail

// Storn, R., Price, K. (1997). Differential evolution-a simple and efficient heuristic for global optimization over
// continuous spaces.
// Journal of global optimization, 11, 341-359.
// See:
// https://www.cp.eng.chula.ac.th/~prabhas//teaching/ec/ec2012/storn_price_de.pdf
template <typename BoundsContainer>
class differential_evolution {
public:
  using BoundType = typename BoundsContainer::value_type;
  using Real = typename BoundType::value_type;
  differential_evolution(BoundsContainer bounds, Real F = 0.65, double crossover_ratio = 0.5, size_t NP = 200,
                         size_t max_generations = 1000, size_t threads = std::thread::hardware_concurrency())
      : bounds_{bounds}, F_{F}, CR_{crossover_ratio}, NP_{NP}, max_generations_{max_generations}, threads_{threads} {
    using std::isfinite;
    using std::isnan;
    std::ostringstream oss;
    // hardware_concurrency() is allowed to return 0:
    if (threads_ == 0) {
       threads_ = 1;
    }
    if (bounds_.size() == 0) {
       oss << __FILE__ << ":" << __LINE__ << ":" << __func__;
       oss << ": The dimension of the problem cannot be zero.";
       throw std::domain_error(oss.str());
    }
    if (bounds_[0].size() != 2) {
       oss << __FILE__ << ":" << __LINE__ << ":" << __func__;
       oss << ": Each element of the bounds container must have two elements.";
       throw std::invalid_argument(oss.str());
    }
    for (auto const &bound : bounds_) {
      if (bound[0] > bound[1]) {
        oss << __FILE__ << ":" << __LINE__ << ":" << __func__;
        oss << ": The upper bound must be greater than or equal to the lower bound, but the upper bound is "
            << bound[1] << " and the lower is " << bound[0] << ".";
        throw std::domain_error(oss.str());
      }
      if (!isfinite(bound[0])) {
        oss << __FILE__ << ":" << __LINE__ << ":" << __func__;
        oss << ": The lower bound must be finite, but got " << bound[1] << "."; 
        oss << " For infinite bounds, emulate with std::numeric_limits<Real>::lower() or use a standard infinite->finite transform.";
        throw std::domain_error(oss.str());
      }
      if (!isfinite(bound[1])) {
        oss << __FILE__ << ":" << __LINE__ << ":" << __func__;
        oss << ": The upper bound must be finite, but got " << bound[1] << "."; 
        oss << " For infinite bounds, emulate with std::numeric_limits<Real>::max() or use a standard infinite->finite transform.";
        throw std::domain_error(oss.str());
      }
    }
    if (NP_ < 4) {
      oss << __FILE__ << ":" << __LINE__ << ":" << __func__;
      oss << ": The population size must be at least 4.";
      throw std::invalid_argument(oss.str());
    }
    if (threads_ > NP_) {
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
    if (isnan(F_) || F_ >= Real(1) || F_ <= Real(0)) {
      oss << __FILE__ << ":" << __LINE__ << ":" << __func__;
      oss << ": F in (0, 1] is required, but got F=" << F << ".";
      throw std::domain_error(oss.str());
    }
    if (max_generations_ < 1) {
       oss << __FILE__ << ":" << __LINE__ << ":" << __func__;
       oss << ": There must be at least one generation.";
       throw std::invalid_argument(oss.str());
    }
  }

  template <typename ArgumentContainer, class Func, class URBG>
  ArgumentContainer
  argmin(const Func cost_function, URBG &g,
           ArgumentContainer* initial_guess= nullptr,
           std::invoke_result_t<Func, ArgumentContainer> target_value = std::numeric_limits<std::invoke_result_t<Func, ArgumentContainer>>::quiet_NaN(),
           std::atomic<bool> *cancellation = nullptr,
           std::vector<std::pair<ArgumentContainer, std::invoke_result_t<Func, ArgumentContainer>>>
               *queries = nullptr,
          std::atomic<std::invoke_result_t<Func, ArgumentContainer>>* current_minimum_cost = nullptr
         ) {
    constexpr bool has_resize = detail::has_resize_v<ArgumentContainer>;

    using ResultType = std::invoke_result_t<Func, ArgumentContainer>;
    using std::clamp;
    using std::round;
    using std::isnan;
    std::vector<ArgumentContainer> population(NP_);
    for (size_t i = 0; i < population.size(); ++i) {
      if constexpr (has_resize) {
        // Resize it to same size as bounds_:
        population[i].resize(bounds_.size());
      } else {
        // Argument type must be known at compile-time; like std::array:
        if (population[i].size() != bounds_.size()) {
          std::ostringstream oss;
          oss << __FILE__ << ":" << __LINE__ << ":" << __func__;
          oss << ": For containers which do not have resize, the default size must be the same as the dimension, ";
          oss << "but the default container size is " << population[i].size() << " and the dimension of the problem is " << bounds_.size() << ".";
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
    using std::uniform_real_distribution;
    uniform_real_distribution<Real> dis(Real(0), Real(1));
    for (size_t i = 0; i < population.size(); ++i) {
      for (size_t j = 0; j < bounds_.size(); ++j) {
        auto const &bound = bounds_[j];
        population[i][j] = bound[0] + dis(g) * (bound[1] - bound[0]);
      }
    }
    if (initial_guess) {
        population[0] = *initial_guess;
    }

    std::atomic<bool> target_attained = false;
    std::vector<ResultType> cost(NP_);
    for (size_t i = 0; i < cost.size(); ++i) {
      cost[i] = cost_function(population[i]);
      if (!isnan(target_value) && cost[i] <= target_value) {
         target_attained = true;
      }
      if (current_minimum_cost && cost[i] < *current_minimum_cost) {
          *current_minimum_cost = cost[i];
      }
      if (queries) {
        queries->push_back(std::make_pair(population[i], cost[i]));
      }
    }

    const size_t dimension = bounds_.size();
    const auto crossover_int = static_cast<decltype(g())>(round(static_cast<double>( (URBG::max)() - (URBG::min)()) * CR_));
    std::vector<URBG> generators(threads_);
    for (size_t i = 0; i < threads_; ++i) {
       generators[i].seed(g());
    }
    for (size_t generation = 0; generation < max_generations_; ++generation) {
      if (cancellation && *cancellation) {
          break;
      }
      if (target_attained) {
         break;
      }
#pragma omp parallel for num_threads(threads_)
      for (size_t i = 0; i < NP_; ++i) {
        #if defined(_OPENMP)
        size_t thread_idx = omp_get_thread_num();
        #else 
        size_t thread_idx = 0;
        #endif
        auto& gen = generators[thread_idx];
        size_t r1, r2, r3;
        do {
          r1 = gen() % NP_;
        } while (r1 == i);
        do {
          r2 = gen() % NP_;
        } while (r2 == i || r2 == r1);
        do {
          r3 = gen() % NP_;
        } while (r3 == i || r3 == r2 || r3 == r1);
        ArgumentContainer trial_vector;
        if constexpr (has_resize) {
          trial_vector.resize(bounds_.size());
        }
        for (size_t j = 0; j < dimension; ++j) {
          // See equation (4) of the reference:
          auto guaranteed_changed_idx = gen() % NP_;
          if (gen() < crossover_int || j == guaranteed_changed_idx) {
            auto tmp = population[r1][j] + F_ * (population[r2][j] - population[r3][j]);
            trial_vector[j] = clamp(tmp, bounds_[j][0], bounds_[j][1]);
          } else {
            trial_vector[j] = population[i][j];
          }
        }
        auto trial_cost = cost_function(trial_vector);
        if (queries) {
          #pragma omp critical
          queries->push_back(std::make_pair(trial_vector, trial_cost));
        }
 
        if (trial_cost < cost[i]) {
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

private:
  BoundsContainer bounds_;
  Real F_;
  double CR_;
  size_t NP_;
  size_t max_generations_;
  size_t threads_;
};

} // namespace boost::math::tools
#endif
