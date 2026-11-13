/*
 * Copyright Nick Thompson, 2026
 * Use, modification and distribution are subject to the
 * Boost Software License, Version 1.0. (See accompanying file
 * LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef BOOST_MATH_OPTIMIZATION_NSGA2_HPP
#define BOOST_MATH_OPTIMIZATION_NSGA2_HPP

#include <boost/math/tools/config.hpp>
#ifndef BOOST_MATH_BUILD_MODULE
#include <atomic>
#endif
#include <boost/math/optimization/detail/common.hpp>
#ifndef BOOST_MATH_BUILD_MODULE
#include <algorithm>
#include <cmath>
#include <iterator>
#include <limits>
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

// Deb, K., Pratap, A., Agarwal, S., Meyarivan, T. (2002). A fast and elitist multiobjective genetic algorithm: NSGA-II.
// IEEE Transactions on Evolutionary Computation, 6(2), 182-197.
// The variation operators (bounded simulated binary crossover and polynomial mutation) follow Deb's reference C code.
// Unlike the paper, children which duplicate a member of the population are not kept.

namespace detail {

// Does a dominate b, all objectives being minimized? A point with a NaN objective is dominated by every point without.
BOOST_MATH_TEST_EXPORT template <class ObjectiveContainer>
bool dominates(ObjectiveContainer const &a, ObjectiveContainer const &b) {
  using std::isnan;
  bool better = false;
  bool worse = false;
  bool b_is_nan = false;
  for (size_t m = 0; m < a.size(); ++m) {
    if (isnan(a[m])) {
      return false;
    }
    b_is_nan = b_is_nan || isnan(b[m]);
    better = better || a[m] < b[m];
    worse = worse || a[m] > b[m];
  }
  return b_is_nan || (better && !worse);
}

// The fast nondominated sort of Section III-A: fronts[0] holds the indices of the nondominated points,
// fronts[1] those dominated only by members of fronts[0], and so on.
BOOST_MATH_TEST_EXPORT template <class ObjectiveContainer>
std::vector<std::vector<size_t>> nondominated_sort(std::vector<ObjectiveContainer> const &values) {
  const size_t n = values.size();
  std::vector<std::vector<size_t>> dominated(n);
  std::vector<size_t> dominators(n, 0);
  std::vector<std::vector<size_t>> fronts(1);
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = i + 1; j < n; ++j) {
      if (dominates(values[i], values[j])) {
        dominated[i].push_back(j);
        ++dominators[j];
      } else if (dominates(values[j], values[i])) {
        dominated[j].push_back(i);
        ++dominators[i];
      }
    }
    if (dominators[i] == 0) {
      fronts[0].push_back(i);
    }
  }
  while (!fronts.back().empty()) {
    std::vector<size_t> next;
    for (size_t p : fronts.back()) {
      for (size_t q : dominated[p]) {
        if (--dominators[q] == 0) {
          next.push_back(q);
        }
      }
    }
    fronts.push_back(std::move(next));
  }
  fronts.pop_back();
  return fronts;
}

// The crowding distance of Section III-B of each member of a front, which is not empty, in the order of the front.
BOOST_MATH_TEST_EXPORT template <class ObjectiveContainer>
auto crowding_distance(std::vector<ObjectiveContainer> const &values, std::vector<size_t> const &front) {
  using Real = typename ObjectiveContainer::value_type;
  using std::isfinite;
  using std::isnan;
  const size_t n = front.size();
  std::vector<Real> distance(n, Real(0));
  std::vector<size_t> order(n);
  for (size_t m = 0; m < values[front[0]].size(); ++m) {
    std::iota(order.begin(), order.end(), size_t(0));
    auto value = [&](size_t i) { return values[front[i]][m]; };
    // NaNs are not ordered; a front holding one holds nothing else, and its distances are of no use.
    if (std::any_of(order.begin(), order.end(), [&](size_t i) { return isnan(value(i)); })) {
      continue;
    }
    // Of equal values, the one that is first in the front is taken to be the lesser:
    std::stable_sort(order.begin(), order.end(), [&](size_t a, size_t b) { return value(a) < value(b); });
    distance[order[0]] = distance[order[n - 1]] = std::numeric_limits<Real>::infinity();
    const Real range = value(order[n - 1]) - value(order[0]);
    if (!(range > 0) || !isfinite(range)) {
      continue;
    }
    for (size_t i = 1; i + 1 < n; ++i) {
      distance[order[i]] += (value(order[i + 1]) - value(order[i - 1])) / range;
    }
  }
  return distance;
}

} // namespace detail

BOOST_MATH_EXPORT template <typename ArgumentContainer> struct nsga2_parameters {
  using Real = typename ArgumentContainer::value_type;
  using DimensionlessReal = decltype(Real()/Real());
  ArgumentContainer lower_bounds;
  ArgumentContainer upper_bounds;
  size_t population_size = 100;
  size_t max_generations = 250;
  // Simulated binary crossover: the probability that a pair of parents is crossed, and the distribution index eta_c.
  DimensionlessReal crossover_probability = static_cast<DimensionlessReal>(0.9);
  DimensionlessReal crossover_distribution_index = static_cast<DimensionlessReal>(15);
  // Polynomial mutation: the probability that a variable is mutated (NaN: 1/dimension), and the distribution index
  // eta_m.
  DimensionlessReal mutation_probability = std::numeric_limits<DimensionlessReal>::quiet_NaN();
  DimensionlessReal mutation_distribution_index = static_cast<DimensionlessReal>(20);
  // Up to population_size members of the first generation; the rest are uniformly distributed in the bounds:
  std::vector<ArgumentContainer> const *initial_population = nullptr;
  unsigned threads = std::thread::hardware_concurrency();
};

BOOST_MATH_EXPORT template <typename ArgumentContainer>
void validate_nsga2_parameters(nsga2_parameters<ArgumentContainer> const &params) {
  using std::isnan;
  std::ostringstream oss;
  detail::validate_bounds(params.lower_bounds, params.upper_bounds);
  if (params.population_size < 4) {
    oss << __FILE__ << ":" << __LINE__ << ":" << __func__;
    oss << ": The population size must be at least 4, but requested population size of " << params.population_size
        << ".";
    throw std::invalid_argument(oss.str());
  }
  if (params.max_generations < 1) {
    oss << __FILE__ << ":" << __LINE__ << ":" << __func__;
    oss << ": There must be at least one generation.";
    throw std::invalid_argument(oss.str());
  }
  if (!(params.crossover_probability >= 0 && params.crossover_probability <= 1)) {
    oss << __FILE__ << ":" << __LINE__ << ":" << __func__;
    oss << ": The crossover probability must be in [0, 1], but got " << params.crossover_probability << ".";
    throw std::domain_error(oss.str());
  }
  if (!isnan(params.mutation_probability) && !(params.mutation_probability >= 0 && params.mutation_probability <= 1)) {
    oss << __FILE__ << ":" << __LINE__ << ":" << __func__;
    oss << ": The mutation probability must be in [0, 1], or NaN for 1/dimension, but got "
        << params.mutation_probability << ".";
    throw std::domain_error(oss.str());
  }
  using std::isfinite;
  if (!(params.crossover_distribution_index >= 0 && isfinite(params.crossover_distribution_index)) ||
      !(params.mutation_distribution_index >= 0 && isfinite(params.mutation_distribution_index))) {
    oss << __FILE__ << ":" << __LINE__ << ":" << __func__;
    oss << ": The distribution indices must be nonnegative and finite, but got " << params.crossover_distribution_index
        << " for the crossover and " << params.mutation_distribution_index << " for the mutation.";
    throw std::domain_error(oss.str());
  }
  if (params.initial_population) {
    if (params.initial_population->size() > params.population_size) {
      oss << __FILE__ << ":" << __LINE__ << ":" << __func__;
      oss << ": The initial population has " << params.initial_population->size()
          << " members, but the population size is " << params.population_size << ".";
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

// NSGA-II in ask and tell form, for objectives the caller evaluates: in batches, on a GPU or a cluster, or between
// checkpoints. ask() returns the points to evaluate, first the initial population and then the children of each
// generation, and tell(values) takes their objectives, in the same order. front() is the nondominated members so far.
// The max_generations and threads of the parameters are not used: the caller decides when to stop.
// The random numbers come from the generator given at construction, of which the optimizer keeps a copy unless URBG
// is a reference type.
BOOST_MATH_EXPORT template <typename ArgumentContainer, typename ObjectiveContainer, class URBG = std::mt19937_64>
class nsga2_optimizer {
public:
  nsga2_optimizer(nsga2_parameters<ArgumentContainer> const &params, URBG gen)
      : params_(params), gen_(std::forward<URBG>(gen)), permutation_(params.population_size),
        drawn_(params.population_size) {
    // These two are not used here, so are not to be objected to:
    params_.threads = 1;
    params_.max_generations = 1;
    validate_nsga2_parameters(params_);
    if (params_.initial_population) {
      pending_ = *params_.initial_population;
      params_.initial_population = nullptr;
    }
    std::iota(permutation_.begin(), permutation_.end(), size_t(0));
  }

  // The same points are returned until their objectives are told.
  std::vector<ArgumentContainer> const &ask() {
    using Real = typename ArgumentContainer::value_type;
    using DimensionlessReal = decltype(Real()/Real());
    using std::clamp;
    using std::expm1;
    using std::isnan;
    using std::log1p;
    using std::pow;
    if (asked_) {
      return pending_;
    }
    const size_t dimension = params_.lower_bounds.size();
    const size_t N = params_.population_size;
    if (population_.empty()) {
      auto initial = detail::random_initial_population(params_.lower_bounds, params_.upper_bounds, N, gen_);
      std::copy(pending_.begin(), pending_.end(), initial.begin());
      pending_ = std::move(initial);
      asked_ = true;
      return pending_;
    }
    const DimensionlessReal mutation_probability =
        isnan(params_.mutation_probability) ? 1 / static_cast<DimensionlessReal>(dimension) : params_.mutation_probability;
    const DimensionlessReal crossover_exponent = params_.crossover_distribution_index + 1;
    const DimensionlessReal mutation_exponent = params_.mutation_distribution_index + 1;
    std::uniform_real_distribution<DimensionlessReal> unif01(DimensionlessReal(0), DimensionlessReal(1));
    auto coin = [&]() { return 2 * unif01(gen_) <= 1; };
    // Binary tournaments between successive members of random permutations, so each member takes part in two per
    // generation:
    auto draw = [&]() {
      if (drawn_ == N) {
        std::shuffle(permutation_.begin(), permutation_.end(), gen_);
        drawn_ = 0;
      }
      return permutation_[drawn_++];
    };
    auto tournament = [&]() -> ArgumentContainer const & {
      const size_t a = draw();
      const size_t b = draw();
      return population_[rank_[a] < rank_[b] || (rank_[a] == rank_[b] && crowding_[a] >= crowding_[b]) ? a : b];
    };
    size_t copies = 0;
    while (pending_.size() < N) {
      ArgumentContainer child1 = tournament();
      ArgumentContainer child2 = tournament();
      if (unif01(gen_) <= params_.crossover_probability) {
        for (size_t k = 0; k < dimension; ++k) {
          if (!coin() || child1[k] == child2[k]) {
            continue;
          }
          const Real y1 = (std::min)(child1[k], child2[k]);
          const Real y2 = (std::max)(child1[k], child2[k]);
          auto const &lb = params_.lower_bounds[k];
          auto const &ub = params_.upper_bounds[k];
          const DimensionlessReal u = unif01(gen_);
          // The spread factor, whose distribution is cut off so that the child stays in the bounds:
          auto spread = [&](DimensionlessReal beta) {
            const DimensionlessReal alpha = 2 - pow(beta, -crossover_exponent);
            return pow(u * alpha <= 1 ? u * alpha : 1 / (2 - u * alpha), 1 / crossover_exponent);
          };
          child1[k] = clamp((y1 + y2 - spread(1 + 2 * (y1 - lb) / (y2 - y1)) * (y2 - y1)) / 2, lb, ub);
          child2[k] = clamp((y1 + y2 + spread(1 + 2 * (ub - y2) / (y2 - y1)) * (y2 - y1)) / 2, lb, ub);
          if (coin()) {
            std::swap(child1[k], child2[k]);
          }
        }
      }
      for (auto *child : {&child1, &child2}) {
        for (size_t k = 0; k < dimension; ++k) {
          auto const &lb = params_.lower_bounds[k];
          auto const &ub = params_.upper_bounds[k];
          if (!(unif01(gen_) <= mutation_probability) || !(ub > lb)) {
            continue;
          }
          auto &y = (*child)[k];
          const DimensionlessReal u = unif01(gen_);
          // (1 - w + w (1 - d)^e)^(1/e) - 1, for the distance d to the bound the step is towards and the exponent e,
          // written so that it is accurate, and not zero, when d is below the machine epsilon:
          auto step = [&](DimensionlessReal w, DimensionlessReal d) {
            return expm1(log1p(w * expm1(mutation_exponent * log1p(-d))) / mutation_exponent);
          };
          const DimensionlessReal delta = 2 * u <= 1 ? step(1 - 2 * u, (y - lb) / (ub - lb))
                                                     : -step(2 * u - 1, (ub - y) / (ub - lb));
          y = clamp(y + delta * (ub - lb), lb, ub);
        }
      }
      // A copy of a parent or of another child would take a place and an evaluation for nothing, so it is dropped,
      // unless copies are all that the generation yields:
      for (auto *child : {&child1, &child2}) {
        auto in = [&](auto const &members) { return std::find(members.begin(), members.end(), *child) != members.end(); };
        if (pending_.size() < N && ((!in(population_) && !in(pending_)) || ++copies > N)) {
          pending_.push_back(std::move(*child));
        }
      }
    }
    asked_ = true;
    return pending_;
  }

  // The objectives of the points last asked for, in their order. Keeps the population_size best of these and the
  // population: whole fronts while they fit, then the least crowded members of the next.
  void tell(std::vector<ObjectiveContainer> values) {
    const size_t objectives = values_.empty() ? (values.empty() ? 0 : values[0].size()) : values_[0].size();
    auto const bad = std::find_if(values.begin(), values.end(), [&](auto const &v) { return v.size() != objectives; });
    if (!asked_ || values.size() != pending_.size() || objectives == 0 || bad != values.end()) {
      std::ostringstream oss;
      oss << __FILE__ << ":" << __LINE__ << ":" << __func__;
      oss << ": tell takes the objectives of the " << (asked_ ? pending_.size() : 0) << " points last asked for, "
          << objectives << " (not zero) for each, but got those of " << values.size() << " points";
      if (bad != values.end()) {
        oss << ", " << bad->size() << " for one of them";
      }
      oss << ".";
      throw std::domain_error(oss.str());
    }
    // The next generation is built aside and swapped in, so that an exception leaves the optimizer as it was.
    auto all = population_;
    all.insert(all.end(), pending_.begin(), pending_.end());
    values.insert(values.begin(), values_.begin(), values_.end());
    const size_t N = params_.population_size;
    auto const fronts = detail::nondominated_sort(values);
    std::vector<ArgumentContainer> next_population;
    std::vector<ObjectiveContainer> next_values;
    std::vector<size_t> next_rank;
    std::vector<typename ObjectiveContainer::value_type> next_crowding;
    for (size_t r = 0; next_population.size() < N; ++r) {
      auto const &front = fronts[r];
      auto const distance = detail::crowding_distance(values, front);
      std::vector<size_t> order(front.size());
      std::iota(order.begin(), order.end(), size_t(0));
      if (next_population.size() + front.size() > N) {
        std::stable_sort(order.begin(), order.end(), [&](size_t a, size_t b) { return distance[a] > distance[b]; });
        order.resize(N - next_population.size());
      }
      for (size_t i : order) {
        next_rank.push_back(r);
        next_crowding.push_back(distance[i]);
        next_population.push_back(std::move(all[front[i]]));
        next_values.push_back(std::move(values[front[i]]));
      }
    }
    population_.swap(next_population);
    values_.swap(next_values);
    rank_.swap(next_rank);
    crowding_.swap(next_crowding);
    pending_.clear();
    asked_ = false;
  }

  // The nondominated members of the population as (argument, objectives) pairs, in lexicographic order of the
  // objectives, those with equal objectives once. Empty until the first tell.
  std::vector<std::pair<ArgumentContainer, ObjectiveContainer>> front() const {
    using std::isnan;
    std::vector<std::pair<ArgumentContainer, ObjectiveContainer>> front;
    for (size_t i = 0; i < population_.size(); ++i) {
      if (rank_[i] == 0) {
        front.emplace_back(population_[i], values_[i]);
      }
    }
    // There is no ordering NaNs, which are only here if every member has one.
    if (!front.empty() &&
        !std::any_of(front[0].second.begin(), front[0].second.end(), [](auto const &v) { return isnan(v); })) {
      std::stable_sort(front.begin(), front.end(), [](auto const &a, auto const &b) {
        return std::lexicographical_compare(a.second.begin(), a.second.end(), b.second.begin(), b.second.end());
      });
      auto same = [](auto const &a, auto const &b) { return a.second == b.second; };
      front.erase(std::unique(front.begin(), front.end(), same), front.end());
    }
    return front;
  }

  // The whole population after the last tell: its members, their objectives, and the front each is in (0: nondominated).
  std::vector<ArgumentContainer> const &population() const { return population_; }
  std::vector<ObjectiveContainer> const &objectives() const { return values_; }
  std::vector<size_t> const &ranks() const { return rank_; }

private:
  nsga2_parameters<ArgumentContainer> params_;
  URBG gen_;
  // The population_size survivors with their objectives, fronts and crowding distances, and the points asked for:
  std::vector<ArgumentContainer> population_;
  std::vector<ObjectiveContainer> values_;
  std::vector<size_t> rank_;
  std::vector<typename ObjectiveContainer::value_type> crowding_;
  std::vector<ArgumentContainer> pending_;
  bool asked_ = false;
  std::vector<size_t> permutation_;
  size_t drawn_;
};

// Minimizes every component of objectives(x) over the bounds: asks an nsga2_optimizer for max_generations + 1
// batches, evaluates each on params.threads threads, and returns its front.
// All random numbers are drawn on the calling thread, so the result does not depend on the number of threads.
BOOST_MATH_EXPORT template <typename ArgumentContainer, class Func, class URBG>
std::vector<std::pair<ArgumentContainer, std::invoke_result_t<Func, ArgumentContainer>>>
nsga2(const Func objectives, nsga2_parameters<ArgumentContainer> const &params, URBG &gen,
      std::atomic<bool> *cancellation = nullptr,
      std::vector<std::pair<ArgumentContainer, std::invoke_result_t<Func, ArgumentContainer>>> *queries = nullptr) {
  using ObjectiveContainer = std::invoke_result_t<Func, ArgumentContainer>;
  validate_nsga2_parameters(params);
  nsga2_optimizer<ArgumentContainer, ObjectiveContainer, URBG &> optimizer(params, gen);
  for (size_t generation = 0;; ++generation) {
    auto const &points = optimizer.ask();
    std::vector<ObjectiveContainer> values(points.size());
    // The calling thread is one of the workers.
    const size_t threads = params.threads;
    auto work = [&](size_t j) {
      for (size_t i = j; i < points.size(); i += threads) {
        values[i] = objectives(points[i]);
      }
    };
    std::vector<std::thread> thread_pool;
    for (size_t j = 1; j < threads; ++j) {
      thread_pool.emplace_back(work, j);
    }
    work(0);
    for (auto &thread : thread_pool) {
      thread.join();
    }
    for (size_t i = 0; queries && i < points.size(); ++i) {
      queries->emplace_back(points[i], values[i]);
    }
    optimizer.tell(std::move(values));
    if (generation == params.max_generations || (cancellation && *cancellation)) {
      return optimizer.front();
    }
  }
}

} BOOST_MATH_NAMESPACE_END // namespace boost::math::optimization
#endif
