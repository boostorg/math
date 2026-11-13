/*
 * Copyright Nick Thompson, 2026
 * Use, modification and distribution are subject to the
 * Boost Software License, Version 1.0. (See accompanying file
 * LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
 */

#ifndef BOOST_MATH_BUILD_MODULE
#include <boost/math/optimization/nsga2.hpp>
#else
import boost.math;
#endif

#include "math_unit_test.hpp"
#include "test_functions_for_optimization.hpp"
#ifndef BOOST_MATH_BUILD_MODULE
#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <limits>
#include <random>
#include <stdexcept>
#include <utility>
#include <vector>
#endif
using boost::math::optimization::nsga2;
using boost::math::optimization::nsga2_optimizer;
using boost::math::optimization::nsga2_parameters;
using boost::math::optimization::validate_nsga2_parameters;
using boost::math::optimization::detail::crowding_distance;
using boost::math::optimization::detail::dominates;
using boost::math::optimization::detail::nondominated_sort;

// The random number generators of the standard library are fixed, but its distributions are not, so the same seed
// gives different runs on different standard libraries. The thresholds below hold, with a wide margin, for the 60
// seeds we tried; the comments give the worst of them.

// The objectives of a point, evaluated twice. Where the compiler contracts floating point expressions, two copies of
// the same code need not agree to the last bit.
template <class Objectives> bool agree(Objectives const &a, Objectives const &b) {
  using Real = typename Objectives::value_type;
  using std::abs;
  bool close = a.size() == b.size();
  for (size_t m = 0; close && m < a.size(); ++m) {
    const Real scale = (std::max)({Real(1), abs(a[m]), abs(b[m])});
    close = a[m] == b[m] || abs(a[m] - b[m]) <= 16 * std::numeric_limits<Real>::epsilon() * scale;
  }
  return close;
}

// Every point is in the bounds, no point dominates another, no two share their objectives, and the objectives are
// those of the argument.
template <class ArgType, class Func>
void check_front(std::vector<std::pair<ArgType, std::invoke_result_t<Func, ArgType>>> const &front, Func f,
                 nsga2_parameters<ArgType> const &params) {
  CHECK_GE(front.size(), size_t(1));
  CHECK_LE(front.size(), params.population_size);
  bool valid = true;
  for (auto const &p : front) {
    for (size_t i = 0; i < params.lower_bounds.size(); ++i) {
      valid = valid && params.lower_bounds[i] <= p.first[i] && p.first[i] <= params.upper_bounds[i];
    }
    valid = valid && agree(f(p.first), p.second);
    for (auto const &q : front) {
      valid = valid && !dominates(p.second, q.second) && (&p == &q || p.second != q.second);
    }
  }
  CHECK_TRUE(valid);
}

// The area dominated by a front of two objectives, in lexicographic order, inside the box below (r1, r2).
template <class Front> double hypervolume(Front const &front, double r1 = 1.1, double r2 = 1.1) {
  double volume = 0;
  double previous = r2;
  for (auto const &p : front) {
    if (p.second[0] < r1 && p.second[1] < previous) {
      volume += (r1 - p.second[0]) * (previous - p.second[1]);
      previous = p.second[1];
    }
  }
  return volume;
}

void test_dominates() {
  using A = std::array<double, 2>;
  double nan = std::numeric_limits<double>::quiet_NaN();
  CHECK_TRUE(dominates(A{1, 2}, A{1, 3}));
  CHECK_TRUE(dominates(A{0, 2}, A{1, 3}));
  CHECK_TRUE(!dominates(A{1, 3}, A{1, 2}));
  CHECK_TRUE(!dominates(A{1, 2}, A{1, 2}));
  CHECK_TRUE(!dominates(A{0, 3}, A{1, 2}));
  CHECK_TRUE(!dominates(A{1, 2}, A{0, 3}));
  // A NaN is dominated by anything else, and dominates nothing:
  CHECK_TRUE(dominates(A{5, 5}, A{nan, 0}));
  CHECK_TRUE(dominates(A{5, 5}, A{0, nan}));
  CHECK_TRUE(!dominates(A{nan, 0}, A{5, 5}));
  CHECK_TRUE(!dominates(A{0, nan}, A{nan, 0}));
  CHECK_TRUE(!dominates(A{nan, nan}, A{nan, nan}));
}

template <class Real> void test_nondominated_sort() {
  using A = std::array<Real, 2>;
  Real nan = std::numeric_limits<Real>::quiet_NaN();
  // Three staircases, a duplicate, and a NaN, which is alone at the back:
  std::vector<A> values = {{3, 3}, {1, 4}, {nan, 0}, {2, 2}, {4, 1}, {5, 5}, {2, 5}, {2, 2}, {4, 4}};
  auto fronts = nondominated_sort(values);
  for (auto &front : fronts) {
    std::sort(front.begin(), front.end());
  }
  std::vector<std::vector<size_t>> expected = {{1, 3, 4, 7}, {0, 6}, {8}, {5}, {2}};
  CHECK_TRUE(fronts == expected);
  CHECK_TRUE(nondominated_sort(std::vector<A>{}).empty());

  // Against the definition: a front is what nothing left dominates. Few distinct values, so that there are ties.
  std::mt19937_64 gen(1234);
  std::uniform_int_distribution<int> dis(0, 12);
  for (size_t M : {1, 2, 3, 4}) {
    std::vector<std::vector<Real>> random(200, std::vector<Real>(M));
    for (auto &v : random) {
      for (auto &x : v) {
        x = static_cast<Real>(dis(gen));
      }
    }
    fronts = nondominated_sort(random);
    std::vector<int> assigned(random.size(), 0);
    bool agree = true;
    for (auto const &front : fronts) {
      agree = agree && !front.empty();
      std::vector<int> in_front(random.size(), 0);
      for (size_t i : front) {
        in_front[i] = 1;
      }
      for (size_t i = 0; i < random.size(); ++i) {
        bool dominated = false;
        for (size_t j = 0; j < random.size(); ++j) {
          dominated = dominated || (!assigned[j] && dominates(random[j], random[i]));
        }
        agree = agree && (assigned[i] || in_front[i] == !dominated);
      }
      for (size_t i : front) {
        agree = agree && !assigned[i];
        assigned[i] = 1;
      }
    }
    CHECK_TRUE(agree);
    CHECK_EQUAL(std::count(assigned.begin(), assigned.end(), 1), std::ptrdiff_t(random.size()));
  }
}

template <class Real> void test_crowding_distance() {
  using A = std::array<Real, 2>;
  Real inf = std::numeric_limits<Real>::infinity();
  // The front is (0, 10), (1, 6), (3, 3), (4, 0), whose ranges are 4 and 10; the other points are not in it.
  std::vector<A> values = {{3, 3}, {9, 9}, {4, 0}, {0, 10}, {7, 7}, {1, 6}};
  auto d = crowding_distance(values, std::vector<size_t>{0, 2, 3, 5});
  CHECK_EQUAL(d.size(), size_t(4));
  CHECK_ULP_CLOSE(Real(3) / 4 + Real(6) / 10, d[0], 2);
  CHECK_EQUAL(d[1], inf);
  CHECK_EQUAL(d[2], inf);
  CHECK_ULP_CLOSE(Real(3) / 4 + Real(7) / 10, d[3], 2);
  // One or two members are all boundary:
  CHECK_EQUAL(crowding_distance(values, std::vector<size_t>{4})[0], inf);
  d = crowding_distance(values, std::vector<size_t>{2, 3});
  CHECK_EQUAL(d[0], inf);
  CHECK_EQUAL(d[1], inf);
  // An objective of zero range, or of infinite range, adds nothing to the interior:
  values = {{1, 0}, {1, 1}, {1, 3}, {1, 4}, {1, inf}};
  d = crowding_distance(values, std::vector<size_t>{0, 1, 2, 3});
  CHECK_ULP_CLOSE(Real(3) / 4, d[1], 2);
  CHECK_ULP_CLOSE(Real(3) / 4, d[2], 2);
  d = crowding_distance(values, std::vector<size_t>{0, 1, 2, 3, 4});
  CHECK_EQUAL(d[1], Real(0));
  CHECK_EQUAL(d[4], inf);
  // NaNs are not sorted:
  Real nan = std::numeric_limits<Real>::quiet_NaN();
  values = {{nan, 1}, {2, nan}, {nan, nan}};
  d = crowding_distance(values, std::vector<size_t>{0, 1, 2});
  CHECK_TRUE(d == std::vector<Real>(3, Real(0)));
}

// Schaffer's problem, with an argument type that cannot be resized and objectives of another type.
// The Pareto set is 0 <= x <= 2.
void test_schaffer() {
  using ArgType = std::array<float, 1>;
  auto f = [](ArgType const &x) {
    double t = x[0];
    return std::array<double, 2>{t * t, (t - 2) * (t - 2)};
  };
  nsga2_parameters<ArgType> params;
  params.lower_bounds = {-10};
  params.upper_bounds = {10};
  params.population_size = 40;
  params.max_generations = 50;
  std::mt19937_64 gen(12345);
  auto front = nsga2(f, params, gen);
  check_front(front, f, params);
  CHECK_GE(front.size(), size_t(30));
  // In lexicographic order; the ends are x = 0 and x = 2. The worst of 60 seeds are 0.004 away.
  auto less = [](auto const &a, auto const &b) { return a.second < b.second; };
  CHECK_TRUE(std::is_sorted(front.begin(), front.end(), less));
  CHECK_ABSOLUTE_ERROR(0.0f, front.front().first[0], 0.02f);
  CHECK_ABSOLUTE_ERROR(2.0f, front.back().first[0], 0.02f);
}

template <class Real> void test_zdt(int problem, double threshold) {
  using ArgType = std::vector<Real>;
  auto f = [problem](ArgType const &x) { return zdt(x, problem); };
  nsga2_parameters<ArgType> params;
  params.lower_bounds.assign(problem == 6 ? 10 : 30, Real(0));
  params.upper_bounds.assign(problem == 6 ? 10 : 30, Real(1));
  std::mt19937_64 gen(problem);
  auto front = nsga2(f, params, gen);
  check_front(front, f, params);
  CHECK_GE(hypervolume(front), threshold);
}

template <class Real> void test_dtlz2() {
  using ArgType = std::vector<Real>;
  nsga2_parameters<ArgType> params;
  params.lower_bounds.assign(12, Real(0));
  params.upper_bounds.assign(12, Real(1));
  std::mt19937_64 gen(77);
  auto front = nsga2(dtlz2<ArgType>, params, gen);
  check_front(front, dtlz2<ArgType>, params);
  // The front is the unit sphere in the positive octant. Over 60 seeds, the member farthest from it is 0.12 away, and
  // the largest value of each objective is at least 0.997.
  std::array<Real, 3> largest = {0, 0, 0};
  for (auto const &p : front) {
    CHECK_EQUAL(p.second.size(), size_t(3));
    using std::hypot;
    CHECK_ABSOLUTE_ERROR(Real(1), hypot(p.second[0], p.second[1], p.second[2]), Real(0.3));
    for (size_t m = 0; m < 3; ++m) {
      largest[m] = (std::max)(largest[m], p.second[m]);
    }
  }
  for (size_t m = 0; m < 3; ++m) {
    CHECK_GE(largest[m], Real(0.95));
  }
}

// The front is known in closed form, and the Pareto set is the diagonal -1 <= sqrt(3) x[i] <= 1. Over 60 seeds, the
// member farthest from the front is 0.036 away, and the argument farthest from the diagonal 0.16: the arguments are
// less accurate than the objectives.
template <class Real> void test_fonseca_fleming() {
  using ArgType = std::vector<Real>;
  using std::exp;
  using std::log1p;
  using std::sqrt;
  nsga2_parameters<ArgType> params;
  params.lower_bounds.assign(3, Real(-4));
  params.upper_bounds.assign(3, Real(4));
  std::mt19937_64 gen(31);
  auto front = nsga2(fonseca_fleming<ArgType>, params, gen);
  check_front(front, fonseca_fleming<ArgType>, params);
  CHECK_ABSOLUTE_ERROR(Real(0), front.front().second[0], Real(0.02));
  CHECK_ABSOLUTE_ERROR(Real(0), front.back().second[1], Real(0.02));
  for (auto const &p : front) {
    Real t = 2 - sqrt(-log1p(-p.second[0]));
    CHECK_ABSOLUTE_ERROR(1 - exp(-t * t), p.second[1], Real(0.1));
    Real mean = (p.first[0] + p.first[1] + p.first[2]) / 3;
    for (Real xi : p.first) {
      CHECK_ABSOLUTE_ERROR(mean, xi, Real(0.5));
    }
  }
}

// The hypervolume below (-14, 1) is between 36.97 and 37.06 over 60 seeds in each of float and double.
template <class Real> void test_kursawe() {
  using ArgType = std::array<Real, 3>;
  nsga2_parameters<ArgType> params;
  params.lower_bounds = {-5, -5, -5};
  params.upper_bounds = {5, 5, 5};
  std::mt19937_64 gen(9);
  auto front = nsga2(kursawe<ArgType>, params, gen);
  check_front(front, kursawe<ArgType>, params);
  CHECK_GE(hypervolume(front, -14, 1), 36.5);
}

// The random numbers are all drawn on the calling thread.
void test_thread_count_independence() {
  using ArgType = std::vector<double>;
  auto f = [](ArgType const &x) { return zdt(x, 3); };
  nsga2_parameters<ArgType> params;
  params.lower_bounds.assign(30, 0.0);
  params.upper_bounds.assign(30, 1.0);
  params.population_size = 37;
  params.max_generations = 20;
  std::vector<std::vector<std::pair<ArgType, std::array<double, 2>>>> fronts;
  std::vector<std::vector<std::pair<ArgType, std::array<double, 2>>>> queries(4);
  for (unsigned threads : {1u, 2u, 5u, 64u}) {
    params.threads = threads;
    std::mt19937_64 gen(2002);
    fronts.push_back(nsga2(f, params, gen, nullptr, &queries[fronts.size()]));
  }
  for (size_t i = 1; i < fronts.size(); ++i) {
    CHECK_TRUE(fronts[i] == fronts[0]);
    CHECK_TRUE(queries[i] == queries[0]);
  }
}

void test_queries_and_cancellation() {
  using ArgType = std::vector<double>;
  nsga2_parameters<ArgType> params;
  params.lower_bounds.assign(5, 0.0);
  params.upper_bounds.assign(5, 1.0);
  // An odd population, which leaves one child of the last pair out:
  params.population_size = 31;
  params.max_generations = 12;
  std::atomic<size_t> calls{0};
  std::atomic<bool> cancellation{false};
  auto f = [&](ArgType const &x) {
    // Asks to stop during the fourth generation of children, which is then finished:
    if (++calls == 4 * params.population_size + 1) {
      cancellation = true;
    }
    return zdt(x, 1);
  };
  std::vector<std::pair<ArgType, std::array<double, 2>>> queries;
  std::mt19937_64 gen(31);
  auto front = nsga2(f, params, gen, nullptr, &queries);
  check_front(front, [](ArgType const &x) { return zdt(x, 1); }, params);
  CHECK_EQUAL(queries.size(), params.population_size * (params.max_generations + 1));
  CHECK_EQUAL(calls.load(), queries.size());
  bool valid = true;
  for (auto const &q : queries) {
    valid = valid && agree(zdt(q.first, 1), q.second);
  }
  CHECK_TRUE(valid);

  calls = 0;
  cancellation = false;
  queries.clear();
  front = nsga2(f, params, gen, &cancellation, &queries);
  CHECK_EQUAL(queries.size(), 5 * params.population_size);
  CHECK_EQUAL(calls.load(), queries.size());

  // Asking and telling by hand:
  using Objectives = std::array<double, 2>;
  params.max_generations = 6;
  std::mt19937_64 gen2(5);
  nsga2_optimizer<ArgType, Objectives> optimizer(params, gen2);
  CHECK_TRUE(optimizer.front().empty());
  CHECK_THROW(optimizer.tell({}), std::domain_error);
  for (size_t generation = 0; generation <= params.max_generations; ++generation) {
    auto const &points = optimizer.ask();
    CHECK_EQUAL(points.size(), params.population_size);
    // The same points until they are told, and only their objectives are taken:
    CHECK_TRUE(&optimizer.ask() == &points);
    std::vector<Objectives> told;
    CHECK_THROW(optimizer.tell(told), std::domain_error);
    CHECK_TRUE(&optimizer.ask() == &points);
    for (auto const &x : points) {
      told.push_back(zdt(x, 1));
    }
    optimizer.tell(told);
    CHECK_THROW(optimizer.tell(told), std::domain_error);
  }
  // The optimizer has its own copy of the generator, where nsga2 advances the caller's:
  CHECK_TRUE(gen2 == std::mt19937_64(5));
  front = optimizer.front();
  check_front(front, [](ArgType const &x) { return zdt(x, 1); }, params);
  // The whole population, front by front:
  CHECK_EQUAL(optimizer.population().size(), params.population_size);
  CHECK_TRUE(std::is_sorted(optimizer.ranks().begin(), optimizer.ranks().end()));
  CHECK_GE(size_t(std::count(optimizer.ranks().begin(), optimizer.ranks().end(), size_t(0))), front.size());
  for (size_t i = 0; i < params.population_size; ++i) {
    CHECK_TRUE(agree(zdt(optimizer.population()[i], 1), optimizer.objectives()[i]));
  }
  // No child is a copy of another of its generation, or of a member of the one before:
  bool distinct = true;
  for (size_t i = params.population_size; i < queries.size(); ++i) {
    for (size_t j = i - i % params.population_size - params.population_size; j < i; ++j) {
      distinct = distinct && queries[i].first != queries[j].first;
    }
  }
  CHECK_TRUE(distinct);
  // Cancelled from the start, the first generation is still evaluated:
  queries.clear();
  front = nsga2(f, params, gen, &cancellation, &queries);
  CHECK_EQUAL(queries.size(), params.population_size);
  CHECK_GE(front.size(), size_t(1));
  // A box of one point, where every child is a copy:
  params.lower_bounds = params.upper_bounds;
  cancellation = false;
  calls = 4 * params.population_size + 1;
  queries.clear();
  front = nsga2(f, params, gen, &cancellation, &queries);
  CHECK_EQUAL(queries.size(), params.population_size * (params.max_generations + 1));
  CHECK_EQUAL(front.size(), size_t(1));
}

template <class Real> void test_warm_start() {
  using ArgType = std::vector<Real>;
  auto f = [](ArgType const &x) { return zdt(x, 1); };
  nsga2_parameters<ArgType> params;
  params.lower_bounds.assign(30, Real(0));
  params.upper_bounds.assign(30, Real(1));
  params.population_size = 20;
  params.max_generations = 30;
  // Members of the Pareto set. Nothing dominates them, and the two at the ends of the front, which are the least in
  // an objective, have infinite crowding distance, so they are never dropped:
  std::vector<ArgType> initial_population(7, ArgType(30, Real(0)));
  for (size_t i = 0; i < initial_population.size(); ++i) {
    initial_population[i][0] = static_cast<Real>(i) / 6;
  }
  params.initial_population = &initial_population;
  std::vector<std::pair<ArgType, std::array<Real, 2>>> queries;
  std::mt19937_64 gen(5);
  auto front = nsga2(f, params, gen, nullptr, &queries);
  check_front(front, f, params);
  for (size_t i = 0; i < initial_population.size(); ++i) {
    CHECK_TRUE(queries[i].first == initial_population[i]);
  }
  CHECK_TRUE(front.front().first == initial_population.front());
  CHECK_TRUE(front.back().first == initial_population.back());
  // A whole population, here of duplicates, which are returned once:
  initial_population.assign(20, ArgType(30, Real(0)));
  front = nsga2(f, params, gen);
  check_front(front, f, params);
}

void test_nan_and_objective_count() {
  using ArgType = std::vector<double>;
  double nan = std::numeric_limits<double>::quiet_NaN();
  nsga2_parameters<ArgType> params;
  params.lower_bounds.assign(3, 0.0);
  params.upper_bounds.assign(3, 1.0);
  params.population_size = 24;
  params.max_generations = 10;
  params.threads = 1;
  std::mt19937_64 gen(8);
  // A point with a NaN objective is not returned while there is any other:
  auto partial = [nan](ArgType const &x) {
    auto y = zdt(x, 2);
    if (x[0] > 0.75) {
      y[x[1] > 0.5] = nan;
    }
    return y;
  };
  auto front = nsga2(partial, params, gen);
  check_front(front, partial, params);
  // Nor do NaNs do any harm when there is nothing else:
  auto all = nsga2([nan](ArgType const &x) { return std::vector<double>{nan, x[0]}; }, params, gen);
  CHECK_EQUAL(all.size(), params.population_size);
  // Infinite objectives are ordered like any other:
  auto infinite = [](ArgType const &x) {
    auto y = zdt(x, 2);
    if (x[0] > 0.25) {
      y[1] = std::numeric_limits<double>::infinity();
    }
    return y;
  };
  check_front(nsga2(infinite, params, gen), infinite, params);

  CHECK_THROW(nsga2([](ArgType const &) { return std::vector<double>{}; }, params, gen), std::domain_error);
  CHECK_THROW(nsga2([](ArgType const &x) { return std::vector<double>(x[0] < 0.5 ? 2 : 3, x[1]); }, params, gen),
              std::domain_error);
}

void test_validation() {
  using ArgType = std::vector<double>;
  double nan = std::numeric_limits<double>::quiet_NaN();
  nsga2_parameters<ArgType> good;
  good.lower_bounds = {0, 0};
  good.upper_bounds = {1, 1};
  CHECK_NO_THROW(validate_nsga2_parameters(good));
  using P = nsga2_parameters<ArgType>;
  auto rejects = [&good](void (*change)(P &, double), double bad = 0) {
    auto p = good;
    change(p, bad);
    try {
      validate_nsga2_parameters(p);
    } catch (std::logic_error const &) {
      return true;
    }
    return false;
  };
  CHECK_TRUE(rejects([](P &p, double) { p.upper_bounds = {1, -1}; }));
  CHECK_TRUE(rejects([](P &p, double) { p.upper_bounds = {1}; }));
  CHECK_TRUE(rejects([](P &p, double) { p.population_size = 3; }));
  CHECK_TRUE(rejects([](P &p, double) { p.max_generations = 0; }));
  CHECK_TRUE(rejects([](P &p, double) { p.threads = 0; }));
  CHECK_TRUE(rejects([](P &p, double) { p.mutation_distribution_index = std::numeric_limits<double>::infinity(); }));
  for (double bad : {-0.1, 1.1, nan}) {
    CHECK_TRUE(rejects([](P &p, double b) { p.crossover_probability = b; }, bad));
    CHECK_TRUE(bad != bad || rejects([](P &p, double b) { p.mutation_probability = b; }, bad));
    CHECK_TRUE(rejects([](P &p, double b) { p.crossover_distribution_index = b > 0 ? -b : b; }, bad));
    CHECK_TRUE(rejects([](P &p, double b) { p.mutation_distribution_index = b > 0 ? -b : b; }, bad));
  }
  auto params = good;
  params.crossover_probability = 0;
  params.mutation_probability = 1;
  params.crossover_distribution_index = 0;
  params.mutation_distribution_index = 0;
  CHECK_NO_THROW(validate_nsga2_parameters(params));
  // The initial population: too many members, a member of the wrong dimension, one out of the bounds, one not finite.
  params = good;
  params.population_size = 4;
  std::vector<ArgType> initial_population(5, ArgType{0.5, 0.5});
  params.initial_population = &initial_population;
  CHECK_THROW(validate_nsga2_parameters(params), std::invalid_argument);
  initial_population.resize(4);
  CHECK_NO_THROW(validate_nsga2_parameters(params));
  initial_population[3] = {0.5};
  CHECK_THROW(validate_nsga2_parameters(params), std::domain_error);
  initial_population[3] = {0.5, 1.5};
  CHECK_THROW(validate_nsga2_parameters(params), std::domain_error);
  initial_population[3] = {nan, 0.5};
  CHECK_THROW(validate_nsga2_parameters(params), std::domain_error);
  // And nsga2 validates:
  std::mt19937_64 gen(0);
  CHECK_THROW(nsga2([](ArgType const &x) { return x; }, params, gen), std::domain_error);
}

int main() {
  test_dominates();
  test_nondominated_sort<float>();
  test_nondominated_sort<double>();
  test_crowding_distance<float>();
  test_crowding_distance<double>();
  test_schaffer();
  // The hypervolumes of the Pareto fronts are 0.8767, 0.5433, 1.3317 and 0.5077; over 60 seeds in each of float and
  // double, the least found are 0.8689, 0.5354, 1.2444 and 0.4937. About one run in sixty misses one of the five
  // pieces of the front of ZDT3, which costs 0.083.
  test_zdt<float>(1, 0.86);
  test_zdt<double>(2, 0.52);
  test_zdt<double>(3, 1.2);
  test_zdt<float>(6, 0.48);
  test_dtlz2<float>();
  test_dtlz2<double>();
  test_fonseca_fleming<float>();
  test_fonseca_fleming<double>();
  test_kursawe<float>();
  test_kursawe<double>();
  test_thread_count_independence();
  test_queries_and_cancellation();
  test_warm_start<float>();
  test_warm_start<double>();
  test_nan_and_objective_count();
  test_validation();
  return boost::math::test::report_errors();
}
