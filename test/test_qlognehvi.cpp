/*
 * Copyright Nick Thompson, 2026
 * Use, modification and distribution are subject to the
 * Boost Software License, Version 1.0. (See accompanying file
 * LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
 */

#ifndef BOOST_MATH_BUILD_MODULE
#include <boost/math/optimization/qlognehvi.hpp>
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
using boost::math::optimization::qlognehvi;
using boost::math::optimization::qlognehvi_optimizer;
using boost::math::optimization::qlognehvi_parameters;
using boost::math::optimization::validate_qlognehvi_parameters;
using boost::math::optimization::detail::gaussian_process;
using boost::math::optimization::detail::gaussian_process_parameters;
using boost::math::optimization::detail::hypervolume;
using boost::math::optimization::detail::mark_nondominated;
using boost::math::optimization::detail::qlognehvi_acquisition;

// The random number generators of the standard library are fixed, but its distributions are not, so the same seed
// gives different runs on different standard libraries. The thresholds below hold, with a margin, for the 30 seeds
// we tried; the comments give the worst of them.

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

// Every point is in the bounds, its objectives are finite and those of the argument, no point dominates another, no
// two share their objectives, and they are in lexicographic order.
template <class ArgType, class Func>
void check_front(std::vector<std::pair<ArgType, std::invoke_result_t<Func, ArgType>>> const &front, Func f,
                 qlognehvi_parameters<ArgType> const &params) {
  using std::isfinite;
  CHECK_GE(front.size(), size_t(1));
  bool valid = true;
  for (auto const &p : front) {
    for (size_t i = 0; i < params.lower_bounds.size(); ++i) {
      valid = valid && params.lower_bounds[i] <= p.first[i] && p.first[i] <= params.upper_bounds[i];
    }
    valid = valid && agree(f(p.first), p.second);
    for (auto const &q : front) {
      bool better = false, worse = false;
      for (size_t m = 0; m < p.second.size(); ++m) {
        valid = valid && isfinite(p.second[m]);
        better = better || p.second[m] < q.second[m];
        worse = worse || p.second[m] > q.second[m];
      }
      valid = valid && (&p == &q || better) && (&p == &q || worse);
    }
  }
  CHECK_TRUE(valid);
  CHECK_TRUE(std::is_sorted(front.begin(), front.end(), [](auto const &a, auto const &b) { return a.second < b.second; }));
}

// The hypervolume of the objectives of a front, or of queries, below a reference point.
template <class Pairs> double volume(Pairs const &pairs, std::vector<double> const &reference) {
  std::vector<std::vector<double>> points;
  for (auto const &p : pairs) {
    points.emplace_back(p.second.begin(), p.second.end());
  }
  return hypervolume(points, reference);
}

void test_mark_nondominated() {
  // Against the definition, with few distinct values, so that there are ties, and a reference which cuts some off.
  std::mt19937_64 gen(1234);
  std::uniform_int_distribution<int> dis(0, 9);
  for (size_t m : {2, 3, 4}) {
    size_t const count = 60;
    std::vector<double> f(count * m), reference(m, 8.0);
    for (auto &x : f) {
      x = dis(gen);
    }
    std::vector<size_t> order;
    std::vector<char> kept(count, 0);
    mark_nondominated(f.data(), count, m, reference.data(), order, kept);
    bool right = true;
    for (size_t i = 0; i < count; ++i) {
      bool expected = true;
      for (size_t k = 0; k < m; ++k) {
        expected = expected && f[i * m + k] < reference[k];
      }
      for (size_t j = 0; j < count && expected; ++j) {
        bool better = false, worse = false;
        for (size_t k = 0; k < m; ++k) {
          better = better || f[j * m + k] < f[i * m + k];
          worse = worse || f[j * m + k] > f[i * m + k];
        }
        expected = !(better && !worse);
      }
      right = right && expected == (kept[i] != 0);
    }
    CHECK_TRUE(right);
  }
}

// The gradient of the acquisition function against central differences, with a conditioning set of training points
// and of one other, as when a batch is half chosen.
void test_acquisition_gradient() {
  using std::abs;
  using std::cos;
  using std::sqrt;
  for (size_t m : {2, 3, 4}) {
    size_t const d = 5, n = 20, samples = 16, nb = 9;
    std::mt19937_64 gen(7 + m);
    std::uniform_real_distribution<double> unif(0, 1);
    std::normal_distribution<double> normal(0, 1);
    std::vector<double> X(n * d);
    for (auto &x : X) {
      x = unif(gen);
    }
    std::vector<gaussian_process> models(m, gaussian_process(gaussian_process_parameters::botorch(d)));
    for (size_t k = 0; k < m; ++k) {
      std::vector<double> y(n);
      double mean = 0, variance = 0;
      for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < d; ++j) {
          y[i] += cos(3 * X[i * d + j] + static_cast<double>(k + j)) * (k + 1 == j ? 2 : 0.5);
        }
        mean += y[i] / static_cast<double>(n);
      }
      for (double t : y) {
        variance += (t - mean) * (t - mean) / static_cast<double>(n - 1);
      }
      for (auto &t : y) {
        t = (t - mean) / sqrt(variance);
      }
      CHECK_TRUE(models[k].fit(X.data(), y.data(), n, d));
    }
    std::vector<double> B(X.begin(), X.begin() + 8 * d), zeta(m * samples * nb), z(m * samples);
    for (size_t j = 0; j < d; ++j) {
      B.push_back(unif(gen));
    }
    for (auto &t : zeta) {
      t = normal(gen);
    }
    for (auto &t : z) {
      t = normal(gen);
    }
    qlognehvi_acquisition acquisition(models, std::vector<double>(m, 1.5), samples);
    CHECK_TRUE(acquisition.condition(B, nb, zeta, z, 2));
    qlognehvi_acquisition::workspace ws;
    // The worst over these points and 30 seeds is 6e-7.
    double worst = 0;
    for (int trial = 0; trial < 12; ++trial) {
      std::vector<double> x(d), g(d);
      for (auto &t : x) {
        t = unif(gen);
      }
      if (trial == 11) {
        // Next to a point of the conditioning set, where the conditional variance is nearly zero:
        x.assign(B.begin(), B.begin() + d);
        x[0] += 0.01;
      }
      double const value = acquisition(x.data(), g.data(), ws);
      CHECK_TRUE(std::isfinite(value));
      // The value does not depend on whether the gradient is asked for:
      CHECK_ULP_CLOSE(value, acquisition(x.data(), nullptr, ws), 4);
      for (size_t i = 0; i < d; ++i) {
        double const h = 1e-6;
        auto xp = x, xm = x;
        xp[i] += h;
        xm[i] -= h;
        double const fd = (acquisition(xp.data(), nullptr, ws) - acquisition(xm.data(), nullptr, ws)) / (2 * h);
        worst = (std::max)(worst, abs(fd - g[i]) / (std::max)(1.0, abs(fd)));
      }
    }
    CHECK_LE(worst, 1e-4);
  }
}

// Schaffer's problem, with an argument type that cannot be resized and objectives of another type. The Pareto set
// is 0 <= x <= 2. The hypervolume below (5, 5) is at least 21.4 after 20 evaluations, over 30 seeds.
void test_schaffer() {
  using ArgType = std::array<float, 1>;
  auto f = [](ArgType const &x) {
    double t = x[0];
    return std::array<double, 2>{t * t, (t - 2) * (t - 2)};
  };
  qlognehvi_parameters<ArgType> params;
  params.lower_bounds = {-10};
  params.upper_bounds = {10};
  params.max_evaluations = 20;
  params.batch_size = 2;
  std::mt19937_64 gen(12345);
  auto front = qlognehvi(f, params, gen);
  check_front(front, f, params);
  CHECK_GE(volume(front, {5, 5}), 20.0);
}

// Forty evaluations, in batches of four after the initial design. The hypervolume of the Pareto front of ZDT1 below
// (1.1, 1.1) is 0.877; over 30 seeds this finds between 0.66 and 0.81, where the best of 30 uniform random searches
// has 0.11.
void test_zdt1() {
  using ArgType = std::vector<double>;
  auto f = [](ArgType const &x) { return zdt(x, 1); };
  qlognehvi_parameters<ArgType> params;
  params.lower_bounds.assign(5, 0.0);
  params.upper_bounds.assign(5, 1.0);
  params.max_evaluations = 40;
  params.batch_size = 4;
  std::mt19937_64 gen(1);
  auto front = qlognehvi(f, params, gen);
  check_front(front, f, params);
  CHECK_GE(volume(front, {1.1, 1.1}), 0.5);
}

// Three objectives of four variables in fifty evaluations; the front is the part of the unit sphere in the positive
// octant. The hypervolume below (1.1, 1.1, 1.1) is between 0.54 and 0.68 over 30 seeds; the median of random search
// is 0.44, and its best 0.49.
void test_dtlz2() {
  using ArgType = std::vector<double>;
  qlognehvi_parameters<ArgType> params;
  params.lower_bounds.assign(4, 0.0);
  params.upper_bounds.assign(4, 1.0);
  params.max_evaluations = 50;
  params.batch_size = 4;
  std::mt19937_64 gen(77);
  auto front = qlognehvi(dtlz2<ArgType>, params, gen);
  check_front(front, dtlz2<ArgType>, params);
  CHECK_EQUAL(front[0].second.size(), size_t(3));
  CHECK_GE(volume(front, {1.1, 1.1, 1.1}), 0.46);
}

// The random numbers are all drawn on the calling thread, and the work of each thread does not depend on how many
// there are.
void test_thread_count_independence() {
  using ArgType = std::vector<double>;
  auto f = [](ArgType const &x) { return zdt(x, 3); };
  qlognehvi_parameters<ArgType> params;
  params.lower_bounds.assign(4, 0.0);
  params.upper_bounds.assign(4, 1.0);
  params.max_evaluations = 25;
  params.batch_size = 3;
  std::vector<std::vector<std::pair<ArgType, std::array<double, 2>>>> fronts;
  std::vector<std::vector<std::pair<ArgType, std::array<double, 2>>>> queries(4);
  for (unsigned threads : {1u, 2u, 5u, 64u}) {
    params.threads = threads;
    std::mt19937_64 gen(2002);
    fronts.push_back(qlognehvi(f, params, gen, nullptr, &queries[fronts.size()]));
  }
  for (size_t i = 1; i < fronts.size(); ++i) {
    CHECK_TRUE(fronts[i] == fronts[0]);
    CHECK_TRUE(queries[i] == queries[0]);
  }
  // 10 points of the initial design, then batches of three: the last has two.
  CHECK_EQUAL(queries[0].size(), size_t(25));
}

void test_queries_and_cancellation() {
  using ArgType = std::vector<double>;
  using Objectives = std::array<double, 2>;
  qlognehvi_parameters<ArgType> params;
  params.lower_bounds.assign(2, 0.0);
  params.upper_bounds.assign(2, 1.0);
  params.initial_samples = 7;
  params.batch_size = 3;
  params.max_evaluations = 19;
  std::atomic<size_t> calls{0};
  std::atomic<bool> cancellation{false};
  auto f = [&](ArgType const &x) {
    // Asks to stop during the second batch, which is then finished:
    if (++calls == 11) {
      cancellation = true;
    }
    return zdt(x, 1);
  };
  std::vector<std::pair<ArgType, Objectives>> queries;
  std::mt19937_64 gen(31);
  auto front = qlognehvi(f, params, gen, nullptr, &queries);
  check_front(front, [](ArgType const &x) { return zdt(x, 1); }, params);
  CHECK_EQUAL(queries.size(), params.max_evaluations);
  CHECK_EQUAL(calls.load(), queries.size());
  bool valid = true;
  for (auto const &q : queries) {
    valid = valid && agree(zdt(q.first, 1), q.second);
  }
  CHECK_TRUE(valid);
  // The front is that of the queries:
  CHECK_ABSOLUTE_ERROR(volume(queries, {1.1, 10.0}), volume(front, {1.1, 10.0}), 1e-12);

  calls = 0;
  cancellation = false;
  queries.clear();
  front = qlognehvi(f, params, gen, &cancellation, &queries);
  CHECK_EQUAL(queries.size(), size_t(13));
  CHECK_EQUAL(calls.load(), queries.size());
  // Cancelled from the start, the initial design is still evaluated:
  queries.clear();
  front = qlognehvi(f, params, gen, &cancellation, &queries);
  CHECK_EQUAL(queries.size(), size_t(7));
  CHECK_GE(front.size(), size_t(1));
  // Fewer evaluations than the initial design has points:
  cancellation = false;
  calls = 100;
  params.max_evaluations = 4;
  queries.clear();
  front = qlognehvi(f, params, gen, &cancellation, &queries);
  CHECK_EQUAL(queries.size(), size_t(4));

  // Asking and telling by hand:
  std::mt19937_64 gen2(5);
  qlognehvi_optimizer<ArgType, Objectives> optimizer(params, gen2);
  CHECK_TRUE(optimizer.front().empty());
  CHECK_TRUE(optimizer.reference_point().empty());
  CHECK_THROW(optimizer.tell({}), std::domain_error);
  size_t expected = 0;
  std::vector<double> reference;
  for (size_t batch = 0; batch < 4; ++batch) {
    // A batch of another size, once:
    auto const &points = batch == 2 ? optimizer.ask(5) : optimizer.ask();
    CHECK_EQUAL(points.size(), batch == 0 ? params.initial_samples : batch == 2 ? size_t(5) : params.batch_size);
    // The same points until they are told, and only their objectives are taken:
    CHECK_TRUE(&optimizer.ask() == &points);
    CHECK_TRUE(&optimizer.ask(2) == &points);
    std::vector<Objectives> told;
    CHECK_THROW(optimizer.tell(told), std::domain_error);
    CHECK_TRUE(&optimizer.ask() == &points);
    bool inside = true;
    for (auto const &x : points) {
      inside = inside && x.size() == 2 && 0 <= x[0] && x[0] <= 1 && 0 <= x[1] && x[1] <= 1;
      told.push_back(zdt(x, 1));
    }
    CHECK_TRUE(inside);
    expected += points.size();
    optimizer.tell(told);
    CHECK_THROW(optimizer.tell(told), std::domain_error);
    CHECK_EQUAL(optimizer.points().size(), expected);
    CHECK_EQUAL(optimizer.objectives().size(), expected);
    if (batch == 0) {
      // The reference point: the greatest value of each objective in the initial design, and a tenth of its range.
      double least[2] = {told[0][0], told[0][1]}, greatest[2] = {told[0][0], told[0][1]};
      for (auto const &y : told) {
        for (size_t k = 0; k < 2; ++k) {
          least[k] = (std::min)(least[k], y[k]);
          greatest[k] = (std::max)(greatest[k], y[k]);
        }
      }
      reference = optimizer.reference_point();
      CHECK_EQUAL(reference.size(), size_t(2));
      for (size_t k = 0; k < 2; ++k) {
        CHECK_ULP_CLOSE(greatest[k] + (greatest[k] - least[k]) / 10, optimizer.reference_point()[k], 4);
      }
    }
  }
  for (size_t i = 0; i < expected; ++i) {
    CHECK_TRUE(agree(zdt(optimizer.points()[i], 1), optimizer.objectives()[i]));
  }
  // It is not moved by what comes later:
  CHECK_TRUE(optimizer.reference_point() == reference);
  // The optimizer has its own copy of the generator, where qlognehvi advances the caller's:
  CHECK_TRUE(gen2 == std::mt19937_64(5));
  check_front(optimizer.front(), [](ArgType const &x) { return zdt(x, 1); }, params);
}

// The initial design starts with the points given, a Latin hypercube makes up the number, and a variable whose
// bounds are equal stays put.
template <class Real> void test_initial_population() {
  using ArgType = std::vector<Real>;
  auto f = [](ArgType const &x) { return zdt(x, 1); };
  qlognehvi_parameters<ArgType> params;
  params.lower_bounds = {0, Real(0.25), 0};
  params.upper_bounds = {1, Real(0.25), 1};
  params.max_evaluations = 14;
  params.batch_size = 3;
  std::vector<ArgType> initial_population = {{0, Real(0.25), 0}, {Real(0.5), Real(0.25), 0}, {1, Real(0.25), 0}};
  params.initial_population = &initial_population;
  std::vector<std::pair<ArgType, std::array<Real, 2>>> queries;
  std::mt19937_64 gen(5);
  auto front = qlognehvi(f, params, gen, nullptr, &queries);
  check_front(front, f, params);
  CHECK_EQUAL(queries.size(), size_t(14));
  for (size_t i = 0; i < initial_population.size(); ++i) {
    CHECK_TRUE(queries[i].first == initial_population[i]);
  }
  // The Latin hypercube: of the five other points of the eight, one in each fifth of the first coordinate.
  std::vector<int> slices(5, 0);
  for (size_t i = 3; i < 8; ++i) {
    ++slices[(std::min)(size_t(4), static_cast<size_t>(5 * queries[i].first[0]))];
  }
  CHECK_TRUE(slices == std::vector<int>(5, 1));
  for (auto const &q : queries) {
    CHECK_EQUAL(q.first[1], Real(0.25));
  }
  // The ends of the front were given, and nothing dominates them:
  CHECK_TRUE(front.front().first == initial_population.front());
  CHECK_TRUE(front.back().first == initial_population.back());
  // More points than initial_samples: all are evaluated first, and no others.
  params.initial_samples = 2;
  params.max_evaluations = 5;
  queries.clear();
  qlognehvi(f, params, gen, nullptr, &queries);
  CHECK_EQUAL(queries.size(), size_t(5));
  for (size_t i = 0; i < initial_population.size(); ++i) {
    CHECK_TRUE(queries[i].first == initial_population[i]);
  }
}

void test_nan_and_objective_count() {
  using ArgType = std::vector<double>;
  double nan = std::numeric_limits<double>::quiet_NaN();
  double inf = std::numeric_limits<double>::infinity();
  qlognehvi_parameters<ArgType> params;
  params.lower_bounds.assign(2, 0.0);
  params.upper_bounds.assign(2, 1.0);
  params.max_evaluations = 21;
  params.batch_size = 3;
  params.threads = 1;
  std::mt19937_64 gen(8);
  // A point with an objective that is not finite is not returned, and the others still make a front:
  auto partial = [nan, inf](ArgType const &x) {
    auto y = zdt(x, 2);
    if (x[0] > 0.75) {
      y[x[1] > 0.5] = nan;
    }
    if (x[0] < 0.25) {
      y[1] = x[1] > 0.5 ? inf : -inf;
    }
    return y;
  };
  std::vector<std::pair<ArgType, std::array<double, 2>>> queries;
  auto front = qlognehvi(partial, params, gen, nullptr, &queries);
  check_front(front, partial, params);
  CHECK_EQUAL(queries.size(), params.max_evaluations);
  // Nor do they do any harm when there is nothing else: the points are then random, and the front is empty.
  std::vector<std::pair<ArgType, std::vector<double>>> all_queries;
  auto all = qlognehvi([nan](ArgType const &x) { return std::vector<double>{nan, x[0]}; }, params, gen, nullptr,
                       &all_queries);
  CHECK_TRUE(all.empty());
  CHECK_EQUAL(all_queries.size(), params.max_evaluations);
  // An objective that is constant gives a model nothing to fit, and the other objective is still minimized:
  auto constant = qlognehvi([](ArgType const &x) { return std::vector<double>{1.0, x[0] + x[1]}; }, params, gen);
  CHECK_EQUAL(constant.size(), size_t(1));
  CHECK_LE(constant[0].second[1], 0.3);

  CHECK_THROW(qlognehvi([](ArgType const &) { return std::vector<double>{}; }, params, gen), std::domain_error);
  CHECK_THROW(qlognehvi([](ArgType const &x) { return std::vector<double>{x[0]}; }, params, gen), std::domain_error);
  CHECK_THROW(qlognehvi([](ArgType const &x) { return std::vector<double>(x[0] < 0.5 ? 2 : 3, x[1]); }, params, gen),
              std::domain_error);
  // A reference point for another number of objectives:
  params.reference_point = {1, 1, 1};
  CHECK_THROW(qlognehvi([](ArgType const &x) { return zdt(x, 1); }, params, gen), std::domain_error);
  // One which nothing is below: there is nothing to improve on, and no harm done.
  params.reference_point = {-1, -1};
  front = qlognehvi([](ArgType const &x) { return zdt(x, 1); }, params, gen);
  CHECK_GE(front.size(), size_t(1));
  // An exception of the objectives, on whichever thread, reaches the caller:
  params.threads = 2;
  auto throwing = [](ArgType const &x) { return x[0] > 0.5 ? throw std::runtime_error("objective") : zdt(x, 1); };
  CHECK_THROW(qlognehvi(throwing, params, gen), std::runtime_error);
}

// With a reference point given, the points proposed are below it: here the part of the front of ZDT1 with f1 < 0.5.
// The least hypervolume below (0.5, 1.1) over 30 seeds is 0.27.
void test_reference_point() {
  using ArgType = std::vector<double>;
  auto f = [](ArgType const &x) { return zdt(x, 1); };
  qlognehvi_parameters<ArgType> params;
  params.lower_bounds.assign(3, 0.0);
  params.upper_bounds.assign(3, 1.0);
  params.max_evaluations = 36;
  params.batch_size = 4;
  params.reference_point = {0.5, 1.1};
  std::mt19937_64 gen(11);
  std::vector<std::pair<ArgType, std::array<double, 2>>> queries;
  auto front = qlognehvi(f, params, gen, nullptr, &queries);
  check_front(front, f, params);
  CHECK_GE(volume(front, params.reference_point), 0.2);
  size_t below = 0;
  for (size_t i = 8; i < queries.size(); ++i) {
    below += queries[i].second[0] < 0.5;
  }
  // All 28 of them, in each of 30 seeds.
  CHECK_GE(below, size_t(20));
}

void test_validation() {
  using ArgType = std::vector<double>;
  double nan = std::numeric_limits<double>::quiet_NaN();
  qlognehvi_parameters<ArgType> good;
  good.lower_bounds = {0, 0};
  good.upper_bounds = {1, 1};
  CHECK_NO_THROW(validate_qlognehvi_parameters(good));
  using P = qlognehvi_parameters<ArgType>;
  auto rejects = [&good](void (*change)(P &)) {
    auto p = good;
    change(p);
    try {
      validate_qlognehvi_parameters(p);
    } catch (std::logic_error const &) {
      return true;
    }
    return false;
  };
  CHECK_TRUE(rejects([](P &p) { p.upper_bounds = {1, -1}; }));
  CHECK_TRUE(rejects([](P &p) { p.upper_bounds = {1}; }));
  CHECK_TRUE(rejects([](P &p) { p.lower_bounds = {}; p.upper_bounds = {}; }));
  CHECK_TRUE(rejects([](P &p) { p.batch_size = 0; }));
  CHECK_TRUE(rejects([](P &p) { p.max_evaluations = 0; }));
  CHECK_TRUE(rejects([](P &p) { p.mc_samples = 0; }));
  CHECK_TRUE(rejects([](P &p) { p.restarts = 0; }));
  CHECK_TRUE(rejects([](P &p) { p.raw_samples = 3; }));
  CHECK_TRUE(rejects([](P &p) { p.threads = 0; }));
  CHECK_TRUE(rejects([](P &p) { p.reference_point = {1}; }));
  CHECK_TRUE(rejects([](P &p) { p.reference_point = {1, std::numeric_limits<double>::quiet_NaN()}; }));
  CHECK_TRUE(rejects([](P &p) { p.reference_point = {1, std::numeric_limits<double>::infinity()}; }));
  auto params = good;
  params.initial_samples = 1;
  params.mc_samples = 1;
  params.restarts = 1;
  params.raw_samples = 1;
  params.reference_point = {1, 2, 3};
  CHECK_NO_THROW(validate_qlognehvi_parameters(params));
  // The initial population: more members than evaluations, a member of the wrong dimension, one out of the bounds,
  // one not finite.
  params = good;
  params.max_evaluations = 4;
  std::vector<ArgType> initial_population(5, ArgType{0.5, 0.5});
  params.initial_population = &initial_population;
  CHECK_THROW(validate_qlognehvi_parameters(params), std::invalid_argument);
  initial_population.resize(4);
  CHECK_NO_THROW(validate_qlognehvi_parameters(params));
  initial_population[3] = {0.5};
  CHECK_THROW(validate_qlognehvi_parameters(params), std::domain_error);
  initial_population[3] = {0.5, 1.5};
  CHECK_THROW(validate_qlognehvi_parameters(params), std::domain_error);
  initial_population[3] = {nan, 0.5};
  CHECK_THROW(validate_qlognehvi_parameters(params), std::domain_error);
  // And both qlognehvi and the optimizer validate:
  std::mt19937_64 gen(0);
  CHECK_THROW(qlognehvi([](ArgType const &x) { return x; }, params, gen), std::domain_error);
  using Optimizer = qlognehvi_optimizer<ArgType, ArgType>;
  CHECK_THROW(Optimizer(params, gen), std::domain_error);
}

int main() {
  test_mark_nondominated();
  test_acquisition_gradient();
  test_schaffer();
  test_zdt1();
  test_dtlz2();
  test_thread_count_independence();
  test_queries_and_cancellation();
  test_initial_population<float>();
  test_initial_population<double>();
  test_nan_and_objective_count();
  test_reference_point();
  test_validation();
  return boost::math::test::report_errors();
}
