/*
 * Copyright Nick Thompson, 2026
 * Use, modification and distribution are subject to the
 * Boost Software License, Version 1.0. (See accompanying file
 * LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
 */

#ifndef BOOST_MATH_BUILD_MODULE
#include <boost/math/optimization/detail/hypervolume.hpp>
#else
import boost.math;
#endif

#include "math_unit_test.hpp"
#ifndef BOOST_MATH_BUILD_MODULE
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <random>
#include <vector>
#endif
using boost::math::optimization::detail::box_decomposition;
using boost::math::optimization::detail::fatmin;
using boost::math::optimization::detail::hypervolume;
using boost::math::optimization::detail::hypervolume_improvement;
using boost::math::optimization::detail::log_fatplus;
using boost::math::optimization::detail::log_hypervolume_improvement;
using boost::math::optimization::detail::log_sum_exp;
using boost::math::optimization::detail::nondominated_boxes;

template <class Real> using Points = std::vector<std::vector<Real>>;

// The points which count: those less than the reference in every objective.
template <class Real> Points<Real> below(Points<Real> const &points, std::vector<Real> const &r) {
  Points<Real> kept;
  for (auto const &p : points) {
    bool keep = true;
    for (size_t k = 0; k < r.size(); ++k) {
      keep = keep && p[k] < r[k];
    }
    if (keep) {
      kept.push_back(p);
    }
  }
  return kept;
}

// The hypervolume by slicing along the last of the first m objectives: between two consecutive values of it, the
// cross-section is the hypervolume, in one dimension less, of the points at or below the slice.
template <class Real> Real sliced(Points<Real> points, std::vector<Real> const &r, size_t m) {
  std::sort(points.begin(), points.end(), [&](auto const &a, auto const &b) { return a[m - 1] < b[m - 1]; });
  Real volume = 0;
  for (size_t i = 0; i < points.size(); ++i) {
    const Real next = i + 1 < points.size() ? points[i + 1][m - 1] : r[m - 1];
    if (m == 1) {
      return r[0] - points[0][0];
    }
    if (next > points[i][m - 1]) {
      volume += (next - points[i][m - 1]) * sliced(Points<Real>(points.begin(), points.begin() + i + 1), r, m - 1);
    }
  }
  return volume;
}

template <class Real> Real sliced_hypervolume(Points<Real> const &points, std::vector<Real> const &r) {
  return sliced(below(points, r), r, r.size());
}

// The sum of the volumes of the parts of the boxes which are above lo.
template <class Real> Real clipped_volume(box_decomposition<Real> const &boxes, Real lo) {
  Real total = 0;
  for (size_t i = 0; i < boxes.lower.size(); i += boxes.objectives) {
    Real volume = 1;
    for (size_t k = 0; k < boxes.objectives; ++k) {
      volume *= (std::max)(Real(0), boxes.upper[i + k] - (std::max)(lo, boxes.lower[i + k]));
    }
    total += volume;
  }
  return total;
}

// No two boxes overlap, none is empty or reaches beyond the reference, and no point is strictly below the upper
// corner of a box, which is to say that no interior point of a box is dominated.
template <class Real>
bool boxes_are_sound(box_decomposition<Real> const &boxes, Points<Real> const &points, std::vector<Real> const &r) {
  const size_t m = r.size();
  bool sound = boxes.objectives == m && boxes.upper.size() == boxes.lower.size();
  sound = sound && boxes.size() * m == boxes.lower.size();
  for (size_t i = 0; i < boxes.size(); ++i) {
    for (size_t k = 0; k < m; ++k) {
      sound = sound && boxes.lower[i * m + k] < boxes.upper[i * m + k] && boxes.upper[i * m + k] <= r[k];
    }
    for (size_t j = i + 1; j < boxes.size(); ++j) {
      bool overlap = true;
      for (size_t k = 0; k < m; ++k) {
        overlap = overlap && (std::max)(boxes.lower[i * m + k], boxes.lower[j * m + k]) <
                               (std::min)(boxes.upper[i * m + k], boxes.upper[j * m + k]);
      }
      sound = sound && !overlap;
    }
    for (auto const &p : points) {
      bool inside = true;
      for (size_t k = 0; k < m; ++k) {
        inside = inside && p[k] < boxes.upper[i * m + k];
      }
      sound = sound && !inside;
    }
  }
  return sound;
}

// On the grid of all the coordinates of the points, the reference and lo: every corner of a box is on the grid, a
// cell which a point dominates is in no box and any other cell is in exactly one. Returns the volume of the dominated
// cells, which is the hypervolume, or -1 if the boxes are not such a partition.
template <class Real>
Real grid_hypervolume(box_decomposition<Real> const &boxes, Points<Real> const &points, std::vector<Real> const &r,
                      Real lo) {
  const size_t m = r.size();
  const Points<Real> kept = below(points, r);
  Points<Real> grid(m);
  for (size_t k = 0; k < m; ++k) {
    grid[k] = {lo, r[k]};
    for (auto const &p : kept) {
      grid[k].push_back(p[k]);
    }
    std::sort(grid[k].begin(), grid[k].end());
    grid[k].erase(std::unique(grid[k].begin(), grid[k].end()), grid[k].end());
  }
  bool partition = true;
  for (size_t i = 0; i < boxes.lower.size(); ++i) {
    auto const &g = grid[i % m];
    partition = partition && std::binary_search(g.begin(), g.end(), boxes.upper[i]) &&
                (boxes.lower[i] == -std::numeric_limits<Real>::infinity() ||
                 std::binary_search(g.begin(), g.end(), boxes.lower[i]));
  }
  Real volume = 0;
  std::vector<size_t> cell(m, 0);
  std::vector<Real> centre(m);
  while (cell[m - 1] + 1 < grid[m - 1].size()) {
    Real cell_volume = 1;
    for (size_t k = 0; k < m; ++k) {
      centre[k] = (grid[k][cell[k]] + grid[k][cell[k] + 1]) / 2;
      cell_volume *= grid[k][cell[k] + 1] - grid[k][cell[k]];
    }
    bool dominated = false;
    for (auto const &p : kept) {
      bool below_centre = true;
      for (size_t k = 0; k < m; ++k) {
        below_centre = below_centre && p[k] < centre[k];
      }
      dominated = dominated || below_centre;
    }
    size_t holders = 0;
    for (size_t i = 0; i < boxes.size(); ++i) {
      bool inside = true;
      for (size_t k = 0; k < m; ++k) {
        inside = inside && boxes.lower[i * m + k] < centre[k] && centre[k] < boxes.upper[i * m + k];
      }
      holders += inside;
    }
    partition = partition && holders == (dominated ? 0u : 1u);
    volume += dominated ? cell_volume : 0;
    // The next cell:
    for (size_t k = 0; k < m; ++k) {
      if (++cell[k] + 1 < grid[k].size() || k + 1 == m) {
        break;
      }
      cell[k] = 0;
    }
  }
  return partition ? volume : Real(-1);
}

// kind 0: uniform in the unit cube, so that some points are dominated and some are beyond the reference, 0.9.
// kind 1: integers from 0 to 7, so that there are ties and duplicates, and points on the reference, 6, and beyond it.
// kind 2: on the unit sphere, so that no point dominates another, each third point a copy of the one before, and
// reference 1.1.
template <class Real, class URBG>
Points<Real> random_points(size_t n, size_t m, int kind, URBG &gen, std::vector<Real> &r) {
  std::uniform_real_distribution<Real> uniform(0, 1);
  std::uniform_int_distribution<int> integer(0, 7);
  std::normal_distribution<Real> normal(0, 1);
  using std::abs;
  using std::sqrt;
  r.assign(m, kind == 0 ? Real(0.9) : (kind == 1 ? Real(6) : Real(1.1)));
  Points<Real> points(n, std::vector<Real>(m));
  for (size_t i = 0; i < n; ++i) {
    Real norm = 0;
    for (auto &x : points[i]) {
      x = kind == 0 ? uniform(gen) : (kind == 1 ? static_cast<Real>(integer(gen)) : abs(normal(gen)));
      norm += x * x;
    }
    for (auto &x : points[i]) {
      x = kind == 2 ? x / sqrt(norm) : x;
    }
    if (kind == 2 && i % 3 == 2) {
      points[i] = points[i - 1];
    }
  }
  return points;
}

template <class Real> void test_by_hand() {
  using A2 = std::array<Real, 2>;
  using A3 = std::array<Real, 3>;
  const Real nan = std::numeric_limits<Real>::quiet_NaN();
  const Real tol = 8 * std::numeric_limits<Real>::epsilon();
  // A staircase of three steps, of areas 1, 2 and 3:
  std::vector<A2> stairs = {{3, 1}, {1, 3}, {2, 2}};
  const A2 r2 = {4, 4};
  CHECK_ABSOLUTE_ERROR(Real(6), hypervolume(stairs, r2), 6 * tol);
  auto boxes = nondominated_boxes(stairs, r2);
  CHECK_EQUAL(size_t(4), boxes.size());
  CHECK_ABSOLUTE_ERROR(Real(10), hypervolume_improvement(A2{0, 0}, boxes), 16 * tol);
  CHECK_ABSOLUTE_ERROR(Real(0.25), hypervolume_improvement(A2{2.5, 1.5}, boxes), tol);
  CHECK_ABSOLUTE_ERROR(Real(15), hypervolume_improvement(A2{-1, 1}, boxes) + hypervolume(stairs, r2), 16 * tol);
  // A point of the front, a dominated point, a point on the reference, and a NaN add nothing:
  for (A2 y : {A2{2, 2}, A2{3, 3}, A2{4, 0}, A2{0, 4}, A2{5, 5}, A2{nan, 0}, A2{0, nan}}) {
    CHECK_EQUAL(Real(0), hypervolume_improvement(y, boxes));
  }
  // The same staircase with a duplicate, a dominated point, a point on the reference and a NaN:
  stairs.insert(stairs.end(), {{2, 2}, {3, 3}, {0, 4}, {4, 0}, {nan, 0}, {0, nan}});
  CHECK_ABSOLUTE_ERROR(Real(6), hypervolume(stairs, r2), 6 * tol);
  CHECK_EQUAL(size_t(4), nondominated_boxes(stairs, r2).size());

  // No points: nothing is dominated, and everything below the reference is one box.
  CHECK_EQUAL(Real(0), hypervolume(std::vector<A3>{}, A3{1, 2, 3}));
  CHECK_EQUAL(Real(0), hypervolume(std::vector<A3>{{1, 2, 3}, {nan, 0, 0}, {0, 0, 4}}, A3{1, 2, 3}));
  auto empty = nondominated_boxes(std::vector<A3>{}, A3{1, 2, 3});
  CHECK_EQUAL(size_t(1), empty.size());
  CHECK_ABSOLUTE_ERROR(Real(0.75), hypervolume_improvement(A3{0.5, 1.5, 0}, empty), tol);
  CHECK_ABSOLUTE_ERROR(Real(1), hypervolume(std::vector<A3>{{1, 1, 1}}, A3{2, 2, 2}), tol);

  // Three boxes of volume 6, which meet in pairs in volume 2, and all together in volume 1: 18 - 6 + 1.
  std::vector<A3> three = {{1, 2, 3}, {2, 3, 1}, {3, 1, 2}};
  const A3 r3 = {4, 4, 4};
  CHECK_ABSOLUTE_ERROR(Real(13), hypervolume(three, r3), 13 * tol);
  auto boxes3 = nondominated_boxes(three, r3);
  CHECK_ABSOLUTE_ERROR(Real(51), hypervolume_improvement(A3{0, 0, 0}, boxes3), 64 * tol);
  // The cube above (2, 2, 2) has volume 8, of which 12 - 6 + 1 is dominated:
  CHECK_ABSOLUTE_ERROR(Real(1), hypervolume_improvement(A3{2, 2, 2}, boxes3), 8 * tol);
  three.insert(three.end(), {{3, 3, 3}, {1, 2, 3}, {0, 0, 4}, {nan, nan, nan}});
  CHECK_ABSOLUTE_ERROR(Real(13), hypervolume(three, r3), 13 * tol);

  // One objective:
  using A1 = std::array<Real, 1>;
  CHECK_ABSOLUTE_ERROR(Real(3), hypervolume(std::vector<A1>{{2}, {1}, {5}}, A1{4}), tol);
  const auto boxes1 = nondominated_boxes(std::vector<A1>{{2}, {1}, {5}}, A1{4});
  CHECK_ABSOLUTE_ERROR(Real(0.5), hypervolume_improvement(A1{0.5}, boxes1), tol);
}

// Small sets, so that every cell of the grid of their coordinates can be visited.
template <class Real> void test_partition() {
  std::mt19937_64 gen(271828);
  const Real tol = 500 * std::numeric_limits<Real>::epsilon();
  const size_t sizes[] = {0, 0, 40, 25, 12, 8, 6};
  for (size_t m = 2; m <= 6; ++m) {
    for (int kind = 0; kind < 3; ++kind) {
      for (size_t trial = 0; trial < 4; ++trial) {
        std::vector<Real> r;
        // The first trial of each is of very few points:
        auto points = random_points(trial == 0 ? trial + size_t(kind) : sizes[m], m, kind, gen, r);
        auto boxes = nondominated_boxes(points, r);
        const Real lo = -0.5;
        Real bounding = 1;
        for (size_t k = 0; k < m; ++k) {
          bounding *= r[k] - lo;
        }
        const Real expected = grid_hypervolume(boxes, points, r, lo);
        CHECK_GE(expected, Real(0));
        CHECK_TRUE(boxes_are_sound(boxes, below(points, r), r));
        CHECK_ABSOLUTE_ERROR(expected, hypervolume(points, r), tol * bounding);
        CHECK_ABSOLUTE_ERROR(expected, sliced_hypervolume(points, r), tol * bounding);
        CHECK_ABSOLUTE_ERROR(bounding, expected + clipped_volume(boxes, lo), tol * bounding);
      }
    }
  }
}

// Larger sets, against the hypervolume by slicing, and the improvements against differences of hypervolumes.
template <class Real> void test_random() {
  std::mt19937_64 gen(314159);
  const Real tol = 2000 * std::numeric_limits<Real>::epsilon();
  const size_t sizes[] = {0, 0, 60, 40, 30, 24, 14};
  for (size_t m = 2; m <= 6; ++m) {
    for (int kind = 0; kind < 3; ++kind) {
      std::vector<Real> r;
      auto points = random_points(sizes[m], m, kind, gen, r);
      points[sizes[m] / 2][m - 1] = std::numeric_limits<Real>::quiet_NaN();
      auto boxes = nondominated_boxes(points, r);
      CHECK_TRUE(boxes_are_sound(boxes, below(points, r), r));
      const Real lo = kind == 1 ? Real(-2.5) : Real(-0.25);
      Real bounding = 1;
      for (size_t k = 0; k < m; ++k) {
        bounding *= r[k] - lo;
      }
      const Real volume = hypervolume(points, r);
      CHECK_ABSOLUTE_ERROR(sliced_hypervolume(points, r), volume, tol * bounding);
      CHECK_ABSOLUTE_ERROR(bounding, volume + clipped_volume(boxes, lo), tol * bounding);

      // Candidates from the same distribution, moved towards the origin so that more of them improve on the front,
      // some of them copies of members of the front:
      std::vector<Real> unused;
      auto candidates = random_points(48, m, kind, gen, unused);
      for (auto &y : candidates) {
        for (auto &x : y) {
          x = kind == 1 ? x - 2 : x * Real(0.6);
        }
      }
      candidates[5] = points[0];
      candidates[11] = points[1];
      size_t improving = 0;
      auto joined = points;
      for (auto const &y : candidates) {
        joined.push_back(y);
        const Real improvement = hypervolume_improvement(y, boxes);
        improving += improvement > 0;
        CHECK_ABSOLUTE_ERROR(sliced_hypervolume(joined, r) - volume, improvement, tol * bounding);
        joined.pop_back();
      }
      CHECK_GE(improving, size_t(2));
    }
  }
}

// n*m consecutive values of the linear congruential generator x -> (1103515245 x + 12345) mod 2^31 from the seed,
// divided by 2^31.
Points<double> congruential_points(size_t n, size_t m, std::uint64_t x) {
  Points<double> points(n, std::vector<double>(m));
  for (auto &p : points) {
    for (auto &y : p) {
      x = (1103515245u * x + 12345u) % 2147483648u;
      y = static_cast<double>(x) / 2147483648.0;
    }
  }
  return points;
}

// Against BoTorch 0.18.1, which maximizes: FastNondominatedPartitioning(ref_point=-r, Y=-Y).compute_hypervolume(),
// and the number of cells of its get_hypercell_bounds(), which ours does not exceed. The points are
// congruential_points(n, m, seed); r is 0.9 in each objective.
void test_botorch() {
  struct reference { size_t m; size_t n; std::uint64_t seed; double hypervolume; size_t boxes; };
  const reference cases[] = {{2, 30, 1, 0.60925675464725071, 4},  {3, 30, 2, 0.66025243996159855, 10},
                             {4, 30, 3, 0.3144074518964215, 64},  {5, 30, 4, 0.16364217902248246, 339},
                             {5, 60, 5, 0.3457065676990364, 170}, {6, 20, 6, 0.028019035081987942, 353}};
  for (auto const &c : cases) {
    const auto points = congruential_points(c.n, c.m, c.seed);
    const std::vector<double> r(c.m, 0.9);
    CHECK_ABSOLUTE_ERROR(c.hypervolume, hypervolume(points, r), 1e-14);
    CHECK_LE(nondominated_boxes(points, r).size(), c.boxes);
  }
}

template <class Real> void test_smooth() {
  using std::abs;
  using std::exp;
  using std::log;
  using std::sqrt;
  const Real eps = std::numeric_limits<Real>::epsilon();
  const Real tol = 16 * eps;
  Real d = 0, d_minus = 0, d_plus = 0;
  for (Real tau : {Real(1e-6), Real(1e-2), Real(4)}) {
    // log(x) far above zero, log(tau (log(2) + 1/10)) at zero, and log(tau^3/(10 x^2)) far below it:
    CHECK_ABSOLUTE_ERROR(log(Real(1e9) * tau), log_fatplus(Real(1e9) * tau, tau, d), 2 * tol);
    CHECK_ABSOLUTE_ERROR(Real(1), d * Real(1e9) * tau, tol);
    CHECK_ABSOLUTE_ERROR(log(tau * (log(Real(2)) + Real(0.1))), log_fatplus(Real(0), tau, d), 16 * tol);
    CHECK_ABSOLUTE_ERROR(Real(0.5), d * tau * (log(Real(2)) + Real(0.1)), tol);
    for (Real x : {Real(-1e3), Real(-1e9), Real(-1e15), -(std::numeric_limits<Real>::max)() / 8}) {
      x *= tau;
      CHECK_ABSOLUTE_ERROR(log(tau / 10) - 2 * log(-x / tau), log_fatplus(x, tau, d), Real(1e-5));
      CHECK_ABSOLUTE_ERROR(Real(1), d * x / -2, Real(1e-5));
    }
    CHECK_TRUE(std::isfinite(log_fatplus(tau * ((std::numeric_limits<Real>::max)() / 8), tau, d)) && d > 0);
    // It increases, and the derivative is that of the value, across zero and across both changes of formula:
    const Real h = std::numeric_limits<Real>::digits > 24 ? Real(1e-5) : Real(1e-2);
    for (Real t = -3 / sqrt(eps); t < 50; t = t < -60 ? t / 2 : t + Real(0.37)) {
      const Real value = log_fatplus(t * tau, tau, d), step = h * (1 + abs(t));
      const Real difference = log_fatplus((t + step) * tau, tau, d_plus) - log_fatplus((t - step) * tau, tau, d_minus);
      CHECK_TRUE(d > 0 && difference > 0);
      CHECK_ABSOLUTE_ERROR(d * tau * 2 * step, difference, 100 * eps * (1 + abs(value)) + h * difference);
    }
  }
  // fatmin is the lesser far apart, and the mean less tau log(2), with derivative 1/2, where the two are equal:
  for (Real tau : {Real(1e-3), Real(1e-2), Real(1)}) {
    CHECK_ABSOLUTE_ERROR(Real(3) - tau * log(Real(2)), fatmin(Real(3), Real(3), tau, d), tol);
    CHECK_ABSOLUTE_ERROR(Real(0.5), d, tol);
    CHECK_ABSOLUTE_ERROR(Real(-2), fatmin(Real(-2), Real(3e9), tau, d), tol);
    CHECK_ABSOLUTE_ERROR(Real(1), d, tol);
    CHECK_ABSOLUTE_ERROR(Real(-2), fatmin(Real(3e9), Real(-2), tau, d), tol);
    CHECK_TRUE(d >= 0 && d < tol);
    const Real h = std::numeric_limits<Real>::digits > 24 ? Real(1e-5) : Real(1e-2);
    for (Real a = -40 * tau; a < 40 * tau; a += Real(0.73) * tau) {
      const Real value = fatmin(a, Real(0.01) * tau, tau, d);
      CHECK_LE(value, (std::min)(a, Real(0.01) * tau));
      const Real difference = fatmin(a + h * tau, Real(0.01) * tau, tau, d_plus) -
                              fatmin(a - h * tau, Real(0.01) * tau, tau, d_minus);
      CHECK_ABSOLUTE_ERROR(d * 2 * h * tau, difference, tau * (200 * eps + h * h * h));
      CHECK_ABSOLUTE_ERROR(value, fatmin(Real(0.01) * tau, a, tau, d_minus), 100 * eps * tau);
      CHECK_ABSOLUTE_ERROR(Real(1), d + d_minus, tol);
    }
  }
  // Terms which would overflow or underflow; the rescaled weights are exp(x_i - largest):
  for (Real shift : {Real(0), Real(-1e4), Real(1e4)}) {
    log_sum_exp<Real> sum;
    Real first = 0, rescale = 0;
    for (Real x : {Real(-3), Real(-2), Real(-40), Real(1), Real(1), Real(0)}) {
      const Real weight = sum.add(x + shift, rescale);
      first = x == -3 ? weight : first * rescale;
    }
    const Real expected = exp(Real(-3)) + exp(Real(-2)) + exp(Real(-40)) + 2 * exp(Real(1)) + 1;
    CHECK_ABSOLUTE_ERROR(shift + log(expected), sum.value(), tol * (1 + abs(shift)));
    CHECK_ABSOLUTE_ERROR(exp(Real(-3)) / expected, first / sum.sum, tol);
  }
}

// The improvement itself: each side u - max(l, y) is min(max(u - y, 0), u - l), which the smooth one exceeds by less
// than (log(2) + 1/10) tau_relu, before fatmin takes off a factor of at most 2^tau_max.
template <class Real> void test_log_improvement() {
  using std::exp;
  using std::log;
  using A2 = std::array<Real, 2>;
  const Real eps = std::numeric_limits<Real>::epsilon();
  std::mt19937_64 gen(161803);
  std::uniform_real_distribution<Real> uniform(-0.1, 1);
  for (size_t m = 2; m <= 5; ++m) {
    for (int kind = 0; kind < 3; kind += 2) {
      std::vector<Real> r, y(m), gradient(m);
      const auto boxes = nondominated_boxes(random_points(40, m, kind, gen, r), r);
      size_t improving = 0;
      for (size_t trial = 0; trial < 60; ++trial) {
        for (auto &x : y) {
          x = uniform(gen);
        }
        const Real tau_relu = trial % 2 ? Real(1e-6) : Real(1e-4);
        const Real tau_max = trial % 3 ? Real(1e-2) : Real(1e-5);
        const Real exact = hypervolume_improvement(y, boxes);
        Real most = 0;
        for (size_t i = 0; i < boxes.lower.size(); i += m) {
          Real volume = 1;
          for (size_t k = 0; k < m; ++k) {
            volume *= (std::min)((std::max)(boxes.upper[i + k] - y[k], Real(0)) + Real(0.8) * tau_relu,
                                 boxes.upper[i + k] - (std::max)(boxes.lower[i + k], Real(-1e8)));
          }
          most += volume;
        }
        const Real value = log_hypervolume_improvement(y.data(), boxes, gradient.data(), tau_relu, tau_max);
        CHECK_ULP_CLOSE(value, log_hypervolume_improvement(y.data(), boxes, nullptr, tau_relu, tau_max), 4);
        CHECK_LE(value, log(most) + 100 * eps);
        CHECK_GE(value, log(exact) - m * tau_max * log(Real(2)) - 100 * eps);
        if (exact > Real(1e-3)) {
          // The bounds are then close to one another:
          ++improving;
          CHECK_ABSOLUTE_ERROR(log(exact), value, m * tau_max * log(Real(2)) + boxes.size() * tau_relu / exact);
        }
        for (auto g : gradient) {
          CHECK_TRUE(std::isfinite(g) && g <= 0);
        }
      }
      CHECK_GE(improving, size_t(5));
    }
  }

  // Where nothing improves, on the front, behind it, on the reference and beyond it, the value is finite and small,
  // less than log(4 boxes * 0.8 tau_relu * 4), and since less is better, the gradient is negative:
  const auto stairs = nondominated_boxes(std::vector<A2>{{3, 1}, {1, 3}, {2, 2}}, A2{4, 4});
  A2 gradient;
  for (A2 y : {A2{2, 2}, A2{3, 1}, A2{3, 3}, A2{4, 0}, A2{0, 4}, A2{5, 5}, A2{Real(2.5), 2}}) {
    CHECK_EQUAL(Real(0), hypervolume_improvement(y, stairs));
    const Real value = log_hypervolume_improvement(y.data(), stairs, gradient.data());
    CHECK_TRUE(std::isfinite(value) && value < -11);
    CHECK_TRUE(std::isfinite(gradient[0]) && gradient[0] < 0 && std::isfinite(gradient[1]) && gradient[1] < 0);
  }
  // It decreases all the way from far below the front, across it, to far beyond the reference:
  Real previous = std::numeric_limits<Real>::infinity();
  for (Real s = -1e6; s < 1e6; s = s < -1 || s > 8 ? (s < 0 ? s / 2 : s * 2) : s + Real(0.0625)) {
    const A2 y = {s, s / 2 + 1};
    const Real value = log_hypervolume_improvement(y.data(), stairs, gradient.data());
    CHECK_TRUE(std::isfinite(value) && value < previous);
    CHECK_TRUE(gradient[0] < 0 && gradient[1] < 0);
    previous = value;
  }
  CHECK_LE(previous, Real(-100));
}

// The gradient against fourth order central differences in double precision, and the single precision gradient
// against that: near corners of the boxes, where the scale of the function is less than tau_relu, among them, and
// far away. The error is that of the differences, which a smaller step reduces there, at the expense of the small
// derivatives far away.
void test_log_improvement_gradient() {
  std::mt19937_64 gen(141421);
  std::uniform_real_distribution<double> uniform(-1, 1);
  for (size_t m = 2; m <= 5; ++m) {
    for (int kind = 0; kind < 3; kind += 2) {
      std::vector<double> r, y(m), gradient(m);
      const auto points = random_points(kind && m == 2 ? 0 : 30, m, kind, gen, r);
      const auto boxes = nondominated_boxes(points, r);
      box_decomposition<float> boxes_float;
      boxes_float.objectives = m;
      boxes_float.lower.assign(boxes.lower.begin(), boxes.lower.end());
      boxes_float.upper.assign(boxes.upper.begin(), boxes.upper.end());
      box_decomposition<double> rounded;
      rounded.objectives = m;
      rounded.lower.assign(boxes_float.lower.begin(), boxes_float.lower.end());
      rounded.upper.assign(boxes_float.upper.begin(), boxes_float.upper.end());
      for (size_t trial = 0; trial < 200; ++trial) {
        // The distance from a corner of a box, or from the origin:
        const double scale = std::pow(10.0, trial % 4 == 3 ? 3 * uniform(gen) : -4 + 4 * uniform(gen));
        const size_t box = gen() % boxes.size();
        for (size_t k = 0; k < m; ++k) {
          const double lower = boxes.lower[box * m + k], upper = boxes.upper[box * m + k];
          y[k] = scale * uniform(gen) + (trial % 4 == 3 ? 0 : (trial % 4 == 2 || lower < 0 ? upper : lower));
        }
        const double tau_relu = trial % 5 ? 1e-6 : 1e-3, tau_max = trial % 7 ? 1e-2 : 1e-3;
        const double value = log_hypervolume_improvement(y.data(), boxes, gradient.data(), tau_relu, tau_max);
        CHECK_TRUE(std::isfinite(value));
        for (size_t k = 0; k < m; ++k) {
          const double yk = y[k], h = 1e-9 * (std::max)(1.0, 100 * std::abs(yk));
          double difference = 0;
          for (double step : {-2.0, -1.0, 1.0, 2.0}) {
            y[k] = yk + step * h;
            difference += (std::abs(step) > 1 ? -step / 24 : step * 2 / 3) *
                          log_hypervolume_improvement(y.data(), boxes, nullptr, tau_relu, tau_max);
          }
          y[k] = yk;
          const double error = std::abs(difference / h - gradient[k]) / (std::abs(gradient[k]) + 1e-3);
          CHECK_LE(error, 1e-3);
        }
        std::vector<float> y_float(y.begin(), y.end()), gradient_float(m);
        y.assign(y_float.begin(), y_float.end());
        const float tau_relu_float = static_cast<float>(tau_relu), tau_max_float = static_cast<float>(tau_max);
        const double value_float =
          log_hypervolume_improvement(y_float.data(), boxes_float, gradient_float.data(), tau_relu_float,
                                      tau_max_float);
        CHECK_ABSOLUTE_ERROR(log_hypervolume_improvement(y.data(), rounded, gradient.data(), tau_relu, tau_max),
                             value_float, 1e-3);
        for (size_t k = 0; k < m; ++k) {
          const double error = std::abs(gradient_float[k] - gradient[k]) / (std::abs(gradient[k]) + 1e-3);
          CHECK_LE(error, 1e-2);
        }
      }
    }
  }
}

// Nothing overflows, underflows to zero or is NaN, whatever the candidate and the temperatures.
template <class Real> void test_log_improvement_extremes() {
  std::mt19937_64 gen(173205);
  std::uniform_real_distribution<Real> uniform(-1, 1);
  using std::pow;
  const Real taus[][2] = {{Real(1e-6), Real(1e-2)}, {Real(1e-9), Real(1e-4)}, {Real(1e-3), Real(1)}, {1, 1}};
  for (size_t m = 1; m <= 9; ++m) {
    for (int kind = 0; kind < 3; ++kind) {
      std::vector<Real> r, y(m), gradient(m);
      // With no points there is a single box, without lower corner:
      const auto boxes = nondominated_boxes(random_points(kind == 2 || m > 5 ? 0 : 30, m, kind, gen, r), r);
      CHECK_TRUE((kind != 2 && m <= 5) || boxes.size() == 1);
      for (size_t trial = 0; trial < 400; ++trial) {
        for (auto &x : y) {
          x = (trial % 3 ? uniform(gen) : Real(1)) * pow(Real(10), 6 * uniform(gen));
          x = trial % 8 == 0 ? Real(trial % 16 ? 1e6 : -1e6) : x;
        }
        auto const &tau = taus[trial % 4];
        const Real value = log_hypervolume_improvement(y.data(), boxes, gradient.data(), tau[0], tau[1]);
        CHECK_TRUE(std::isfinite(value));
        for (auto g : gradient) {
          CHECK_TRUE(std::isfinite(g) && g <= 0);
        }
      }
    }
  }
}

// Against BoTorch 0.18.1, default tau_relu = 1e-6, tau_max = 1e-2 and fat = True: the points are
// congruential_points(30, m, m - 1), r is 0.9, and the value is
//   qLogExpectedHypervolumeImprovement(model, ref_point=-r, partitioning=FastNondominatedPartitioning(-r, -Y))
//     ._compute_log_qehvi(samples=-y.view(1, 1, 1, m)),
// by its fused kernel and without it, and gradient0 is minus the derivative of that with respect to -y[0]. The
// candidates improve by much, improve on one point of the front by 0.01, are that point, are behind another, on the
// reference, beyond it, and far below the front.
void test_botorch_log_improvement() {
  struct reference { size_t m; double y[4]; double value; double gradient0; };
  const reference cases[] = {
    {2, {0.05, 0.05}, -2.178256993816055, -7.5082636171907895},
    {2, {0.07389880321919919, 0.37964712107554077}, -4.942631615543717, -72.94321043067649},
    {2, {0.08389880321919918, 0.3896471210755408}, -14.400588340316105, -458075.7420023379},
    {2, {0.2972257970646024, 0.3880707192979753}, -39.4272435316422, -39.82240860855442},
    {2, {0.9, 0.1}, -42.74008143038837, -113.31211050820197},
    {2, {1.5, 0.1}, -43.99131958271342, -1.5363485235275507},
    {2, {-3.0, -3.0}, 2.681072348397625, -0.26710973086812906},
    {3, {0.05, 0.05, 0.05}, -4.780847609522647, -49.948744117147456},
    {3, {-0.002103490047156811, 0.10526243900880218, 0.6121016384661198}, -5.862066475064785, -80.12489376193763},
    {3, {0.00789650995284319, 0.11526243900880218, 0.6221016384661198}, -15.353506946939909, -500618.1213945685},
    {3, {0.03114998897537589, 0.8359149736166, 0.029696674067527056}, -38.05584386672447, -44.53870290880194},
    {3, {0.9, 0.1, 0.1}, -43.88573121079869, -147.3872813139434},
    {3, {1.5, 0.1, 0.1}, -44.94910427947569, -1.3738055097665296},
    {3, {-3.0, -3.0, -3.0}, 4.071736529487127, -0.2592963621620229},
    {4, {0.05, 0.05, 0.05, 0.05}, -1.50653267191317, -2.7697654059040455},
    {4, {0.07166811242699624, 0.0835700391419232, 0.6583647064492106, 0.7766595122031867}, -7.247601136672095,
     -34.066011585353365},
    {4, {0.08166811242699623, 0.09357003914192319, 0.6683647064492106, 0.7866595122031868}, -16.79303692404261,
     -203623.96235112107},
    {4, {0.15730478845536708, 0.292003332618624, 0.25546307254582645, 0.2088661808334291}, -36.17385986662908,
     -34.30490955552053},
    {4, {0.9, 0.1, 0.1, 0.1}, -16.16879975265499, -630400.021856345},
    {4, {1.5, 0.1, 0.1, 0.1}, -44.31774863586113, -2.610837760661149},
    {4, {-3.0, -3.0, -3.0, -3.0}, 5.442546026289447, -0.2567592037031373}};
  // The cells of its partitioning for two and three objectives, negated: upper is -get_hypercell_bounds()[0].
  const double inf = std::numeric_limits<double>::infinity();
  box_decomposition<double> cells[4];
  cells[2].lower = {0.3191331308335066, -inf, 0.2772257970646024, -inf, 0.08389880321919918, -inf, -inf, -inf};
  cells[2].upper = {0.9, 0.059359821025282145, 0.3191331308335066, 0.3680707192979753,
                    0.2772257970646024, 0.3896471210755408, 0.08389880321919918, 0.9};
  cells[3].lower = {-inf, -inf, -inf, -inf, 0.11526243900880218, 0.6221016384661198,
                    -inf, 0.021475194953382015, -inf, -inf, 0.021475194953382015, 0.020409866701811552,
                    -inf, 0.11526243900880218, 0.020409866701811552, -inf, 0.6963285580277443, 0.020409866701811552,
                    -inf, 0.6963285580277443, 0.31248870911076665, -inf, 0.8159149736166, 0.31248870911076665,
                    -inf, 0.8159149736166, -inf, -inf, 0.8159149736166, 0.009696674067527056};
  cells[3].upper = {0.9, 0.021475194953382015, 0.9, 0.00789650995284319, 0.9, 0.9,
                    0.9, 0.8159149736166, 0.020409866701811552, 0.07002052338793874, 0.11526243900880218, 0.9,
                    0.07002052338793874, 0.6963285580277443, 0.6221016384661198,
                    0.07002052338793874, 0.8159149736166, 0.31248870911076665,
                    0.02773440768942237, 0.8159149736166, 0.6221016384661198,
                    0.01114998897537589, 0.9, 0.6221016384661198, 0.9, 0.9, 0.009696674067527056,
                    0.01114998897537589, 0.9, 0.31248870911076665};
  for (auto const &c : cases) {
    double gradient[4];
    if (c.m < 4) {
      // The same boxes give the same value:
      cells[c.m].objectives = c.m;
      CHECK_ABSOLUTE_ERROR(c.value, log_hypervolume_improvement(c.y, cells[c.m], gradient), 1e-13);
      CHECK_ABSOLUTE_ERROR(1.0, gradient[0] / c.gradient0, 1e-13);
    }
    // Our boxes are not BoTorch's, and the value depends on them: little where y improves on the front, since
    // both are then within m tau_max log(2) of the logarithm of the improvement, but elsewhere it is a sum over the
    // boxes of the tails of log_fatplus, which changes by a factor of order one with the decomposition. So it does
    // between BoTorch's FastNondominatedPartitioning and NondominatedPartitioning.
    const auto boxes = nondominated_boxes(congruential_points(30, c.m, c.m - 1), std::vector<double>(c.m, 0.9));
    const bool improving = hypervolume_improvement(std::vector<double>(c.y, c.y + c.m), boxes) > 0;
    const double difference = std::abs(log_hypervolume_improvement(c.y, boxes) - c.value);
    CHECK_LE(difference, improving ? c.m * 1e-2 * std::log(2.0) : 1.0);
  }
}

int main() {
  test_by_hand<float>();
  test_by_hand<double>();
  test_partition<float>();
  test_partition<double>();
  test_random<float>();
  test_random<double>();
  test_botorch();
  test_smooth<float>();
  test_smooth<double>();
  test_log_improvement<float>();
  test_log_improvement<double>();
  test_log_improvement_gradient();
  test_log_improvement_extremes<float>();
  test_log_improvement_extremes<double>();
  test_botorch_log_improvement();
  return boost::math::test::report_errors();
}
