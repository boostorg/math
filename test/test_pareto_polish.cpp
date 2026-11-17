/*
 * Copyright Nick Thompson, 2026
 * Use, modification and distribution are subject to the
 * Boost Software License, Version 1.0. (See accompanying file
 * LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
 */

#include <boost/math/optimization/pareto_polish.hpp>

#include "math_unit_test.hpp"
#include "test_functions_for_optimization.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <mutex>
#include <random>
#include <set>
#include <stdexcept>
#include <thread>
#include <utility>
#include <vector>
using namespace boost::math::optimization;
using V = std::vector<double>;
using cell_class = pareto_polish_cell_class;
using status = pareto_polish_status;
using termination = pareto_polish_termination;
using event_kind = pareto_polish_event_kind;

// Where there are no threads, as in WebAssembly without pthreads, the same paths run on the calling thread.
// So they do on MinGW, where destroying the thread_local reverse mode tape of a worker corrupts the heap.
#if (defined(__EMSCRIPTEN__) && !defined(__EMSCRIPTEN_PTHREADS__)) || defined(BOOST_MATH_NO_THREAD_LOCAL_WITH_NON_TRIVIAL_TYPES)
const unsigned several_threads = 1;
#else
const unsigned several_threads = 3;
#endif

// In [0, 1), and the same everywhere: the distributions of the standard library differ from one to another.
double unit(std::mt19937_64 &gen) { return static_cast<double>(gen() >> 11) * 0x1p-53; }

// The subproblem is exact: its primal value equals its dual value at the multipliers returned, the step is
// feasible, and no other multipliers give a larger dual value.
template <class Real> void test_subproblem(size_t m, bool cone) {
  using std::abs;
  const Real inf = std::numeric_limits<Real>::infinity();
  std::mt19937_64 gen(1234 + 10 * m + cone);
  auto uniform = [&gen]() { return static_cast<Real>(unit(gen)); };
  Real worst_gap = 0, worst_feasibility = 0, worst_sample = 0;
  for (int t = 0; t < 77; ++t) {
    const size_t n = 1 + static_cast<size_t>(t) % 7;
    std::vector<Real> J(m * n), c(m, 0), lo(n), hi(n), lam, d;
    for (auto &v : J) {
      v = 4 * uniform() - 2 + (t % 3 == 0 ? 1 : 0);
    }
    if (t % 7 == 0) {
      std::copy(J.begin(), J.begin() + n, J.end() - n); // two equal rows
    }
    if (t % 11 == 0) {
      std::fill(J.begin(), J.begin() + n, Real(0)); // an objective with no gradient
    }
    for (size_t k = 0; k < n; ++k) {
      lo[k] = t % 5 == 0 ? 0 : (t % 13 == 0 ? -inf : -uniform());
      hi[k] = t % 4 == 0 ? inf : (t % 9 == 0 ? 0 : uniform());
    }
    for (size_t i = 0; i < m && t % 2; ++i) {
      c[i] = -uniform();
    }
    const Real t1 = uniform(), t2 = uniform(), tau = Real(1) / 10 + 10 * t1 * t2;
    detail::pareto_polish_subproblem(J, m, c, lo, hi, tau, cone, lam, d);
    auto dual = [&](std::vector<Real> const &l) {
      Real q = 0;
      for (size_t i = 0; i < m && !cone; ++i) {
        q += c[i] * l[i];
      }
      for (size_t k = 0; k < n; ++k) {
        Real w = 0;
        for (size_t i = 0; i < m; ++i) {
          w += l[i] * J[i * n + k];
        }
        const Real u = (std::min)((std::max)(-tau * w, lo[k]), hi[k]);
        q += w * u + u * u / (2 * tau);
      }
      return q;
    };
    Real primal = cone ? 0 : -inf, sum = 0, norm2 = 0, scale = 1;
    for (size_t k = 0; k < n; ++k) {
      norm2 += d[k] * d[k];
      worst_feasibility = (std::max)(worst_feasibility, (std::max)(lo[k] - d[k], d[k] - hi[k]));
    }
    for (size_t i = 0; i < m; ++i) {
      Real gd = 0;
      for (size_t k = 0; k < n; ++k) {
        gd += J[i * n + k] * d[k];
        scale += abs(J[i * n + k] * d[k]);
      }
      sum += lam[i];
      if (cone) {
        primal += gd;
        worst_feasibility = (std::max)(worst_feasibility, (std::max)(gd / scale, 1 - lam[i]));
      } else {
        primal = (std::max)(primal, c[i] + gd);
        worst_feasibility = (std::max)(worst_feasibility, -lam[i]);
      }
    }
    if (!cone) {
      worst_feasibility = (std::max)(worst_feasibility, abs(sum - 1));
    }
    primal += norm2 / (2 * tau);
    worst_gap = (std::max)(worst_gap, abs(primal - dual(lam)) / scale);
    std::vector<Real> l(m);
    for (int sample = 0; sample < 50; ++sample) {
      Real total = 0;
      for (auto &v : l) {
        v = uniform() * uniform() + std::numeric_limits<Real>::min();
        total += v;
      }
      for (auto &v : l) {
        v = cone ? 1 + 3 * v : v / total;
      }
      worst_sample = (std::max)(worst_sample, (dual(l) - primal) / scale);
    }
  }
  // Measured: 1e-11 for the simplex and 3e-9 for the cone in double precision.
  const Real tol = 100 * std::sqrt(std::numeric_limits<Real>::epsilon());
  CHECK_LE(worst_gap, tol);
  CHECK_LE(worst_feasibility, tol);
  CHECK_LE(worst_sample, tol);
}

// Squared distances to three sites: the Pareto set is their triangle, in the interior of the box.
const double site[3][4] = {{0.2, 0.3, 0.7, 0.4}, {0.8, 0.2, 0.5, 0.6}, {0.5, 0.9, 0.3, 0.5}};
struct sites {
  template <class C> std::vector<typename C::value_type> operator()(C const &x) const {
    std::vector<typename C::value_type> y;
    for (auto const &s : site) {
      typename C::value_type acc = (x[0] - s[0]) * (x[0] - s[0]);
      for (size_t k = 1; k < 4; ++k) {
        acc += (x[k] - s[k]) * (x[k] - s[k]);
      }
      y.push_back(acc);
    }
    return y;
  }
};

// Points of the box [0.05, 0.95]^n, the same everywhere.
template <class Real, class F> auto scattered_front(F f, size_t n, size_t count, unsigned seed) {
  std::mt19937_64 gen(seed);
  std::vector<std::pair<std::vector<Real>, decltype(f(std::vector<Real>()))>> front;
  for (size_t i = 0; i < count; ++i) {
    std::vector<Real> x(n);
    for (auto &v : x) {
      v = static_cast<Real>(0.05 + 0.9 * unit(gen));
    }
    front.emplace_back(x, f(x));
  }
  return front;
}

template <class V> pareto_polish_parameters<V> unit_box(size_t n, typename V::value_type refinement_tolerance) {
  pareto_polish_parameters<V> params;
  params.lower_bounds = V(n, 0);
  params.upper_bounds = V(n, 1);
  params.refinement_tolerance = refinement_tolerance;
  return params;
}

// The front of f at the points given.
template <class F, class X = V> auto front_at(F f, std::initializer_list<X> points) {
  std::vector<std::pair<X, decltype(f(X()))>> front;
  for (auto const &x : points) {
    front.emplace_back(x, f(x));
  }
  return front;
}

// The same points, to the last bit.
template <class Result> void check_same(Result const &p, Result const &q) {
  CHECK_EQUAL(p.size(), q.size());
  for (size_t i = 0; i < p.size() && i < q.size(); ++i) {
    CHECK_TRUE(p[i].x == q[i].x && p[i].f == q[i].f);
  }
}

// What holds of every result. Every point of the front with finite objectives is weakly dominated by a point of
// the result, one whose objectives are not finite is returned as it is, and no point of the result dominates
// another. The ranges are those of the points returned, a point which is said to be stationary is, and the
// counts, the achieved tolerance and the status agree with the cells.
template <class Front, class Result, class Params>
void check_guarantee(Front const &front, Result const &result, Params const &params) {
  using std::isfinite;
  using Real = typename Params::Real;
  auto no_worse = [](auto const &a, auto const &b) {
    bool ok = true;
    auto q = std::begin(b);
    for (auto const &v : a) {
      ok = ok && v <= *q;
      ++q;
    }
    return ok;
  };
  for (size_t i = 0; i < front.size(); ++i) {
    bool covered = false;
    for (auto const &r : result) {
      covered = covered || no_worse(r.f, front[i].second) ||
                (r.origin == i && r.termination == termination::non_finite && r.x == front[i].first);
    }
    CHECK_TRUE(covered);
  }
  const size_t m = result.ranges.size();
  std::vector<Real> low(m, std::numeric_limits<Real>::infinity()), high(m, -std::numeric_limits<Real>::infinity());
  size_t created = 0;
  for (auto const &r : result) {
    CHECK_TRUE(r.origin < front.size() || r.origin == pareto_polish_created);
    created += r.origin == pareto_polish_created;
    for (auto const &other : result) {
      CHECK_TRUE(&r == &other || !no_worse(other.f, r.f));
    }
    if (r.termination == termination::stationarity_tolerance) {
      CHECK_LE(r.weak_kkt_residual, 100 * params.tolerance);
    }
    bool finite = true;
    for (auto const &v : r.f) {
      finite = finite && isfinite(v);
    }
    auto v = std::begin(r.f);
    for (size_t i = 0; i < m && finite; ++i, ++v) {
      low[i] = (std::min)(low[i], *v);
      high[i] = (std::max)(high[i], *v);
    }
  }
  CHECK_EQUAL(created, result.new_points);
  for (size_t i = 0; i < m; ++i) {
    CHECK_EQUAL(result.ranges[i], high[i] > low[i] ? high[i] - low[i] : Real(1));
  }
  size_t gaps = 0, unresolved = 0, boundary = 0;
  Real worst = 0;
  bool within = true;
  for (auto const &c : result.cells) {
    // Triangles, and the edges where they end, for three objectives; edges otherwise.
    CHECK_TRUE(c.vertices.size() == 2 || (m == 3 && c.vertices.size() == 3));
    for (size_t v : c.vertices) {
      CHECK_LE(v + 1, result.size());
    }
    gaps += c.kind == cell_class::gap;
    unresolved += c.kind == cell_class::unresolved;
    boundary += c.kind == cell_class::boundary;
    CHECK_TRUE(c.kind != cell_class::boundary || m == 3);
    // Unresolved: above the tolerance, or longer than a cell may be (a quarter, at most, for three objectives).
    CHECK_TRUE(c.kind != cell_class::unresolved || !(c.error <= params.refinement_tolerance) ||
               params.max_edge_length < 1 || m >= 3);
    CHECK_TRUE(!c.bound_under_convexity || c.kind == cell_class::convex || c.kind == cell_class::concave);
    // A triangle which is presumed outside the front is as a gap is, and an edge where the triangles end has an
    // estimate, which the achieved tolerance leaves out.
    if (c.kind != cell_class::gap && !(c.kind == cell_class::boundary && c.vertices.size() == 3)) {
      worst = c.kind == cell_class::boundary ? worst : (std::max)(worst, c.error);
      within = within && (c.kind == cell_class::unresolved || c.error <= params.refinement_tolerance);
    }
    CHECK_TRUE(c.multipliers.empty() || c.multipliers.size() == m * c.vertices.size());
  }
  CHECK_EQUAL(gaps, result.gap_cells);
  CHECK_EQUAL(unresolved, result.unresolved_cells);
  CHECK_EQUAL(boundary, result.boundary_cells);
  CHECK_TRUE(result.chart_direction.empty() ? result.cells.empty() || m != 3
                                            : m == 3 && result.chart_direction.size() == 3);
  CHECK_EQUAL(worst, result.achieved_tolerance);
  CHECK_EQUAL(result.status == status::not_refined, !(params.refinement_tolerance > 0));
  CHECK_TRUE(result.status != status::no_cells || result.cells.empty());
  CHECK_TRUE(result.status != status::discovered_cells_resolved || (within && unresolved == 0));
  CHECK_TRUE(result.status != status::cells_unresolved || (within && unresolved > 0));
}

template <class Trace> size_t count_events(Trace const &trace, event_kind k, size_t from = 0) {
  size_t count = 0;
  for (size_t i = from; i < trace.size(); ++i) {
    count += trace[i].kind == k;
  }
  return count;
}

// The trace ends with the cells of the result, once.
template <class Trace, class Result> void check_trace(Trace const &trace, Result const &result) {
  size_t last_round = 0;
  for (size_t i = 0; i < trace.size(); ++i) {
    last_round = trace[i].kind == event_kind::round ? i : last_round;
  }
  CHECK_EQUAL(count_events(trace, event_kind::cell, last_round), result.cells.size());
}

// The triangle between a chord and the tangents at its ends, for a decreasing parabola, convex or concave.
void test_sandwich() {
  using detail::pareto_polish_edge;
  using detail::pareto_polish_sandwich;
  using detail::pareto_polish_simplex_error;
  for (double c : {1.0, -0.3}) {
    for (double x0 : {0.1, 0.4}) {
      const double x1 = x0 + 0.3, xs = c > 0 ? 1 : 0, slope0 = 2 * c * (x0 - xs), slope1 = 2 * c * (x1 - xs);
      const double A[2] = {x0, c * (x0 - xs) * (x0 - xs) + 2}, B[2] = {x1, c * (x1 - xs) * (x1 - xs) + 2};
      // A normal (s, 1 - s) with slope -s / (1 - s); the tangents meet over the middle of the chord.
      const double sa = -slope0 / (1 - slope0), sb = -slope1 / (1 - slope1), xt = (x0 + x1) / 2;
      const double T[2] = {xt, A[1] + slope0 * (xt - x0)};
      const double dx = B[0] - A[0], dy = B[1] - A[1], L = std::sqrt(dx * dx + dy * dy);
      const double height = std::abs((T[0] - A[0]) * dy - (T[1] - A[1]) * dx) / L;
      pareto_polish_edge<double> e;
      const double la[2] = {sa, sa}, lb[2] = {sb, sb};
      pareto_polish_sandwich(A, B, la, lb, 0.0, 0.0, e);
      CHECK_TRUE(e.kind == (c > 0 ? cell_class::convex : cell_class::concave));
      CHECK_TRUE(e.bound);
      CHECK_ABSOLUTE_ERROR(height, e.error, 1e-12);
      CHECK_ABSOLUTE_ERROR(dx * -dy / L, e.box, 1e-12);
      // The parabola is at most half of that from the chord.
      CHECK_LE(std::abs(c) * 0.3 * 0.3 / 4 * dx / L, e.error);
      // The envelope of a cell of two corners is the same triangle, whatever the direction r, and so it is for a
      // front which does not depend on a third objective.
      for (double r0 : {0.5, 0.2, 0.9}) {
        const double r1 = std::sqrt(1 - r0 * r0), r2 = 0.3, s = std::sqrt(1 - r2 * r2);
        const auto two = pareto_polish_simplex_error(V{A[0], A[1], B[0], B[1]}, V{sa, 1 - sa, sb, 1 - sb}, V{}, V{r0, r1}, 2);
        const auto three = pareto_polish_simplex_error(V{A[0], A[1], 0, B[0], B[1], 0.1, A[0], A[1], 1},
                                                       V{sa, 1 - sa, 0, sb, 1 - sb, 0, sa, 1 - sa, 0}, V{},
                                                       V{s * r0, s * r1, r2}, 3);
        CHECK_EQUAL(two.kind, c > 0 ? 0 : 1);
        CHECK_EQUAL(three.kind, c > 0 ? 0 : 1);
        CHECK_ABSOLUTE_ERROR(e.error, two.error, 1e-13);
        CHECK_ABSOLUTE_ERROR(e.error, three.error, 1e-13);
      }
      // Multipliers which are not unique, as at a kink or an end: the tangent which faces the other end is taken.
      const double wide_a[2] = {sa, (1 + sa) / 2}, wide_b[2] = {sb / 2, sb};
      pareto_polish_edge<double> w;
      pareto_polish_sandwich(A, B, wide_a, wide_b, 0.0, 0.0, w);
      CHECK_ABSOLUTE_ERROR(height, w.error, 1e-12);
      // A change of the sign of the curvature: the box, which needs no tangents.
      pareto_polish_edge<double> mixed;
      pareto_polish_sandwich(A, B, la, la, 0.0, 0.0, mixed);
      CHECK_TRUE(mixed.kind == cell_class::mixed);
      CHECK_TRUE(!mixed.bound);
      CHECK_EQUAL(mixed.error, mixed.box);
      // Within the allowances it is not a change of sign.
      pareto_polish_sandwich(A, B, la, la, L, L, mixed);
      CHECK_TRUE(mixed.bound);
    }
  }
  // Ends which are only near the front, with tangents which are only nearly its own: with the distance by which
  // each tangent cuts into the front for its allowance the error is still a bound, and without it is not.
  std::mt19937_64 gen(5);
  size_t needed = 0, bounds = 0;
  for (int trial = 0; trial < 400; ++trial) {
    const double c = trial % 2 ? 1.0 : -0.5, xs = c > 0 ? 1 : 0;
    auto front = [&](double x) { return c * (x - xs) * (x - xs) + 2; };
    const double x0 = 0.05 + 0.5 * unit(gen), x1 = x0 + 0.05 + 0.3 * unit(gen), up = c > 0 ? 1 : -1;
    const double A[2] = {x0, front(x0) + up * 0.004 * unit(gen)}, B[2] = {x1, front(x1) + up * 0.004 * unit(gen)};
    double s[2][2], allowance[2] = {0, 0}, worst = 0;
    const double dx = B[0] - A[0], dy = B[1] - A[1], L = std::sqrt(dx * dx + dy * dy);
    for (int end = 0; end < 2; ++end) {
      const double slope = 2 * c * ((end ? x1 : x0) - xs) * (1 + 0.1 * (unit(gen) - 0.5));
      s[end][0] = s[end][1] = -slope / (1 - slope);
    }
    for (int i = 0; i <= 2000; ++i) {
      const double x = x0 + (x1 - x0) * i / 2000, y = front(x);
      worst = (std::max)(worst, std::abs((x - A[0]) * dy - (y - A[1]) * dx) / L);
      for (int end = 0; end < 2; ++end) {
        const double *E = end ? B : A;
        allowance[end] = (std::max)(allowance[end], -up * (s[end][0] * (x - E[0]) + (1 - s[end][0]) * (y - E[1])));
      }
    }
    pareto_polish_edge<double> with, without;
    pareto_polish_sandwich(A, B, s[0], s[1], allowance[0], allowance[1], with);
    pareto_polish_sandwich(A, B, s[0], s[1], 0.0, 0.0, without);
    if (with.bound) {
      ++bounds;
      CHECK_LE(worst, with.error * (1 + 1e-9) + 1e-12);
    }
    needed += without.bound && without.error < worst * (1 - 1e-6);
  }
  // Measured: 400 and 64 of the 400.
  CHECK_GE(bounds, size_t(200));
  CHECK_GE(needed, size_t(10));
}

// Two objectives: the cells, their classes, and the bound against the front, which is known.
void test_cells() {
  // Where the five pieces of the front of ZDT3 end and begin, in the first objective.
  const double gap[4][2] = {{0.0830015349, 0.1822287280}, {0.2577623634, 0.4093136748}, {0.4538821041, 0.6183967944},
                            {0.6525117038, 0.8233317983}};
  for (int problem : {1, 2, 3}) {
    auto params = unit_box<V>(6, 0.002);
    auto f = [problem](auto const &x) { return zdt(x, problem); };
    std::vector<std::pair<V, std::array<double, 2>>> front;
    // For ZDT3, on each of the five pieces of the front.
    for (double x0 : {0.02, 0.06, 0.2, 0.24, 0.42, 0.44, 0.63, 0.65, 0.83, 0.85}) {
      V x(6, 0.01);
      x[0] = problem == 3 ? x0 : (x0 - 0.02) / 0.83;
      front.emplace_back(x, f(x));
    }
    std::vector<pareto_polish_event<double>> trace;
    const auto result = pareto_polish(f, front, params, {}, &trace);
    check_guarantee(front, result, params);
    check_trace(trace, result);
    CHECK_EQUAL(result.cells.size() + 1, result.size());
    CHECK_GE(result.objective_evaluations, result.new_points);
    size_t of_class = 0;
    if (problem != 3) {
      CHECK_TRUE(result.status == status::discovered_cells_resolved);
      CHECK_EQUAL(result.gap_cells, size_t(0));
    }
    for (auto const &c : result.cells) {
      auto const &a = result[c.vertices[0]].f, &b = result[c.vertices[1]].f;
      CHECK_LE(a[0], b[0]);
      of_class += c.kind == (problem == 2 ? cell_class::concave : cell_class::convex);
      const double dx = (b[0] - a[0]) / result.ranges[0], dy = (b[1] - a[1]) / result.ranges[1];
      const double L = std::sqrt(dx * dx + dy * dy);
      CHECK_LE(c.dominance_bound, L / 2 * 1.001);
      if (c.kind == cell_class::gap) {
        // A presumed gap is one of the four of the front, from end to end.
        bool is_gap = false;
        for (auto const &g : gap) {
          is_gap = is_gap || (std::abs(a[0] - g[0]) < 1e-3 && std::abs(b[0] - g[1]) < 1e-3);
        }
        CHECK_TRUE(is_gap);
      }
      // A cell is a bound only between points with Jacobians of their own; any other has the box.
      if (!c.bound_under_convexity) {
        CHECK_EQUAL(c.error, c.dominance_bound);
        continue;
      }
      CHECK_TRUE(!result[c.vertices[0]].jacobian_from_inside && !result[c.vertices[1]].jacobian_from_inside);
      CHECK_LE(c.error, c.dominance_bound);
      // The largest distance of the front from the chord, in the scale of the result.
      double worst = 0;
      for (int i = 1; i < 64; ++i) {
        V x(6, 0.0);
        x[0] = a[0] + (b[0] - a[0]) * i / 64;
        const auto y = f(x);
        worst = (std::max)(worst, std::abs((y[0] - a[0]) / result.ranges[0] * dy - (y[1] - a[1]) / result.ranges[1] * dx));
      }
      CHECK_LE(worst / L, c.error * 1.01 + 1e-7);
    }
    // Measured: all but the gaps.
    CHECK_GE(2 * of_class, result.cells.size());
    // The gaps of ZDT3 are found (measured: the four of them), and no point lies in one.
    if (problem == 3) {
      CHECK_GE(result.gap_cells, size_t(1));
      CHECK_LE(result.gap_cells, size_t(4));
      for (auto const &r : result) {
        for (auto const &g : gap) {
          CHECK_TRUE(r.f[0] < g[0] + 1e-3 || r.f[0] > g[1] - 1e-3);
        }
      }
    }
    // A guard is not the criterion: the status says which of them ended the refinement.
    params.max_new_points = 4;
    const auto capped = pareto_polish(f, front, params);
    check_guarantee(front, capped, params);
    CHECK_TRUE(capped.status == status::budget_exhausted);
    CHECK_GE(capped.achieved_tolerance, params.refinement_tolerance);
    CHECK_LE(capped.new_points, size_t(4));
    // With one step for each point few searches converge: what holds of every result still does.
    params.max_new_points = 1000;
    params.max_iterations = 1;
    check_guarantee(front, pareto_polish(f, front, params), params);
  }
}

// Two pieces with a jump between them, f = (x, 1 - x) up to x = 0.3 and (x, 1.2 - x) from x = 0.7, and above both
// in between. The searches into the cell between the pieces stop at the jump, where the derivatives say that
// they could go on: their line searches fail, their radius shrinks, and that is not a gap.
struct jump {
  std::array<double, 2> operator()(std::vector<double> const &x) const {
    return {x[0], (x[0] <= 0.3 ? 1.0 : (x[0] < 0.7 ? 2.0 : 1.2)) - x[0]};
  }
};

void test_failed_search_is_no_gap() {
  auto jacobian = [](jump const &f, V const &x, V &values, V &J) {
    const auto y = f(x);
    values.assign(y.begin(), y.end());
    J = {1.0, -1.0};
  };
  const auto params = unit_box<V>(1, 0.01);
  const auto front = front_at(jump{}, {V{0.1}, V{0.3}, V{0.8}, V{0.9}});
  const auto result = pareto_polish(jump{}, front, params, jacobian);
  check_guarantee(front, result, params);
  CHECK_EQUAL(result.gap_cells, size_t(0));
  CHECK_EQUAL(result.unresolved_cells, size_t(1));
  CHECK_TRUE(result.status == status::cells_unresolved);
  for (auto const &c : result.cells) {
    if (c.kind == cell_class::unresolved) {
      CHECK_EQUAL(result[c.vertices[0]].x[0], 0.3);
      CHECK_ABSOLUTE_ERROR(0.7, result[c.vertices[1]].x[0], 1e-6);
    } else {
      CHECK_LE(c.error, params.refinement_tolerance);
    }
  }
}

// f = (x^2, (x - 1)^2) on [0, 1]: every point is Pareto optimal, and the front is convex.
struct parabolas {
  template <class C> std::array<typename C::value_type, 2> operator()(C const &x) const {
    return {x[0] * x[0], (x[0] - 1) * (x[0] - 1)};
  }
};

void test_parameters() {
  auto params = unit_box<V>(1, 1e-100);
  auto front = front_at(parabolas{}, {V{0.2}, V{0.5}, V{0.8}});
  // A tolerance which no number of points reaches ends at the budget.
  params.max_new_points = 4;
  for (double tolerance : {1e-100, std::numeric_limits<double>::denorm_min()}) {
    params.refinement_tolerance = tolerance;
    const auto result = pareto_polish(parabolas{}, front, params);
    check_guarantee(front, result, params);
    CHECK_TRUE(result.status == status::budget_exhausted);
  }
  // The largest length of a cell is the one asked for, and one which no cell could be divided to is refused.
  params.max_new_points = 1000;
  params.refinement_tolerance = 0.01;
  params.max_edge_length = 0.001;
  CHECK_THROW(pareto_polish(parabolas{}, front, params), std::domain_error);
  params.max_edge_length = 0.004;
  const auto result = pareto_polish(parabolas{}, front, params);
  check_guarantee(front, result, params);
  CHECK_TRUE(result.status == status::discovered_cells_resolved);
  for (auto const &c : result.cells) {
    auto const &a = result[c.vertices[0]].f, &b = result[c.vertices[1]].f;
    const double dx = (b[0] - a[0]) / result.ranges[0], dy = (b[1] - a[1]) / result.ranges[1];
    CHECK_LE(std::sqrt(dx * dx + dy * dy), 0.004 * (1 + 1e-12));
  }
  // A front which needs no new point: its cells are made once.
  std::vector<pareto_polish_event<double>> trace;
  auto line = [](auto const &x) {
    using T = typename std::decay<decltype(x[0])>::type;
    return std::array<T, 2>{x[0], 1 - x[0]};
  };
  params = unit_box<V>(1, 0.01);
  front = front_at(line, {V{0.0}, V{0.5}, V{1.0}});
  const auto straight = pareto_polish(line, front, params, {}, &trace);
  check_guarantee(front, straight, params);
  CHECK_EQUAL(straight.size(), size_t(3));
  CHECK_EQUAL(count_events(trace, event_kind::round), size_t(1));
}

// A curved Pareto set (the second objective is a distance in another metric).
// With a loose tolerance of the polish the points are only nearly stationary: the cells are then no bounds, or
// are bounds made wider by the allowance, and a bound still holds against the front, which is known.
void test_curved_set() {
  const double site[2][4] = {{0.2, 0.3, 0.4, 0.5}, {0.8, 0.65, 0.5, 0.35}};
  auto f = [&site](auto const &x) {
    using T = typename std::decay<decltype(x[0])>::type;
    T a = 0 * x[0], b = 0 * x[0];
    for (size_t k = 0; k < x.size(); ++k) {
      a += (x[k] - site[0][k]) * (x[k] - site[0][k]);
      b += (1 + 3 * k) * (x[k] - site[1][k]) * (x[k] - site[1][k]);
    }
    return std::array<T, 2>{a, b};
  };
  // The front: the minimizers of t f_1 + (1 - t) f_2.
  std::vector<std::array<double, 2>> truth;
  for (int i = 0; i <= 4000; ++i) {
    const double t = i / 4000.0;
    V x(4);
    for (size_t k = 0; k < 4; ++k) {
      const double w = (1 - t) * (1 + 3 * k);
      x[k] = (t * site[0][k] + w * site[1][k]) / (t + w);
    }
    truth.push_back(f(x));
  }
  auto params = unit_box<V>(4, 0.001);
  const auto front = scattered_front<double>(f, 4, 6, 3);
  for (double tolerance : {1.5e-8, 1e-3}) {
    params.tolerance = tolerance;
    const auto result = pareto_polish(f, front, params);
    check_guarantee(front, result, params);
    CHECK_TRUE(result.status == status::discovered_cells_resolved);
    size_t estimates = 0;
    for (auto const &c : result.cells) {
      auto const &a = result[c.vertices[0]].f, &b = result[c.vertices[1]].f;
      const double dx = (b[0] - a[0]) / result.ranges[0], dy = (b[1] - a[1]) / result.ranges[1];
      double worst = 0;
      for (auto const &y : truth) {
        if (y[0] >= a[0] && y[0] <= b[0]) {
          worst = (std::max)(worst, std::abs((y[0] - a[0]) / result.ranges[0] * dy - (y[1] - a[1]) / result.ranges[1] * dx));
        }
      }
      estimates += !c.bound_under_convexity;
      if (c.bound_under_convexity) {
        // Two points of the front may reach the same end to within rounding: the cell between them has no class.
        CHECK_TRUE(c.kind == cell_class::convex || dx * dx + dy * dy < 1e-20);
        CHECK_LE(worst / std::sqrt(dx * dx + dy * dy), c.error * 1.01 + 1e-7);
      }
    }
    if (tolerance < 1e-6) {
      // Measured: every cell is a bound.
      CHECK_LE(2 * estimates, result.cells.size());
    }
  }
  // Points which stop at the cap on their steps far from the front: a cell is a bound only between points which
  // are stationary to the hundred tolerances of the documentation and its allowance for rounding.
  params.tolerance = 1.5e-8;
  params.max_iterations = 2;
  params.max_new_points = 20;
  const auto capped = pareto_polish(f, front, params);
  check_guarantee(front, capped, params);
  size_t at_cap = 0;
  for (auto const &r : capped) {
    at_cap += r.termination == termination::iteration_cap;
  }
  CHECK_GE(at_cap, size_t(1));
  for (auto const &c : capped.cells) {
    if (c.bound_under_convexity) {
      CHECK_LE((std::max)(capped[c.vertices[0]].weak_kkt_residual, capped[c.vertices[1]].weak_kkt_residual), 1e-4);
    }
  }
}

// On ZDT2 a point with x[0] = 0 is weakly stationary whatever the others are, and the derivative of f2 in x[0]
// vanishes there; on ZDT1 the Jacobian is not finite there.
void test_bounds() {
  auto params = unit_box<V>(6, 0.0);
  params.upper_bounds[5] = 0.0; // a variable with no room
  for (int problem : {1, 2}) {
    auto f = [problem](auto const &x) { return zdt(x, problem); };
    std::vector<std::pair<V, std::array<double, 2>>> front;
    for (double x0 : {0.0, 0.0, problem == 1 ? 0.3 : 0.0, problem == 1 ? 1.0 : 0.0}) {
      V x(6, 0.25);
      x[0] = x0;
      x[5] = 0;
      front.emplace_back(x, f(x));
    }
    params.refinement_tolerance = 0;
    auto result = pareto_polish(f, front, params);
    check_guarantee(front, result, params);
    // Points which start alike end alike, and one of them is returned.
    CHECK_LE(result.size(), size_t(problem == 1 ? 3 : 1));
    params.refinement_tolerance = 0.01;
    result = pareto_polish(f, front, params);
    check_guarantee(front, result, params);
    double high = 0, widest = 0;
    std::vector<double> f1;
    for (auto const &c : result.cells) {
      CHECK_TRUE(!(c.bound_under_convexity &&
                   (result[c.vertices[0]].jacobian_from_inside || result[c.vertices[1]].jacobian_from_inside)));
    }
    for (auto const &r : result) {
      // The Jacobian of ZDT1 is not finite where x[0] = 0: that of a point just inside is used, and says so.
      CHECK_EQUAL(r.jacobian_from_inside, problem == 1 && r.x[0] == 0);
      // Measured: zero.
      CHECK_LE(r.x[1], 1e-3);
      CHECK_EQUAL(r.x[5], 0.0);
      f1.push_back(r.f[0]);
      high = (std::max)(high, r.f[0]);
    }
    std::sort(f1.begin(), f1.end());
    for (size_t i = 1; i < f1.size(); ++i) {
      widest = (std::max)(widest, f1[i] - f1[i - 1]);
    }
    // The whole front is found from one end of it. Measured: to 1, with no two points farther apart than 0.25.
    CHECK_GE(high, 0.9);
    CHECK_LE(widest, 0.4);
  }
}

void test_edge_cases() {
  using V = std::array<double, 3>;
  using Front = std::vector<std::pair<V, std::array<double, 2>>>;
  // Not finite for x[2] < 0.4:
  auto f = [](auto const &x) {
    using T = typename std::decay<decltype(x[0])>::type;
    using std::sqrt;
    T a = x[0] * x[0] + x[1] * x[1] + sqrt(x[2] - 0.4);
    T b = (x[0] - 1) * (x[0] - 1) + x[1] * x[1] + sqrt(x[2] - 0.4);
    return std::array<T, 2>{a, b};
  };
  pareto_polish_parameters<V> params;
  params.lower_bounds = {0, 0, 0};
  params.upper_bounds = {1, 1, 1};
  CHECK_EQUAL(pareto_polish(f, Front{}, params).size(), size_t(0));
  CHECK_TRUE(pareto_polish(f, Front{}, params).status == status::no_cells);
  params.refinement_tolerance = 0;
  CHECK_TRUE(pareto_polish(f, Front{}, params).status == status::not_refined);
  params.refinement_tolerance = 0.01;
  Front front = front_at(f, {V{0.5, 0.5, 0.9}, V{0.5, 0.5, 0.9}, V{0.2, 0.0, 0.4}, V{0.3, 0.3, 0.1}, V{1, 1, 1}});
  // A point whose objectives are not finite is returned as it is, once, whatever the others do.
  const auto result = pareto_polish(f, front, params);
  check_guarantee(front, result, params);
  size_t unchanged = 0;
  for (auto const &r : result) {
    unchanged += r.origin == 3 && r.x == front[3].first && std::isnan(r.f[0]) && std::isnan(r.f[1]) &&
                 r.termination == termination::non_finite && r.objective_evaluations == 0;
  }
  CHECK_EQUAL(unchanged, size_t(1));
  check_guarantee(Front{front[0]}, pareto_polish(f, Front{front[0]}, params), params);
  // Points with no finite objectives make no cell.
  const Front bad = front_at(f, {front[3].first, V{0.1, 0.1, 0.2}});
  const auto none = pareto_polish(f, bad, params);
  check_guarantee(bad, none, params);
  CHECK_EQUAL(none.size(), size_t(2));
  CHECK_EQUAL(none.cells.size(), size_t(0));
  CHECK_TRUE(none.status == status::no_cells);
  CHECK_EQUAL(none.polished_inputs, size_t(0));
  // A box which is a point.
  auto fixed = params;
  fixed.lower_bounds = fixed.upper_bounds = front[0].first;
  const auto single = pareto_polish(f, Front{front[0]}, fixed);
  CHECK_EQUAL(single.size(), size_t(1));
  CHECK_TRUE(single.status == status::no_cells);
  // No Jacobian anywhere: the cells are those of the box, and nothing is a bound. It is asked for at a point, and
  // once more just inside the bounds for a point on them, and not again.
  auto no_jacobian = [](decltype(f) const &g, std::vector<double> const &x, std::vector<double> &values, std::vector<double> &J) {
    const auto y = g(x);
    values.assign(y.begin(), y.end());
    J.assign(6, std::numeric_limits<double>::quiet_NaN());
  };
  const Front two = front_at(f, {V{0.2, 0.0, 0.5}, V{0.8, 0.5, 0.5}});
  fixed = params;
  fixed.max_new_points = 12;
  const auto blind = pareto_polish(f, two, fixed, no_jacobian);
  check_guarantee(two, blind, fixed);
  CHECK_GE(blind.cells.size(), size_t(1));
  CHECK_TRUE(blind.status == status::budget_exhausted || blind.status == status::discovered_cells_resolved);
  for (auto const &c : blind.cells) {
    CHECK_TRUE(!c.bound_under_convexity && c.kind == cell_class::mixed);
  }
  for (auto const &r : blind) {
    CHECK_TRUE(r.termination == termination::non_finite);
    CHECK_EQUAL(r.jacobian_evaluations, size_t(r.x[1] == 0 ? 2 : 1));
  }
  // What the objectives throw arrives here, from whichever thread.
  auto thrower = [](auto const &x) {
    using T = typename std::decay<decltype(x[0])>::type;
    if (x.size() == 3) {
      throw std::runtime_error("objective");
    }
    return std::array<T, 2>{x[0], x[1]};
  };
  for (unsigned threads : {1u, several_threads}) {
    params.threads = threads;
    CHECK_THROW(pareto_polish(thrower, front, params), std::runtime_error);
  }
  params.tolerance = 0;
  CHECK_THROW(pareto_polish(f, front, params), std::domain_error);
  params.tolerance = 1e-8;
  params.threads = 0;
  CHECK_THROW(pareto_polish(f, front, params), std::invalid_argument);
  params.threads = 1;
  front[0].first[0] = 2;
  CHECK_THROW(pareto_polish(f, front, params), std::domain_error);
}

// The middle of the first cell is not finite: the new point starts there once, and is not asked for again.
void test_start_not_finite() {
  struct hole {
    std::vector<double> *asked;
    std::array<double, 2> operator()(V const &x) const {
      if (x[1] == 0.5) {
        asked->push_back(x[0]);
        return {std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN()};
      }
      return {x[0] * x[0], (x[0] - 1) * (x[0] - 1)};
    }
  };
  auto jacobian = [](hole const &, V const &x, V &values, V &J) {
    values = {x[0] * x[0], (x[0] - 1) * (x[0] - 1)};
    J = {2 * x[0], 0.0, 2 * (x[0] - 1), 0.0};
  };
  std::vector<double> asked;
  const hole f{&asked};
  const auto params = unit_box<V>(2, 0.01);
  const auto front = front_at(f, {V{0.25, 0.25}, V{0.75, 0.75}});
  const auto result = pareto_polish(f, front, params, jacobian);
  check_guarantee(front, result, params);
  CHECK_GE(asked.size(), size_t(1));
  std::sort(asked.begin(), asked.end());
  CHECK_TRUE(std::adjacent_find(asked.begin(), asked.end()) == asked.end());
}

// Two branches of stationary points, one of which dominates the other: what a front which is not mature, or a
// Pareto set of several components, looks like from inside.
void test_dominated_inputs() {
  auto f = [](auto const &x) {
    using T = typename std::decay<decltype(x[0])>::type;
    T well = (x[1] * x[1] - 1) * (x[1] * x[1] - 1) + (x[1] + 1) / 4; // minima near -1 (0) and 1 (1/2)
    return std::array<T, 2>{x[0] * x[0], (x[0] - 1) * (x[0] - 1) + well};
  };
  pareto_polish_parameters<V> params;
  params.lower_bounds = {0.0, -1.5};
  params.upper_bounds = {1.0, 1.5};
  params.refinement_tolerance = 0.01;
  auto front = front_at(f, {V{0.1, -0.9}, V{0.5, -1.1}, V{0.9, -1.0}, V{0.3, 0.9}, V{0.7, 1.1}});
  const auto result = pareto_polish(f, front, params);
  check_guarantee(front, result, params);
  // Measured: all five reach the tolerance.
  CHECK_GE(result.polished_inputs, size_t(3));
  CHECK_LE(result.dominated_inputs, (std::min)(result.polished_inputs, size_t(2)));
  for (auto const &r : result) {
    CHECK_LE(r.x[1], -0.5);
  }
  CHECK_TRUE(result.status == status::discovered_cells_resolved);
  // One point, twice: its copy is left out, not counted as dominated, and the front is found from the other.
  front.assign(2, front[1]);
  const auto twice = pareto_polish(f, front, params);
  check_guarantee(front, twice, params);
  CHECK_TRUE(twice.status == status::discovered_cells_resolved);
  CHECK_LE(twice.size(), twice.new_points + 1);
  CHECK_GE(twice.size(), size_t(3));
  CHECK_TRUE(twice.polished_inputs == 0 || twice.polished_inputs == 2);
  CHECK_TRUE(twice.dominated_inputs == 0 || twice.size() == twice.new_points);
}

// Two sites, the second in a diagonal metric with weights from 1 to 30: strictly convex quadratics with a curved
// Pareto set, on which the first order steps are slow. With `warp` the same in the squares of the variables,
// which is not convex.
struct bent {
  bool warp = false;
  template <class C> std::array<typename C::value_type, 2> operator()(C const &x) const {
    using T = typename C::value_type;
    T a = 0 * x[0], b = 0 * x[0];
    for (size_t k = 0; k < x.size(); ++k) {
      const double s0 = 0.3 + 0.05 * static_cast<double>(k % 5), s1 = 0.75 - 0.07 * static_cast<double>(k % 4);
      const T y = warp ? T(x[k] * x[k]) : T(x[k]);
      a += (y - s0) * (y - s0);
      b += weight(k) * (y - s1) * (y - s1);
    }
    return {a, b};
  }
  static double weight(size_t k) { return static_cast<double>(1 + 4 * k); }
};

// The same objectives for vectors of double alone, with the derivatives which the caller then supplies.
struct bent_double {
  std::array<double, 2> operator()(std::vector<double> const &x) const { return bent{}(x); }
};
struct bent_jacobian {
  void operator()(bent_double const &f, std::vector<double> const &x, std::vector<double> &values,
                  std::vector<double> &J) const {
    const auto y = f(x);
    const size_t n = x.size();
    values.assign(y.begin(), y.end());
    J.resize(2 * n);
    for (size_t k = 0; k < n; ++k) {
      J[k] = 2 * (x[k] - (0.3 + 0.05 * static_cast<double>(k % 5)));
      J[n + k] = 2 * bent::weight(k) * (x[k] - (0.75 - 0.07 * static_cast<double>(k % 4)));
    }
  }
};
struct bent_hvp {
  void operator()(bent_double const &, std::vector<double> const &x, std::vector<double> const &weights,
                  std::vector<double> const &directions, std::vector<double> &products) const {
    const size_t n = x.size();
    products.resize(directions.size());
    for (size_t j = 0; j < directions.size(); ++j) {
      products[j] = 2 * (weights[0] + weights[1] * bent::weight(j % n)) * directions[j];
    }
  }
};

// The corrector of second order: the same front, with the same guarantee. What it saves is for a benchmark to
// measure, and is not tested.
void test_corrector() {
  const size_t n = 8;
  auto params = unit_box<V>(n, 0.0005);
  const auto front = scattered_front<double>(bent{}, n, 24, 11);
  const pareto_polish_reverse_mode_policy<double> jacobian;
  std::vector<pareto_polish_event<double>> trace;
  const auto first = pareto_polish(bent{}, front, params, jacobian);
  check_guarantee(front, first, params);
  CHECK_TRUE(first.status == status::discovered_cells_resolved);
  CHECK_EQUAL(first.hvp_evaluations, size_t(0));
  // Without products nothing changes, to the last bit.
  check_same(first, pareto_polish(bent{}, front, params, jacobian, nullptr, pareto_polish_no_hvp()));
  // The multipliers of a cell are the tangents it was compared with. The front is convex, so every point is on or
  // above each of them, and they meet the chord in the triangle whose height is the error of the cell, but for
  // the allowances. Measured: no point is below a tangent, and the errors agree to 2e-8.
  for (auto const &c : first.cells) {
    CHECK_TRUE(c.bound_under_convexity && c.multipliers.size() == 4);
    auto const &a = first[c.vertices[0]], &b = first[c.vertices[1]];
    double const *la = &c.multipliers[0], *lb = &c.multipliers[2];
    auto above = [&](double const *l, std::array<double, 2> const &at, std::array<double, 2> const &y) {
      const double h = l[0] * (y[0] - at[0]) / first.ranges[0] + l[1] * (y[1] - at[1]) / first.ranges[1];
      return h / std::hypot(l[0], l[1]);
    };
    for (auto const &r : first) {
      CHECK_GE(above(la, a.f, r.f), -1e-8);
      CHECK_GE(above(lb, b.f, r.f), -1e-8);
    }
    const double L = std::hypot((b.f[0] - a.f[0]) / first.ranges[0], (b.f[1] - a.f[1]) / first.ranges[1]);
    const double sine = std::abs(la[0] * lb[1] - la[1] * lb[0]) / (std::hypot(la[0], la[1]) * std::hypot(lb[0], lb[1]));
    CHECK_ABSOLUTE_ERROR(above(la, a.f, b.f) * above(lb, b.f, a.f) / (L * sine), c.error, 1e-7);
  }
  for (auto const &r : first) {
    CHECK_ABSOLUTE_ERROR(1.0, r.multipliers[0] + r.multipliers[1], 1e-14);
  }
  for (int variant = 0; variant < 3; ++variant) {
    trace.clear();
    const auto second = variant == 0   ? pareto_polish(bent{}, front, params, jacobian, &trace,
                                                       pareto_polish_reverse_mode_hvp_policy<double>())
                        : variant == 1 ? pareto_polish(bent{}, front, params, jacobian, &trace,
                                                       pareto_polish_jacobian_difference_hvp())
                                       : pareto_polish(bent_double{}, front, params, bent_jacobian{}, &trace, bent_hvp{});
    check_guarantee(front, second, params);
    CHECK_TRUE(second.status == status::discovered_cells_resolved);
    // Measured: 168 steps of the corrector, and every point stationary to 4e-9.
    CHECK_GE(count_events(trace, event_kind::corrector_step), size_t(1));
    size_t products = 0, steps = 0, stationary = 0, traced = 0;
    for (auto const &r : second) {
      products += r.hvp_evaluations;
      steps += r.corrector_steps;
      stationary += r.weak_kkt_residual <= 1e-5;
    }
    for (auto const &e : trace) {
      traced = (std::max)(traced, e.hvp_evaluations);
    }
    CHECK_GE(2 * stationary, second.size());
    CHECK_LE(steps, count_events(trace, event_kind::corrector_step));
    CHECK_LE(products, second.hvp_evaluations);
    CHECK_EQUAL(second.hvp_evaluations > 0, variant != 1);
    CHECK_EQUAL(traced > 0, variant != 1);
    CHECK_LE(traced, second.hvp_evaluations);
    if (variant == 0) {
      // On several threads, the same.
      auto threaded = params;
      threaded.threads = several_threads;
      check_same(second, pareto_polish(bent{}, front, threaded, jacobian, nullptr, pareto_polish_reverse_mode_hvp_policy<double>()));
    }
  }
}

// Objectives which are not convex: the Hessian of the corrector is indefinite away from the front.
void test_corrector_indefinite() {
  const size_t n = 6;
  const auto params = unit_box<V>(n, 0.002);
  bent f;
  f.warp = true;
  const auto front = scattered_front<double>(f, n, 20, 5);
  std::vector<pareto_polish_event<double>> trace;
  const auto second = pareto_polish(f, front, params, pareto_polish_reverse_mode_policy<double>(), &trace,
                                    pareto_polish_reverse_mode_hvp_policy<double>());
  check_guarantee(front, second, params);
  // Measured: 263 steps of the corrector, and 25 iterations which it left to the first order step.
  CHECK_GE(count_events(trace, event_kind::corrector_step), size_t(1));
  for (auto const &e : trace) {
    CHECK_TRUE(e.kind != event_kind::corrector_fallback || e.first <= 3);
  }
  size_t fallbacks = 0;
  for (auto const &r : second) {
    fallbacks += r.corrector_fallbacks;
    for (double v : r.x) {
      CHECK_TRUE(v >= 0 && v <= 1);
    }
  }
  CHECK_LE(fallbacks, count_events(trace, event_kind::corrector_fallback));
}

// The calling thread is one of the workers, and with one thread the only one. The built-in policies record on
// its reverse mode tapes after what is there, and take off them what they recorded.
void test_threads_and_tapes() {
  namespace rdiff = boost::math::differentiation::reverse_mode;
  auto &tape1 = rdiff::get_active_tape<double, 1>();
  auto &tape2 = rdiff::get_active_tape<double, 2>();
  rdiff::rvar<double, 1> a = 2.0;
  rdiff::rvar<double, 1> z = a * a * a;
  rdiff::rvar<double, 2> b = 3.0;
  rdiff::rvar<double, 2> w = b * b;
  // The derivative of w is itself recorded, on the tape of the first order: it must outlive the policies.
  tape2.zero_grad();
  w.backward();
  tape1.add_checkpoint();
  const auto size1 = tape1.end() - tape1.begin(), size2 = tape2.end() - tape2.begin();
  std::mutex mutex;
  std::set<std::thread::id> threads;
  auto f = [&mutex, &threads](auto const &x) {
    {
      std::lock_guard<std::mutex> lock(mutex);
      threads.insert(std::this_thread::get_id());
    }
    return bent{}(x);
  };
  const size_t n = 4;
  auto params = unit_box<V>(n, 0.01);
  const auto front = scattered_front<double>(bent{}, n, 12, 2);
  auto caller_recordings = [&]() {
    CHECK_EQUAL(tape1.end() - tape1.begin(), size1);
    CHECK_EQUAL(tape2.end() - tape2.begin(), size2);
    tape1.zero_grad();
    z.backward();
    CHECK_EQUAL(a.adjoint(), 12.0);
    CHECK_EQUAL(w.item(), 9.0);
    CHECK_EQUAL(b.adjoint().item(), 6.0);
    tape1.zero_grad();
    b.adjoint().backward();
    CHECK_EQUAL(b.get_value().adjoint(), 2.0);
    // The last checkpoint is the caller's again.
    rdiff::rvar<double, 1> more = a * z;
    CHECK_EQUAL(more.item(), 16.0);
    tape1.rewind_to_last_checkpoint();
    CHECK_EQUAL(tape1.end() - tape1.begin(), size1);
  };
  const auto one = pareto_polish(f, front, params, pareto_polish_reverse_mode_policy<double>(), nullptr,
                                   pareto_polish_reverse_mode_hvp_policy<double>());
  check_guarantee(front, one, params);
  CHECK_EQUAL(threads.size(), size_t(1));
  CHECK_EQUAL(threads.count(std::this_thread::get_id()), size_t(1));
  caller_recordings();
  params.threads = several_threads;
  check_same(one, pareto_polish(f, front, params, pareto_polish_reverse_mode_policy<double>(), nullptr,
                                pareto_polish_reverse_mode_hvp_policy<double>()));
  caller_recordings();
  // Policies which the caller has used: their copies record arguments of their own, on the tapes of their threads.
  pareto_polish_reverse_mode_policy<double> jacobian;
  pareto_polish_reverse_mode_hvp_policy<double> hvp;
  V values, J, again, products;
  const V x = front[0].first, weights{1.0, 2.0}, direction{1.0, 0.0, 0.0, 0.0};
  jacobian(f, x, values, J);
  hvp(f, x, weights, direction, products);
  CHECK_EQUAL(products.size(), n);
  CHECK_ABSOLUTE_ERROR(2 * (weights[0] + weights[1] * bent::weight(0)), products[0], 1e-12);
  for (unsigned count : {1u, several_threads}) {
    params.threads = count;
    check_same(one, pareto_polish(f, front, params, jacobian, nullptr, hvp));
  }
  jacobian(f, x, values, again);
  CHECK_TRUE(J == again);
  hvp.release();
  jacobian.release();
  caller_recordings();
  // The policies again with arguments of other dimensions: each recording replaces the one before it, and
  // release() leaves the tapes as they were. With no arguments, or directions which are not whole, no products.
  auto g = [](auto const &y) {
    using T = typename std::decay<decltype(y[0])>::type;
    T sum = 0 * y[0];
    for (auto const &v : y) {
      sum += v * v * v;
    }
    return std::vector<T>{sum, 2 * sum};
  };
  for (size_t dimension : {size_t(1), size_t(3), size_t(2), size_t(3)}) {
    const V y(dimension, 0.5), unit = {1.0, 1.0};
    V directions(dimension * dimension, 0.0);
    for (size_t k = 0; k < dimension; ++k) {
      directions[k * dimension + k] = 1;
    }
    jacobian(g, y, values, again);
    CHECK_EQUAL(again.size(), 2 * dimension);
    CHECK_EQUAL(again[0], 0.75);
    hvp(g, y, unit, directions, products);
    CHECK_EQUAL(products.size(), dimension * dimension);
    CHECK_ABSOLUTE_ERROR(9.0, products[0], 1e-12);
    hvp(g, y, unit, V(dimension + 1, 1.0), products);
    CHECK_EQUAL(products.size(), dimension == 1 ? size_t(2) : size_t(0));
  }
  hvp(g, V(), V{1.0, 1.0}, V(), products);
  CHECK_TRUE(products.empty());
  hvp.release();
  jacobian.release();
  caller_recordings();
  tape1.clear();
  tape2.clear();
}

// The curvature diagnostic: f_2 - f_1 depends on x[0] alone, so the directions which change no objective to first
// order are those of x[1], in which both objectives have a minimum at 1 and at -1 and a maximum at 0.
void test_curvature() {
  auto f = [](auto const &x) {
    using T = typename std::decay<decltype(x[0])>::type;
    T well = (x[1] * x[1] - 1) * (x[1] * x[1] - 1);
    return std::array<T, 2>{x[0] * x[0] + well, (x[0] - 1) * (x[0] - 1) + well};
  };
  pareto_polish_parameters<V> params;
  params.lower_bounds = {0.0, -1.5};
  params.upper_bounds = {1.0, 1.5};
  params.refinement_tolerance = 0;
  params.curvature_diagnostic = true;
  const pareto_polish_reverse_mode_policy<double> jacobian;
  for (int branch = 0; branch < 2; ++branch) {
    // On the saddles, which stay where they are as the derivative in x[1] vanishes there, and near the minima.
    const auto front = front_at(f, {V{0.3, branch ? 1.02 : 0.0}, V{0.6, branch ? -0.97 : 0.0}, V{0.0, branch ? 1.0 : 0.0}});
    for (int variant = 0; variant < 2; ++variant) {
      const auto result = variant ? pareto_polish(f, front, params, jacobian, nullptr, pareto_polish_jacobian_difference_hvp())
                                  : pareto_polish(f, front, params, jacobian, nullptr, pareto_polish_reverse_mode_hvp_policy<double>());
      check_guarantee(front, result, params);
      CHECK_EQUAL(result.size(), size_t(3));
      size_t computed = 0;
      for (auto const &r : result) {
        // Measured: to 1e-9.
        CHECK_ABSOLUTE_ERROR(branch ? 1.0 : 0.0, std::abs(r.x[1]), 1e-4);
        if (r.origin == 2 || r.termination != termination::stationarity_tolerance) {
          // Where the point is not stationary the diagnostic says nothing; and on the bound of x[0] the gradients
          // of both objectives vanish in the other variable, so that the multipliers are not unique.
          CHECK_NAN(r.least_reduced_curvature);
          continue;
        }
        // The second derivative of the well at the point, -4 or near 8, times the sum of the multipliers over the
        // ranges. The multipliers of the objectives over their ranges are as (1 - x[0]) / r_2 to x[0] / r_1.
        const double a = r.x[0] / result.ranges[0], b = (1 - r.x[0]) / result.ranges[1];
        const double expected = (12 * r.x[1] * r.x[1] - 4) * (b / result.ranges[0] + a / result.ranges[1]) / (a + b);
        CHECK_ABSOLUTE_ERROR(expected, r.least_reduced_curvature, 1e-4 * std::abs(expected));
        CHECK_TRUE(variant ? r.hvp_evaluations == 0 : r.hvp_evaluations >= 2);
        ++computed;
      }
      // Measured: both. The saddles do not move, so they are stationary from the start.
      CHECK_GE(computed, size_t(branch ? 0 : 2));
    }
    // Not asked for, or with no products to compute it from.
    for (auto const &r : pareto_polish(f, front, params)) {
      CHECK_NAN(r.least_reduced_curvature);
    }
    auto silent = params;
    silent.curvature_diagnostic = false;
    for (auto const &r : pareto_polish(f, front, silent, jacobian, nullptr, pareto_polish_jacobian_difference_hvp())) {
      CHECK_NAN(r.least_reduced_curvature);
    }
    // What the policy for the products throws arrives here, from whichever thread.
    auto thrower = [](decltype(f) const &, V const &, V const &, V const &, V &) { throw std::runtime_error("products"); };
    auto threaded = params;
    threaded.threads = several_threads;
    CHECK_THROW(pareto_polish(f, front, threaded, jacobian, nullptr, thrower), std::runtime_error);
  }
  // The least eigenvalue of a matrix on the complement of the rows of another.
  using detail::pareto_polish_reduced_curvature;
  const V W{1, 0, 0, 0, -2, 0, 0, 0, 3};
  size_t rank = 7;
  CHECK_ABSOLUTE_ERROR(-2.0, pareto_polish_reduced_curvature(W, V{}, 3, 0.0, rank), 1e-12);
  CHECK_EQUAL(rank, size_t(0));
  CHECK_ABSOLUTE_ERROR(1.0, pareto_polish_reduced_curvature(W, V{0, 5, 0, 0, 5, 0}, 3, 1e-8, rank), 1e-12);
  CHECK_EQUAL(rank, size_t(1));
  // On the complement of (1, 1, 0): the least eigenvalue of [[-1/2, 0], [0, 3]].
  CHECK_ABSOLUTE_ERROR(-0.5, pareto_polish_reduced_curvature(W, V{1, 1, 0}, 3, 1e-8, rank), 1e-12);
  CHECK_TRUE(std::isinf(pareto_polish_reduced_curvature(W, V{1, 0, 0, 0, 1, 0, 0, 0, 1}, 3, 1e-8, rank)));
  CHECK_EQUAL(rank, size_t(3));
  // Entries whose squares overflow: the eigenvalues of [[1, 2], [2, 1]] are -1 and 3.
  CHECK_RELATIVE_ERROR(-1e200, pareto_polish_reduced_curvature(V{1e200, 2e200, 2e200, 1e200}, V{}, 2, 0.0, rank), 1e-12);
  V big = W, row{1e200, 1e200, 0};
  for (auto &v : big) {
    v *= 1e200;
  }
  CHECK_RELATIVE_ERROR(-0.5e200, pareto_polish_reduced_curvature(big, row, 3, 1e-8, rank), 1e-12);
  CHECK_EQUAL(rank, size_t(1));
  big[4] = std::numeric_limits<double>::infinity();
  CHECK_NAN(pareto_polish_reduced_curvature(big, row, 3, 1e-8, rank));
}

// The singular value decomposition: orthogonality, reconstruction and singular values which are known.
void test_svd() {
  using detail::pareto_polish_svd;
  V U, s, W;
  auto check = [&](V const &A, size_t m, size_t n, double tolerance) {
    pareto_polish_svd(A, m, n, U, s, W);
    double largest = 0;
    for (double a : A) {
      largest = (std::max)(largest, std::abs(a));
    }
    for (size_t i = 0; i < m; ++i) {
      CHECK_TRUE(s[i] >= 0 && (i == 0 || s[i] <= s[i - 1]));
      for (size_t j = 0; j < m; ++j) {
        double uu = 0, vv = 0;
        for (size_t k = 0; k < m; ++k) {
          uu += U[i * m + k] * U[j * m + k];
        }
        for (size_t k = 0; k < n; ++k) {
          vv += W[i * n + k] * W[j * n + k];
        }
        CHECK_ABSOLUTE_ERROR(i == j ? 1.0 : 0.0, uu, 1e-13);
        // Rows of V of singular values at the level of rounding are directions of rounding errors.
        if (s[i] > 1e-12 * s[0] && s[j] > 1e-12 * s[0]) {
          CHECK_ABSOLUTE_ERROR(i == j ? 1.0 : 0.0, vv, 1e-10);
        }
      }
      for (size_t j = 0; j < n; ++j) {
        double a = 0;
        for (size_t k = 0; k < m; ++k) {
          a += U[k * m + i] * s[k] * W[k * n + j];
        }
        CHECK_ABSOLUTE_ERROR(A[i * n + j], a, tolerance * largest + 1e-300);
      }
    }
  };
  std::mt19937_64 gen(11);
  // Three orthonormal rows of five columns, and a rotation.
  V Q = {1, 2, 2, 0, 0, 2, 1, -2, 0, 0, 2, -2, 1, 0, 0}, R(9);
  for (auto &q : Q) {
    q /= 3;
  }
  // Rows of lengths 1, 1e-8 and 1e-15: every singular value to full relative accuracy.
  const double graded[3] = {1.0, 1e-8, 1e-15};
  V A(15);
  for (size_t i = 0; i < 3; ++i) {
    for (size_t k = 0; k < 5; ++k) {
      A[i * 5 + k] = graded[i] * Q[i * 5 + k];
    }
  }
  check(A, 3, 5, 1e-14);
  for (size_t i = 0; i < 3; ++i) {
    CHECK_RELATIVE_ERROR(graded[i], s[i], 1e-12);
  }
  // The same rows mixed by a rotation: the second is then known to 1e-16 / 1e-8 of itself, and no better.
  const double c = std::cos(0.7), sn = std::sin(0.7);
  R = {c, sn, 0, -sn, c, 0, 0, 0, 1};
  V B(15, 0.0);
  for (size_t i = 0; i < 3; ++i) {
    for (size_t j = 0; j < 3; ++j) {
      for (size_t k = 0; k < 5; ++k) {
        B[i * 5 + k] += R[i * 3 + j] * A[j * 5 + k];
      }
    }
  }
  check(B, 3, 5, 1e-14);
  CHECK_RELATIVE_ERROR(1.0, s[0], 1e-14);
  CHECK_RELATIVE_ERROR(1e-8, s[1], 1e-6);
  // Dense, and with more rows than columns.
  V D(24), E(15);
  for (auto &d : D) {
    d = unit(gen) - 0.5;
  }
  for (auto &e : E) {
    e = unit(gen) - 0.5;
  }
  check(D, 4, 6, 1e-14);
  CHECK_GE(s[3], 1e-3);
  check(E, 5, 3, 1e-14);
  CHECK_GE(s[2], 1e-3);
  CHECK_LE(s[3], 1e-15);
  check({3, 0, 0, 4}, 2, 2, 1e-15);
  CHECK_RELATIVE_ERROR(4.0, s[0], 1e-15);
  CHECK_RELATIVE_ERROR(3.0, s[1], 1e-15);
  // Rank one, and rank two of three: the left singular vector of zero is the combination of the rows which
  // vanishes.
  check({1, 1, 1, 1}, 2, 2, 1e-15);
  CHECK_RELATIVE_ERROR(2.0, s[0], 1e-15);
  CHECK_LE(s[1], 1e-15);
  check({1, 2, 3, 4, 2, 4, 6, 8, 1, 0, -1, 0}, 3, 4, 1e-14);
  CHECK_LE(s[2], 1e-14);
  CHECK_ABSOLUTE_ERROR(1.0, std::abs(2 * U[6] - U[7]) / std::sqrt(5.0), 1e-13);
  CHECK_ABSOLUTE_ERROR(0.0, U[8], 1e-13);
  check({3, -4, 12}, 1, 3, 1e-15);
  CHECK_RELATIVE_ERROR(13.0, s[0], 1e-15);
  CHECK_ABSOLUTE_ERROR(-4.0 / 13, W[1], 1e-15);
  check(V(6, 0.0), 2, 3, 1e-15);
  CHECK_EQUAL(s[0], 0.0);
  // In single precision.
  std::vector<float> Af(A.begin(), A.end()), Uf, sf, Wf;
  pareto_polish_svd(Af, 3, 5, Uf, sf, Wf);
  CHECK_RELATIVE_ERROR(1e-8f, sf[1], 1e-5f);
  CHECK_RELATIVE_ERROR(1e-15f, sf[2], 1e-5f);
}

// The triangulation: it is one (the count of Euler, every triangle counter-clockwise, every edge at most twice),
// and no point is inside the circle of a triangle.
template <class Real> size_t check_delaunay(std::vector<std::array<Real, 2>> const &points, size_t expected_left_out) {
  detail::pareto_polish_delaunay<Real> mesh;
  for (auto const &p : points) {
    mesh.insert(p[0], p[1]);
  }
  size_t in = 0;
  for (size_t p = 0; p < points.size(); ++p) {
    in += mesh.inserted(p);
  }
  CHECK_EQUAL(points.size() - in, expected_left_out);
  const auto triangles = mesh.triangles();
  std::set<std::pair<size_t, size_t>> edges;
  std::set<size_t> used;
  for (auto const &t : triangles) {
    Real const *a = points[t[0]].data(), *b = points[t[1]].data(), *c = points[t[2]].data();
    CHECK_TRUE((b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]) > 0);
    for (int k = 0; k < 3; ++k) {
      CHECK_TRUE(mesh.inserted(t[k]));
      CHECK_TRUE(edges.insert(std::make_pair(t[k], t[(k + 1) % 3])).second);
      used.insert(t[k]);
    }
    // The centre of the circle, from a.
    const Real bx = b[0] - a[0], by = b[1] - a[1], cx = c[0] - a[0], cy = c[1] - a[1], d = 2 * (bx * cy - by * cx);
    const Real ux = (cy * (bx * bx + by * by) - by * (cx * cx + cy * cy)) / d;
    const Real uy = (bx * (cx * cx + cy * cy) - cx * (bx * bx + by * by)) / d, radius = std::sqrt(ux * ux + uy * uy);
    for (size_t p = 0; p < points.size(); ++p) {
      const Real px = points[p][0] - a[0] - ux, py = points[p][1] - a[1] - uy;
      // A sliver has a circle as large as its rounding errors allow: only the circles of this size are tested.
      if (mesh.inserted(p) && radius < 100) {
        CHECK_GE(std::sqrt(px * px + py * py), radius * (1 - 2000 * std::numeric_limits<Real>::epsilon()));
      }
    }
  }
  size_t hull = 0;
  for (auto const &e : edges) {
    hull += edges.count(std::make_pair(e.second, e.first)) == 0;
  }
  if (!triangles.empty()) {
    CHECK_EQUAL(used.size(), in);
    CHECK_EQUAL(triangles.size() + hull + 2, 2 * in);
  }
  return triangles.size();
}

void test_delaunay() {
  using P = std::vector<std::array<double, 2>>;
  std::mt19937_64 gen(3);
  P scattered, grid, circle, line, mixed;
  for (int i = 0; i < 200; ++i) {
    scattered.push_back({unit(gen), unit(gen)});
  }
  CHECK_GE(check_delaunay(scattered, 0), size_t(300));
  std::vector<std::array<float, 2>> single;
  for (auto const &p : scattered) {
    single.push_back({static_cast<float>(p[0]), static_cast<float>(p[1])});
  }
  check_delaunay(single, 0);
  // A grid: every four neighbours are on one circle. With its points in another order as well.
  for (int i = 0; i < 49; ++i) {
    grid.push_back({double(i % 7), double(i / 7)});
  }
  CHECK_EQUAL(check_delaunay(grid, 0), size_t(72));
  for (size_t i = 0; i < grid.size(); ++i) {
    std::swap(grid[i], grid[(i * 31 + 7) % grid.size()]);
  }
  CHECK_EQUAL(check_delaunay(grid, 0), size_t(72));
  // On one circle, with its centre.
  for (int i = 0; i < 16; ++i) {
    circle.push_back({std::cos(i * 0.39269908169872414), std::sin(i * 0.39269908169872414)});
  }
  CHECK_EQUAL(check_delaunay(circle, 0), size_t(14));
  circle.push_back({0.0, 0.0});
  CHECK_EQUAL(check_delaunay(circle, 0), size_t(16));
  // On one line there is no triangle; a point off it makes one of each segment. Then points which repeat others
  // or are not finite, points on edges, on the hull and outside it.
  for (int i = 0; i < 6; ++i) {
    line.push_back({0.25 * i, 0.5 * i});
  }
  CHECK_EQUAL(check_delaunay(line, 6), size_t(0));
  line.insert(line.begin() + 2, {0.25, 0.5});
  line.push_back({1.0, 0.0});
  CHECK_EQUAL(check_delaunay(line, 1), size_t(5));
  mixed = {{0, 0}, {0, 0}, {1, 0}, {0.5, 0}, {0, 1}, {1, 0}, {std::numeric_limits<double>::quiet_NaN(), 0}, {0.5, 0.5},
           {0.25, 0.25}, {2, 2}, {0.5, 0}, {-1, 0}, {0, 0.5}, {3, 0}, {1, 1}, {std::numeric_limits<double>::infinity(), 1}};
  CHECK_EQUAL(check_delaunay(mixed, 5), size_t(13));
  // Nearly on one line: slivers, and still a triangulation.
  P slivers;
  for (int i = 0; i < 60; ++i) {
    slivers.push_back({i / 59.0, 1e-13 * ((i * 7) % 5) + (i == 30 ? 1e-9 : 0.0)});
  }
  check_delaunay(slivers, 0);
}

// The value of a matrix game, and the error of a cell against closed forms.
void test_simplex_error() {
  using detail::pareto_polish_game_value;
  using detail::pareto_polish_simplex_error;
  V b;
  CHECK_ABSOLUTE_ERROR(0.5, pareto_polish_game_value(V{0, 1, 1, 0}, 2, b), 1e-15);
  CHECK_ABSOLUTE_ERROR(2.0 / 3, pareto_polish_game_value(V{0, 2, 1, 0}, 2, b), 1e-15);
  CHECK_ABSOLUTE_ERROR(2.0 / 3, b[0], 1e-15);
  CHECK_ABSOLUTE_ERROR(1.0 / 3, pareto_polish_game_value(V{1, 0, 0, 0, 1, 0, 0, 0, 1}, 3, b), 1e-15);
  CHECK_ABSOLUTE_ERROR(3.0, pareto_polish_game_value(V{1, 2, 3, 1, 2, 3, 1, 2, 3}, 3, b), 1e-15);
  CHECK_ABSOLUTE_ERROR(1.0, b[2], 1e-15);
  // At any scale: a pivot which is small only against the ones of the constraints is not zero.
  CHECK_RELATIVE_ERROR(2e-15 / 3, pareto_polish_game_value(V{0, 1e-15, 1e-15, 1e-15, 0, 1e-15, 1e-15, 1e-15, 0}, 3, b), 1e-12);
  std::vector<float> bf;
  CHECK_RELATIVE_ERROR(2e-6f / 3, pareto_polish_game_value(std::vector<float>{0, 1e-6f, 1e-6f, 1e-6f, 0, 1e-6f, 1e-6f, 1e-6f, 0}, 3, bf), 1e-5f);
  // Two rows meet on an edge, below the third.
  CHECK_ABSOLUTE_ERROR(0.5, pareto_polish_game_value(V{0, 1, 0, 1, 0, 0, 5, 5, 5}, 3, b), 1e-15);
  CHECK_ABSOLUTE_ERROR(0.0, b[2], 1e-15);
  // A paraboloid over the plane orthogonal to r, g(u) = c |u|^2 / 2, convex or concave: for a cell with its
  // corners at the distance R from u = 0 and that point inside, the tangent planes meet at -c R^2 / 2 and the
  // cell is at c R^2 / 2. In general: the error against a search of the cell.
  const double q = std::sqrt(1.0 / 3), e1[3] = {std::sqrt(0.5), -std::sqrt(0.5), 0}, e2[3] = {q * std::sqrt(0.5), q * std::sqrt(0.5), -2 * q * std::sqrt(0.5)};
  const V r(3, q);
  std::mt19937_64 gen(17);
  for (int trial = 0; trial < 40; ++trial) {
    // The first four with the corners at R from u = 0, the second pair of them at the scale of rounding.
    const double c = trial % 2 ? 0.8 : -0.8, R = trial < 2 ? 0.2 : 1e-8;
    double u[3][2], g[3], slope[3][2];
    V Y(9), N(9);
    for (int i = 0; i < 3; ++i) {
      const double angle = trial < 4 ? 2.0943951023931953 * i + 0.3 : 6.283185307179586 * unit(gen);
      const double radius = trial < 4 ? R : 0.05 + 0.3 * unit(gen);
      u[i][0] = radius * std::cos(angle) + (trial < 4 ? 0.0 : 0.1);
      u[i][1] = radius * std::sin(angle);
      g[i] = c * (u[i][0] * u[i][0] + u[i][1] * u[i][1]) / 2;
      for (int k = 0; k < 3; ++k) {
        Y[3 * i + k] = u[i][0] * e1[k] + u[i][1] * e2[k] + g[i] * r[k];
        N[3 * i + k] = r[k] - c * (u[i][0] * e1[k] + u[i][1] * e2[k]);
      }
      slope[i][0] = c * u[i][0];
      slope[i][1] = c * u[i][1];
    }
    const auto cell = pareto_polish_simplex_error(Y, N, V{}, r, 3);
    CHECK_EQUAL(cell.kind, c > 0 ? 0 : 1);
    // The normal of the cell, for the cosine.
    double normal[3], length = 0, cosine = 0, searched = 0, truth = 0;
    for (int k = 0; k < 3; ++k) {
      const int k1 = (k + 1) % 3, k2 = (k + 2) % 3;
      normal[k] = (Y[3 + k1] - Y[k1]) * (Y[6 + k2] - Y[k2]) - (Y[3 + k2] - Y[k2]) * (Y[6 + k1] - Y[k1]);
      length += normal[k] * normal[k];
      cosine += normal[k] * r[k];
    }
    cosine = std::abs(cosine) / std::sqrt(length);
    const int steps = 120;
    for (int i = 0; i <= steps; ++i) {
      for (int j = 0; i + j <= steps; ++j) {
        const double w[3] = {double(i) / steps, double(j) / steps, double(steps - i - j) / steps};
        double x = 0, y = 0, L = 0, envelope = -1e300;
        for (int k = 0; k < 3; ++k) {
          x += w[k] * u[k][0];
          y += w[k] * u[k][1];
          L += w[k] * g[k];
        }
        for (int k = 0; k < 3; ++k) {
          const double t = g[k] + slope[k][0] * (x - u[k][0]) + slope[k][1] * (y - u[k][1]);
          envelope = (std::max)(envelope, c > 0 ? t : -t);
        }
        searched = (std::max)(searched, (c > 0 ? L : -L) - envelope);
        truth = (std::max)(truth, std::abs(L - c * (x * x + y * y) / 2));
      }
    }
    CHECK_GE(cell.error, cosine * searched * (1 - 1e-12));
    CHECK_LE(cell.error, cosine * searched * 1.05 + 1e-12);
    CHECK_GE(cell.error, cosine * truth);
    CHECK_GE(cell.dominance, cosine * truth);
    if (trial < 4) {
      CHECK_RELATIVE_ERROR(std::abs(c) * R * R, cell.error, 1e-6);
      CHECK_ABSOLUTE_ERROR(1.0 / 3, cell.b[0], 1e-6);
    }
    // The tangent of a concave front at a corner of a convex one: mixed, unless it is within the allowances.
    V wrong = N;
    for (int k = 0; k < 3; ++k) {
      wrong[k] = r[k] + c * (u[0][0] * e1[k] + u[0][1] * e2[k]);
    }
    const auto mixed = pareto_polish_simplex_error(Y, wrong, V{}, r, 3);
    if (trial >= 2) {
      continue;
    }
    CHECK_EQUAL(mixed.kind, 2);
    CHECK_EQUAL(mixed.error, 0.0);
    CHECK_EQUAL(mixed.dominance, cell.dominance);
    CHECK_TRUE(pareto_polish_simplex_error(Y, wrong, V(3, 1.0), r, 3).flat);
  }
  // The unit sphere in the positive octant, a concave front: with the corners at the angle t from r the tangent
  // planes meet at 1 / cos(t) and the cell is at cos(t).
  for (double t : {0.05, 0.3}) {
    V Y(9);
    for (int i = 0; i < 3; ++i) {
      for (int k = 0; k < 3; ++k) {
        Y[3 * i + k] = std::cos(t) * r[k] + std::sin(t) * (std::cos(2.0943951023931953 * i) * e1[k] + std::sin(2.0943951023931953 * i) * e2[k]);
      }
    }
    const auto cell = pareto_polish_simplex_error(Y, Y, V{}, r, 3);
    CHECK_EQUAL(cell.kind, 1);
    CHECK_ABSOLUTE_ERROR(1 / std::cos(t) - std::cos(t), cell.error, 1e-13);
    CHECK_GE(cell.error, 1 - std::cos(t));
    CHECK_GE(cell.dominance, 1 - std::cos(t));
  }
}

// The largest distance from the plane of a cell of the points front(p) which are over it along r, for p among the
// combinations of the parameters P of its corners; every point in objectives divided by their ranges.
template <class Front>
double true_cell_error(std::array<double, 3> const *Y, std::array<double, 2> const *P, std::vector<double> const &r, Front front) {
  auto dot = [](std::array<double, 3> const &a, std::array<double, 3> const &b) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; };
  std::array<double, 3> e1, e2, normal;
  for (int k = 0; k < 3; ++k) {
    e1[k] = Y[1][k] - Y[0][k];
    e2[k] = Y[2][k] - Y[0][k];
  }
  normal = {e1[1] * e2[2] - e1[2] * e2[1], e1[2] * e2[0] - e1[0] * e2[2], e1[0] * e2[1] - e1[1] * e2[0]};
  const double g11 = dot(e1, e1), g12 = dot(e1, e2), g22 = dot(e2, e2), det = g11 * g22 - g12 * g12, length = std::sqrt(dot(normal, normal));
  const double along = (normal[0] * r[0] + normal[1] * r[1] + normal[2] * r[2]) / length;
  double worst = 0;
  const int steps = 16;
  for (int i = -3; i <= steps + 3; ++i) {
    for (int j = -3; i + j <= steps + 3; ++j) {
      const double w[3] = {double(i) / steps, double(j) / steps, 1 - double(i + j) / steps};
      const auto y = front(w[0] * P[0][0] + w[1] * P[1][0] + w[2] * P[2][0], w[0] * P[0][1] + w[1] * P[1][1] + w[2] * P[2][1]);
      std::array<double, 3> d;
      double h = 0;
      for (int k = 0; k < 3; ++k) {
        h += (y[k] - Y[0][k]) * normal[k] / length;
      }
      for (int k = 0; k < 3; ++k) {
        d[k] = y[k] - h / along * r[k] - Y[0][k];
      }
      const double b1 = (dot(d, e1) * g22 - dot(d, e2) * g12) / det, b2 = (dot(d, e2) * g11 - dot(d, e1) * g12) / det;
      if (b1 >= 0 && b2 >= 0 && b1 + b2 <= 1 && std::isfinite(h)) {
        worst = (std::max)(worst, std::abs(h));
      }
    }
  }
  return worst;
}

// The octant of a sphere, a concave front: the first two variables are angles, the third its radius less one.
struct octant {
  template <class C> std::vector<typename C::value_type> operator()(C const &x) const {
    using std::cos;
    using std::sin;
    const auto radius = 1 + (x[2] - 0.5) * (x[2] - 0.5);
    const auto a = 1.5707963267948966 * x[0], b = 1.5707963267948966 * x[1];
    return {radius * cos(a) * cos(b), radius * cos(a) * sin(b), radius * sin(a)};
  }
};

// An objective which is zero on a face of the box, whatever the other variables: there every point is weakly
// Pareto optimal, and the front has x[2] = 1 / 2.
struct weak_face {
  template <class C> std::vector<typename C::value_type> operator()(C const &x) const {
    const auto radius = 1 + (x[2] - 0.5) * (x[2] - 0.5);
    return {radius * x[1], radius * x[0], radius * (1 - x[0])};
  }
};

// Three objectives: the cells are triangles, and where one is a bound it bounds the front, which is known.
void test_triangles() {
  for (int problem = 0; problem < 2; ++problem) {
    const size_t n = problem == 0 ? 4 : 3;
    auto params = unit_box<V>(n, 0.0);
    const auto front = problem == 0 ? scattered_front<double>(sites{}, n, 30, 7) : scattered_front<double>(octant{}, n, 30, 9);
    std::vector<pareto_polish_event<double>> trace;
    auto polish = [&]() {
      trace.clear();
      return problem == 0 ? pareto_polish(sites{}, front, params, {}, &trace) : pareto_polish(octant{}, front, params, {}, &trace);
    };
    const auto measured = polish();
    check_guarantee(front, measured, params);
    CHECK_TRUE(measured.status == status::not_refined);
    // The same on several threads, with and without the refinement.
    params.threads = several_threads;
    check_same(measured, polish());
    params.refinement_tolerance = 0.02;
    const auto result = polish();
    check_guarantee(front, result, params);
    check_trace(trace, result);
    params.threads = 1;
    check_same(result, polish());
    CHECK_TRUE(result.status == status::discovered_cells_resolved || result.status == status::cells_unresolved);
    CHECK_GE(result.new_points, size_t(10));
    CHECK_GE(result.boundary_cells, size_t(3));
    // The parameters of the points: the weights of two sites, or the angles.
    std::vector<std::array<double, 2>> parameters;
    std::vector<char> on_front;
    for (auto const &p : result) {
      if (problem == 1) {
        parameters.push_back({p.x[0], p.x[1]});
        on_front.push_back(std::abs(p.x[2] - 0.5) < 1e-4);
        continue;
      }
      double a11 = 0, a12 = 0, a22 = 0, b1 = 0, b2 = 0, off = 0;
      for (size_t k = 0; k < 4; ++k) {
        const double c1 = site[0][k] - site[2][k], c2 = site[1][k] - site[2][k], c0 = p.x[k] - site[2][k];
        a11 += c1 * c1;
        a12 += c1 * c2;
        a22 += c2 * c2;
        b1 += c1 * c0;
        b2 += c2 * c0;
      }
      const double det = a11 * a22 - a12 * a12, w1 = (b1 * a22 - b2 * a12) / det, w2 = (a11 * b2 - a12 * b1) / det;
      for (size_t k = 0; k < 4; ++k) {
        off = (std::max)(off, std::abs(w1 * site[0][k] + w2 * site[1][k] + (1 - w1 - w2) * site[2][k] - p.x[k]));
      }
      parameters.push_back({w1, w2});
      on_front.push_back(off < 1e-4);
    }
    auto truth = [&](double p1, double p2) {
      V x(n, 0.5);
      const double nan = std::numeric_limits<double>::quiet_NaN();
      if (problem == 1) {
        x[0] = p1;
        x[1] = p2;
      } else {
        for (size_t k = 0; k < 4; ++k) {
          x[k] = p1 * site[0][k] + p2 * site[1][k] + (1 - p1 - p2) * site[2][k];
        }
      }
      const bool inside = problem == 1 ? (p1 >= 0 && p1 <= 1 && p2 >= 0 && p2 <= 1) : (p1 >= 0 && p2 >= 0 && p1 + p2 <= 1);
      const V y = problem == 0 ? sites{}(x) : octant{}(x);
      return std::array<double, 3>{inside ? y[0] / result.ranges[0] : nan, y[1] / result.ranges[1], y[2] / result.ranges[2]};
    };
    size_t bounds = 0, triangles = 0, tested = 0;
    for (auto const &c : result.cells) {
      CHECK_TRUE(std::isfinite(c.error) && c.error >= 0);
      if (c.vertices.size() != 3) {
        CHECK_TRUE(c.kind == cell_class::boundary || c.kind == cell_class::unresolved);
        CHECK_TRUE(!c.bound_under_convexity);
        continue;
      }
      ++triangles;
      CHECK_TRUE(c.vertices[0] != c.vertices[1] && c.vertices[1] != c.vertices[2] && c.vertices[0] != c.vertices[2]);
      CHECK_TRUE(std::isfinite(c.dominance_bound));
      if (!c.bound_under_convexity) {
        continue;
      }
      ++bounds;
      // The front of the squared distances is convex; the octant is concave.
      CHECK_TRUE(c.kind == (problem == 0 ? cell_class::convex : cell_class::concave));
      CHECK_LE(c.error, params.refinement_tolerance);
      std::array<double, 3> Y[3];
      std::array<double, 2> P[3];
      bool known = true;
      for (int k = 0; k < 3; ++k) {
        auto const &p = result[c.vertices[k]];
        Y[k] = {p.f[0] / result.ranges[0], p.f[1] / result.ranges[1], p.f[2] / result.ranges[2]};
        P[k] = parameters[c.vertices[k]];
        known = known && on_front[c.vertices[k]];
      }
      if (known) {
        ++tested;
        CHECK_LE(true_cell_error(Y, P, result.chart_direction, truth), c.error * (1 + 1e-6) + 1e-7);
      }
    }
    // Measured: 229 triangles, 190 of them bounds, for the sites; 262 and 232 for the octant.
    CHECK_GE(triangles, size_t(40));
    CHECK_GE(2 * bounds, triangles);
    CHECK_GE(2 * tested, bounds);
  }
  // A point which is only weakly optimal, on a face where one objective is zero: the steps of the second kind
  // leave it for the front.
  {
    auto params = unit_box<V>(3, 0.0);
    const auto front = front_at(weak_face{}, {V{0.3, 0.0, 0.8}});
    const auto result = pareto_polish(weak_face{}, front, params);
    check_guarantee(front, result, params);
    CHECK_EQUAL(result.size(), size_t(1));
    CHECK_EQUAL(result[0].f[0], 0.0);
    CHECK_ABSOLUTE_ERROR(0.5, result[0].x[2], 1e-3);
  }
  // The same with two objectives: the first is zero on the face x[0] = 0, and its gradient vanishes there, so
  // that the step of the second kind leaves the face and is never acceptable.
  {
    auto weak = [](auto const &x) {
      using T = typename std::decay<decltype(x[0])>::type;
      const T radius = 1 + (x[1] - 0.5) * (x[1] - 0.5);
      return std::vector<T>{radius * x[0] * x[0], radius * (1 - x[0])};
    };
    auto params = unit_box<V>(2, 0.0);
    const auto front = front_at(weak, {V{0.0, 0.8}});
    const auto result = pareto_polish(weak, front, params);
    check_guarantee(front, result, params);
    CHECK_EQUAL(result[0].f[0], 0.0);
    CHECK_ABSOLUTE_ERROR(0.5, result[0].x[1], 1e-3);
  }
  // linearized_nonworsening_gain, against its value in closed form. Inside the triangle of the sites every
  // multiplier is positive and it is zero. Between two of the sites the multiplier of the third is zero: the
  // gain is the length of the gradient of the third objective less its part along the gradients of the other
  // two, which are parallel. (The range of a front of one point is one.)
  {
    auto params = unit_box<V>(4, 0.0);
    for (int where = 0; where < 2; ++where) {
      V x(4);
      double along = 0, length2 = 0, edge2 = 0;
      for (size_t k = 0; k < 4; ++k) {
        x[k] = where == 0 ? 0.5 * site[0][k] + 0.3 * site[1][k] + 0.2 * site[2][k] : 0.6 * site[0][k] + 0.4 * site[1][k];
        along += (x[k] - site[2][k]) * (site[1][k] - site[0][k]);
        length2 += (x[k] - site[2][k]) * (x[k] - site[2][k]);
        edge2 += (site[1][k] - site[0][k]) * (site[1][k] - site[0][k]);
      }
      const auto result = pareto_polish(sites{}, front_at(sites{}, {x}), params);
      CHECK_TRUE(result[0].x == x);
      CHECK_LE(result[0].weak_kkt_residual, 1e-12);
      if (where == 0) {
        CHECK_LE(result[0].linearized_nonworsening_gain, 1e-9);
      } else {
        CHECK_ABSOLUTE_ERROR(2 * std::sqrt(length2 - along * along / edge2), result[0].linearized_nonworsening_gain, 1e-9);
      }
    }
  }
  // Slivers and less: points of the front of two of the objectives, which are on one curve; one point; two; and a
  // point whose objectives are not finite among the others.
  for (size_t count : {size_t(1), size_t(2), size_t(6)}) {
    auto params = unit_box<V>(4, 0.05);
    std::vector<std::pair<V, V>> front;
    for (size_t i = 0; i < count; ++i) {
      const double t = count > 1 ? double(i) / double(count - 1) : 0.4;
      V x(4);
      for (size_t k = 0; k < 4; ++k) {
        x[k] = (1 - t) * site[0][k] + t * site[1][k];
      }
      front.emplace_back(x, sites{}(x));
    }
    if (count == 6) {
      front[2].second[1] = std::numeric_limits<double>::quiet_NaN();
    }
    const auto result = pareto_polish(sites{}, front, params);
    check_guarantee(front, result, params);
    for (auto const &c : result.cells) {
      CHECK_TRUE(std::isfinite(c.error) && c.error >= 0);
    }
    params.max_new_points = 4;
    check_guarantee(front, pareto_polish(sites{}, front, params), params);
  }
  // No triangle: two points, and points on one line of the space of the objectives. The cells are edges, with
  // estimates, as for more objectives; the status is not no_cells.
  {
    auto line = [](auto const &x) {
      using T = typename std::decay<decltype(x[0])>::type;
      return std::vector<T>{x[0] + x[1] * x[1], 1 - x[0] + x[1] * x[1], 0 * x[0] + 0.5};
    };
    for (size_t count : {size_t(2), size_t(5)}) {
      for (double tolerance : {0.0, 0.05}) {
        auto params = unit_box<V>(2, tolerance);
        std::vector<std::pair<V, V>> front;
        for (size_t i = 0; i < count; ++i) {
          const V x = {0.1 + 0.8 * double(i) / double(count - 1), 0.0};
          front.emplace_back(x, line(x));
        }
        const auto result = pareto_polish(line, front, params);
        check_guarantee(front, result, params);
        CHECK_GE(result.cells.size() + 1, result.size());
        CHECK_TRUE(result.status != status::no_cells);
        for (auto const &c : result.cells) {
          CHECK_EQUAL(c.vertices.size(), size_t(2));
          CHECK_TRUE(!c.bound_under_convexity);
        }
      }
    }
  }
  // A point whose projection repeats that of a corner is in no triangle: an unresolved edge joins it to the corner.
  {
    using detail::pareto_polish_edge;
    using detail::pareto_polish_state;
    detail::pareto_polish_triangulation<double> mesh;
    std::vector<pareto_polish_state<double>> states(5);
    const double r = std::sqrt(1.0 / 3), y[5][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {0.4, 0.3, 0.3}, {0.4 + 0.1 * r, 0.3 + 0.1 * r, 0.3 + 0.1 * r}};
    for (size_t i = 0; i < 5; ++i) {
      states[i].f.assign(y[i], y[i] + 3);
      states[i].x = states[i].f;
      states[i].lam = V(3, 1.0 / 3);
      states[i].id = i;
    }
    std::vector<pareto_polish_edge<double>> cells;
    mesh.cells(states, {0, 1, 2, 3, 4}, V(3, 1.0), 0.01, 0.25, 1e300, cells);
    size_t triangles = 0, joined = 0;
    for (auto const &c : cells) {
      triangles += c.k != pareto_polish_created;
      if (c.k == pareto_polish_created && c.kind == cell_class::unresolved) {
        ++joined;
        CHECK_TRUE(c.i == 4 && c.j == 3);
        CHECK_ABSOLUTE_ERROR(0.1, c.error, 1e-12);
      }
    }
    CHECK_EQUAL(triangles, size_t(3));
    CHECK_EQUAL(joined, size_t(1));
  }
  // max_iterations = 0: no step is taken, and the points are certified as they are.
  {
    auto params = unit_box<V>(4, 0.0);
    params.max_iterations = 0;
    const auto front = scattered_front<double>(sites{}, 4, 5, 7);
    const auto result = pareto_polish(sites{}, front, params);
    for (auto const &p : result) {
      CHECK_TRUE(p.termination == termination::iteration_cap);
      CHECK_TRUE(p.x == front[p.origin].first);
      CHECK_EQUAL(p.objective_evaluations, size_t(0));
      CHECK_TRUE(std::isfinite(p.weak_kkt_residual));
    }
  }
  // In single precision.
  {
    auto params = unit_box<std::vector<float>>(4, 0.05f);
    params.threads = several_threads;
    const auto front = scattered_front<float>(sites{}, 4, 20, 3);
    const auto result = pareto_polish(sites{}, front, params);
    check_guarantee(front, result, params);
    CHECK_GE(result.new_points, size_t(1));
    CHECK_GE(result.boundary_cells, size_t(3));
  }
  // A cell without area, corners which repeat one another, and corners which give no tangent plane.
  using detail::pareto_polish_simplex_error;
  const V r(3, std::sqrt(1.0 / 3)), line = {0.3, 0.2, 0.1, 0.2, 0.2, 0.2, 0.1, 0.2, 0.3}, normals = {1, 1, 1, 1, 2, 1, 1, 1, 2};
  for (auto const &Y : {line, V{0.3, 0.2, 0.1, 0.3, 0.2, 0.1, 0.1, 0.2, 0.3}, V(9, 0.25)}) {
    const auto cell = pareto_polish_simplex_error(Y, normals, V{}, r, 3);
    CHECK_TRUE(std::isfinite(cell.error) && cell.error >= 0 && std::isfinite(cell.dominance) && cell.dominance >= 0);
    CHECK_ABSOLUTE_ERROR(1.0, cell.b[0] + cell.b[1] + cell.b[2], 1e-12);
  }
  const V Y = {0.5, 0.1, 0.1, 0.1, 0.5, 0.1, 0.1, 0.1, 0.5}, N = {1, 2, 2, 2, 1, 2, 2, 2, 1};
  const auto all = pareto_polish_simplex_error(Y, N, V{}, r, 3);
  const auto two = pareto_polish_simplex_error(Y, N, V{}, r, 3, std::vector<char>{0, 0, 1});
  const auto none = pareto_polish_simplex_error(Y, N, V{}, r, 3, std::vector<char>(3, 1));
  CHECK_EQUAL(all.kind, 0);
  CHECK_EQUAL(two.kind, 0);
  CHECK_GE(two.error, all.error);
  CHECK_EQUAL(none.kind, 2);
  CHECK_EQUAL(none.dominance, all.dominance);
}

int main() {
  for (size_t m = 2; m <= 5; ++m) {
    test_subproblem<double>(m, false);
    test_subproblem<double>(m, true);
  }
  test_subproblem<float>(3, false);
  test_sandwich();
  test_cells();
  test_failed_search_is_no_gap();
  test_parameters();
  test_curved_set();
  test_bounds();
  test_edge_cases();
  test_start_not_finite();
  test_dominated_inputs();
  test_corrector();
  test_corrector_indefinite();
  test_threads_and_tapes();
  test_curvature();
  test_svd();
  test_delaunay();
  test_simplex_error();
  test_triangles();
  return boost::math::test::report_errors();
}
