/*
 * Copyright Nick Thompson, 2026
 * Use, modification and distribution are subject to the
 * Boost Software License, Version 1.0. (See accompanying file
 * LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef BOOST_MATH_OPTIMIZATION_DETAIL_PARETO_POLISH_CELLS_HPP
#define BOOST_MATH_OPTIMIZATION_DETAIL_PARETO_POLISH_CELLS_HPP
#include <boost/math/tools/config.hpp>
#include <boost/math/tools/precision.hpp>
#include <boost/math/optimization/detail/pareto_polish_subproblem.hpp>
#ifndef BOOST_MATH_BUILD_MODULE
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <map>
#include <utility>
#include <vector>
#endif

BOOST_MATH_NAMESPACE_BEGIN namespace optimization { namespace detail {

// The singular value decomposition A = U' diag(sigma) V of the m by n matrix A, by the one sided Jacobi method of
// Hestenes on its rows: U is m by m and orthogonal, the singular values decrease, and row k of the m by n matrix V
// is the right singular vector of sigma[k], or zero where that is zero. The rows of A are only ever rotated
// against one another, so the small singular values of a matrix with rows of very different lengths keep their
// relative accuracy, which the eigenvalues of A A' do not have.
BOOST_MATH_TEST_EXPORT template <typename Real>
void pareto_polish_svd(std::vector<Real> A, size_t m, size_t n, std::vector<Real> &U, std::vector<Real> &sigma,
                       std::vector<Real> &V) {
  using std::abs;
  using std::sqrt;
  Real scale = 0;
  for (auto const &a : A) {
    scale = (std::max)(scale, abs(a));
  }
  scale = scale > 0 && scale < std::numeric_limits<Real>::infinity() ? scale : Real(1);
  for (auto &a : A) {
    a /= scale;
  }
  std::vector<Real> W(m * m, Real(0));
  for (size_t i = 0; i < m; ++i) {
    W[i * m + i] = 1;
  }
  auto rotate = [](Real *p, Real *q, size_t count, Real c, Real s) {
    for (size_t k = 0; k < count; ++k) {
      const Real a = p[k], b = q[k];
      p[k] = c * a - s * b;
      q[k] = s * a + c * b;
    }
  };
  bool rotated = true;
  for (int sweep = 0; sweep < 60 && rotated; ++sweep) {
    rotated = false;
    for (size_t p = 0; p + 1 < m; ++p) {
      for (size_t q = p + 1; q < m; ++q) {
        Real alpha = 0, beta = 0, gamma = 0;
        for (size_t k = 0; k < n; ++k) {
          alpha += A[p * n + k] * A[p * n + k];
          beta += A[q * n + k] * A[q * n + k];
          gamma += A[p * n + k] * A[q * n + k];
        }
        if (!(abs(gamma) > std::numeric_limits<Real>::epsilon() * sqrt(alpha * beta))) {
          continue;
        }
        const Real zeta = (beta - alpha) / (2 * gamma);
        const Real t = (zeta < 0 ? Real(-1) : Real(1)) / (abs(zeta) + sqrt(1 + zeta * zeta));
        const Real c = 1 / sqrt(1 + t * t);
        rotate(&A[p * n], &A[q * n], n, c, c * t);
        rotate(&W[p * m], &W[q * m], m, c, c * t);
        rotated = true;
      }
    }
  }
  std::vector<std::pair<Real, size_t>> order(m);
  for (size_t i = 0; i < m; ++i) {
    Real norm2 = 0;
    for (size_t k = 0; k < n; ++k) {
      norm2 += A[i * n + k] * A[i * n + k];
    }
    order[i] = {-sqrt(norm2), i};
  }
  std::sort(order.begin(), order.end());
  U.assign(m * m, Real(0));
  V.assign(m * n, Real(0));
  sigma.assign(m, Real(0));
  for (size_t i = 0; i < m; ++i) {
    const size_t from = order[i].second;
    const Real length = -order[i].first;
    sigma[i] = scale * length;
    std::copy(&W[from * m], &W[from * m] + m, &U[i * m]);
    for (size_t k = 0; k < n && length > 0; ++k) {
      V[i * n + k] = A[from * n + k] / length;
    }
  }
}

// The Delaunay triangulation of points of the plane, one added at a time by the method of Bowyer and Watson. The
// edges of the convex hull carry triangles to a vertex at infinity (`ghost`), so that a point outside the hull is
// added as one inside is. The predicates are evaluated in Real, and are not exact: the triangles removed for a
// new point are connected by construction and kept visible from it, so the result is always a triangulation,
// and it is Delaunay up to rounding. A point which repeats another, or is not finite, is left out; so are all of
// them while they lie on one line. Finding the triangle of a new point is a search of all of them: N points take
// O(N^2) operations.
BOOST_MATH_TEST_EXPORT template <typename Real> class pareto_polish_delaunay {
public:
  static size_t ghost() { return static_cast<size_t>(-1); }
  // The points are numbered in the order of the calls.
  void insert(Real x, Real y) {
    using std::isfinite;
    const size_t p = in_.size();
    xy_.push_back(x);
    xy_.push_back(y);
    in_.push_back(0);
    if (!isfinite(x) || !isfinite(y)) {
      return;
    }
    if (!triangles_.empty()) {
      in_[p] = add(p);
      return;
    }
    // The first triangle: of the first point, the first which differs from it and the first off their line.
    waiting_.push_back(p);
    const size_t a = waiting_[0];
    size_t b = a;
    for (size_t q : waiting_) {
      if (b == a && (xy_[2 * q] != xy_[2 * a] || xy_[2 * q + 1] != xy_[2 * a + 1])) {
        b = q;
      }
    }
    // On one line, to rounding: the sine of the angle at a is below ten thousand times the precision.
    using std::abs;
    using std::sqrt;
    Real ab = 0, ap = 0;
    for (size_t c = 0; c < 2; ++c) {
      ab += (xy_[2 * b + c] - xy_[2 * a + c]) * (xy_[2 * b + c] - xy_[2 * a + c]);
      ap += (xy_[2 * p + c] - xy_[2 * a + c]) * (xy_[2 * p + c] - xy_[2 * a + c]);
    }
    if (b == a || b == p || !(abs(orient(a, b, p)) > 10000 * std::numeric_limits<Real>::epsilon() * sqrt(ab * ap))) {
      return;
    }
    const size_t v[3] = {a, orient(a, b, p) > 0 ? b : p, orient(a, b, p) > 0 ? p : b};
    create(v[0], v[1], v[2]);
    for (int k = 0; k < 3; ++k) {
      create(v[(k + 1) % 3], v[k], ghost());
    }
    in_[a] = in_[b] = in_[p] = 1;
    for (size_t q : waiting_) {
      if (q != a && q != b && q != p) {
        in_[q] = add(q);
      }
    }
    waiting_.clear();
  }
  bool inserted(size_t p) const { return in_[p] != 0; }
  // The triangles, each counter-clockwise.
  std::vector<std::array<size_t, 3>> triangles() const {
    std::vector<std::array<size_t, 3>> out;
    for (auto const &t : triangles_) {
      if (t.alive && t.v[2] != ghost()) {
        out.push_back({{t.v[0], t.v[1], t.v[2]}});
      }
    }
    return out;
  }

private:
  struct triangle {
    size_t v[3]; // counter-clockwise; a ghost is (a, b, ghost), the hull to the right of a -> b
    bool alive;
  };
  Real orient(size_t a, size_t b, size_t c) const {
    return (xy_[2 * b] - xy_[2 * a]) * (xy_[2 * c + 1] - xy_[2 * a + 1]) -
           (xy_[2 * b + 1] - xy_[2 * a + 1]) * (xy_[2 * c] - xy_[2 * a]);
  }
  // p is inside the circle through the corners of t; for a ghost, outside its edge of the hull or within it.
  bool conflict(triangle const &t, size_t p) const {
    Real d[3][2] = {{Real(0), Real(0)}, {Real(0), Real(0)}, {Real(0), Real(0)}};
    for (int k = 0; k < 3 && t.v[k] != ghost(); ++k) {
      d[k][0] = xy_[2 * t.v[k]] - xy_[2 * p];
      d[k][1] = xy_[2 * t.v[k] + 1] - xy_[2 * p + 1];
    }
    if (t.v[2] == ghost()) {
      const Real o = orient(t.v[0], t.v[1], p);
      return o > 0 || (o == 0 && d[0][0] * d[1][0] + d[0][1] * d[1][1] < 0);
    }
    Real det = 0;
    for (int k = 0; k < 3; ++k) {
      Real const *a = d[k], *b = d[(k + 1) % 3], *c = d[(k + 2) % 3];
      det += (a[0] * a[0] + a[1] * a[1]) * (b[0] * c[1] - b[1] * c[0]);
    }
    return det > 0;
  }
  void create(size_t a, size_t b, size_t c) {
    // The ghost vertex comes last.
    const size_t v[3] = {a == ghost() ? b : (b == ghost() ? c : a), a == ghost() ? c : (b == ghost() ? a : b),
                         a == ghost() ? a : (b == ghost() ? b : c)};
    for (int k = 0; k < 3; ++k) {
      owner_[std::make_pair(v[k], v[(k + 1) % 3])] = triangles_.size();
    }
    triangles_.push_back({{v[0], v[1], v[2]}, true});
  }
  bool add(size_t p) {
    if (2 * dead_ > triangles_.size()) {
      // The triangles which are no more are dropped, so that the searches are over those which are.
      std::vector<triangle> alive;
      owner_.clear();
      dead_ = 0;
      alive.swap(triangles_);
      for (auto const &t : alive) {
        if (t.alive) {
          create(t.v[0], t.v[1], t.v[2]);
        }
      }
    }
    const size_t none = triangles_.size();
    size_t seed = none;
    for (size_t t = 0; t < none; ++t) {
      triangle const &T = triangles_[t];
      if (!T.alive) {
        continue;
      }
      if (T.v[2] == ghost()) {
        seed = seed == none && conflict(T, p) ? t : seed;
      } else if (orient(T.v[0], T.v[1], p) >= 0 && orient(T.v[1], T.v[2], p) >= 0 && orient(T.v[2], T.v[0], p) >= 0) {
        seed = t;
        break;
      }
    }
    if (seed == none) {
      return false;
    }
    // It repeats a corner, to rounding (ten thousand times the precision of the coordinates).
    for (size_t v : triangles_[seed].v) {
      const Real dx = v != ghost() ? xy_[2 * v] - xy_[2 * p] : Real(1);
      const Real dy = v != ghost() ? xy_[2 * v + 1] - xy_[2 * p + 1] : Real(1);
      const Real size = xy_[2 * p] * xy_[2 * p] + xy_[2 * p + 1] * xy_[2 * p + 1];
      const Real tiny = 10000 * std::numeric_limits<Real>::epsilon();
      if (v != ghost() && dx * dx + dy * dy <= tiny * tiny * size) {
        return false;
      }
    }
    // The cavity: the triangles in conflict with p which are connected to the seed. A triangle with an edge of
    // the cavity which p does not see is left out of it (the seed takes its neighbour in instead).
    std::vector<char> state(none, 0), member(none, 0); // 1: in the cavity, 2: left out, 3: taken in
    std::vector<size_t> cavity;
    std::vector<std::pair<size_t, size_t>> edges;
    for (bool valid = false; !valid;) {
      cavity.assign(1, seed);
      std::fill(member.begin(), member.end(), 0);
      member[seed] = 1;
      for (auto &s : state) {
        s = s == 1 ? 0 : s;
      }
      state[seed] = 1;
      edges.clear();
      valid = true;
      for (size_t q = 0; q < cavity.size() && valid; ++q) {
        size_t const *v = triangles_[cavity[q]].v;
        for (int k = 0; k < 3 && valid; ++k) {
          const size_t a = v[k], b = v[(k + 1) % 3];
          auto const across = owner_.find(std::make_pair(b, a));
          if (across == owner_.end()) {
            return false;
          }
          const size_t t = across->second;
          if (state[t] == 0 && conflict(triangles_[t], p)) {
            state[t] = 1;
          }
          if (state[t] == 1 || state[t] == 3) {
            if (!member[t]) {
              member[t] = 1;
              cavity.push_back(t);
            }
            continue;
          }
          if (a != ghost() && b != ghost() && !(orient(a, b, p) > 0)) {
            if (cavity[q] != seed) {
              state[cavity[q]] = 2;
            } else if (state[t] == 0) {
              state[t] = 3;
            } else {
              return false;
            }
            valid = false;
          }
          edges.emplace_back(a, b);
        }
      }
    }
    for (size_t t : cavity) {
      triangle &T = triangles_[t];
      T.alive = false;
      ++dead_;
      for (int k = 0; k < 3; ++k) {
        owner_.erase(std::make_pair(T.v[k], T.v[(k + 1) % 3]));
      }
    }
    for (auto const &e : edges) {
      create(e.first, e.second, p);
    }
    return true;
  }
  std::vector<Real> xy_;
  std::vector<char> in_;
  std::vector<size_t> waiting_;
  std::vector<triangle> triangles_;
  size_t dead_ = 0;
  std::map<std::pair<size_t, size_t>, size_t> owner_;
};

// The value of the matrix game of the matrix A of k columns and at most as many rows: the largest over b in the
// simplex of the least entry of A b, with a b for it. The value is between the least entry of A b for any b in the
// simplex and the largest entry of p A for any p in the simplex of the rows; both are tried at the points where
// the entries in a set of rows (columns) are equal and b (p) vanishes outside as many columns (rows), for A
// divided by its largest entry. What is returned is the least of the upper bounds, so that it is never below the
// value, whatever rounding makes of the systems, and the b of the largest lower bound.
BOOST_MATH_TEST_EXPORT template <typename Real>
Real pareto_polish_game_value(std::vector<Real> A, size_t k, std::vector<Real> &b) {
  using std::abs;
  const size_t count = A.size() / k;
  const Real infinity = std::numeric_limits<Real>::infinity();
  Real lower = -infinity, upper = infinity, scale = 0;
  for (auto const &v : A) {
    scale = (std::max)(scale, abs(v));
  }
  b.assign(k, Real(1) / static_cast<Real>(k));
  if (!(scale > 0) || !(scale < infinity)) {
    return scale > 0 ? scale : Real(0);
  }
  for (auto &v : A) {
    v /= scale;
  }
  std::vector<Real> a, sol, trial;
  std::vector<size_t> I, K;
  // The strategy `trial` of the columns (the rows, if `dual`): the bound it gives.
  auto bound = [&](bool dual) {
    Real sum = 0, extreme = dual ? -infinity : infinity;
    for (auto &t : trial) {
      t = (std::max)(t, Real(0));
      sum += t;
    }
    if (!(sum > 0) || !(sum < infinity)) {
      return;
    }
    for (size_t i = 0; i < (dual ? k : count); ++i) {
      Real v = 0;
      for (size_t j = 0; j < (dual ? count : k); ++j) {
        v += (dual ? A[j * k + i] : A[i * k + j]) * trial[j] / sum;
      }
      extreme = dual ? (std::max)(extreme, v) : (std::min)(extreme, v);
    }
    if (dual) {
      upper = (std::min)(upper, extreme);
    } else if (extreme > lower) {
      lower = extreme;
      b = trial;
      for (auto &t : b) {
        t /= sum;
      }
    }
  };
  for (size_t rows = 1; rows < (size_t(1) << count); ++rows) {
    for (size_t columns = 1; columns < (size_t(1) << k); ++columns) {
      I.clear();
      K.clear();
      for (size_t i = 0; i < k; ++i) {
        if (i < count && (rows >> i & 1)) {
          I.push_back(i);
        }
        if (columns >> i & 1) {
          K.push_back(i);
        }
      }
      const size_t s = K.size(), width = s + 2;
      if (I.size() != s) {
        continue;
      }
      for (int dual = 0; dual < 2; ++dual) {
        // Unknowns b_K and v: (A b)_i - v = 0 for i in I, and the sum of b is one. Or p_I and v, with A'.
        a.assign((s + 1) * width, Real(0));
        for (size_t i = 0; i < s; ++i) {
          for (size_t j = 0; j < s; ++j) {
            a[i * width + j] = dual ? A[I[j] * k + K[i]] : A[I[i] * k + K[j]];
          }
          a[i * width + s] = -1;
          a[s * width + i] = 1;
        }
        a[s * width + s + 1] = 1;
        trial.assign(dual ? count : k, Real(0));
        if (s == 1) {
          trial[dual ? I[0] : K[0]] = 1;
        } else if (pareto_polish_solve(a, s + 1, sol)) {
          for (size_t j = 0; j < s; ++j) {
            trial[dual ? I[j] : K[j]] = sol[j];
          }
        } else {
          continue;
        }
        bound(dual != 0);
      }
    }
  }
  return scale * upper;
}

// What is known of a cell with k corners in the space of k objectives divided by their ranges.
BOOST_MATH_TEST_EXPORT template <typename Real> struct pareto_polish_simplex {
  // The largest distance from the plane of the cell of a point between that plane and the tangent planes of the
  // corners, over the cell along r: see pareto_polish_simplex_error. Zero for a mixed cell.
  Real error = 0;
  // The same for a point which no corner dominates and which dominates none: no convexity is needed.
  Real dominance = 0;
  int kind = 2;        // 0: convex, 1: concave, 2: mixed
  bool flat = false;   // every height is within the allowances: convex and concave alike
  std::vector<Real> b; // the barycentric coordinates of the point of the cell under (or over) the farthest
};

// The tangent envelope of a cell: Y holds its k corners, row by row; N a normal of the front at each (multipliers
// of the objectives divided by their ranges, not all zero and none negative); e the allowance of each tangent
// plane, as a distance; r a unit vector with positive entries. The front is taken as a graph over the plane
// orthogonal to r, y(u) = u + g(u) r. With L the affine function through the corners and t_i the tangent function
// of corner i, a convex g has max_i t_i <= g <= L over the cell, so that g is within
//   E = max over the cell of (L - max_i (t_i - e_i))
// of L; for a concave g the envelopes change places. In barycentric coordinates b, L - t_i is the affine function
// with the values A_ij = (height of corner j over the tangent plane of corner i, along r) at the corners, and E is
// the value of the matrix game of A. The error is E times the cosine of the angle of r with the normal of the
// cell. A corner with skip set gives no tangent plane; with every corner skipped the cell is mixed.
BOOST_MATH_TEST_EXPORT template <typename Real>
pareto_polish_simplex<Real> pareto_polish_simplex_error(std::vector<Real> const &Y, std::vector<Real> const &N,
                                                        std::vector<Real> const &e, std::vector<Real> const &r,
                                                        size_t k, std::vector<char> const &skip = std::vector<char>()) {
  using std::abs;
  using std::sqrt;
  pareto_polish_simplex<Real> out;
  std::vector<Real> H(k * k, Real(0)), A, along(k, Real(0)), allowance(k), Q, v(k);
  bool convex = true, concave = true;
  Real sum = 0;
  auto used = [&](size_t i) { return i >= skip.size() || !skip[i]; };
  for (size_t i = 0; i < k; ++i) {
    if (!used(i)) {
      continue;
    }
    Real norm = 0;
    for (size_t c = 0; c < k; ++c) {
      norm += N[i * k + c] * N[i * k + c];
      along[i] += N[i * k + c] * r[c];
    }
    norm = sqrt(norm);
    along[i] /= norm;
    allowance[i] = i < e.size() ? e[i] : Real(0);
    for (size_t j = 0; j < k; ++j) {
      for (size_t c = 0; c < k && j != i; ++c) {
        H[i * k + j] += N[i * k + c] * (Y[j * k + c] - Y[i * k + c]) / norm;
      }
      convex = convex && H[i * k + j] >= -allowance[i];
      concave = concave && H[i * k + j] <= allowance[i];
      sum += H[i * k + j];
    }
  }
  // The cosine: r less its part along the edges of the cell, which are made orthonormal first.
  v = r;
  for (size_t j = 1; j < k; ++j) {
    std::vector<Real> q(k);
    for (size_t c = 0; c < k; ++c) {
      q[c] = Y[j * k + c] - Y[c];
    }
    Real norm = 0, start = 0;
    for (int twice = 0; twice < 2; ++twice) {
      for (size_t p = 0; p * k < Q.size(); ++p) {
        Real dot = 0;
        for (size_t c = 0; c < k; ++c) {
          dot += q[c] * Q[p * k + c];
        }
        for (size_t c = 0; c < k; ++c) {
          q[c] -= dot * Q[p * k + c];
        }
      }
      norm = 0;
      for (size_t c = 0; c < k; ++c) {
        norm += q[c] * q[c];
      }
      start = twice == 0 ? norm : start;
    }
    if (!(norm > start / 1000000) || !(norm > 0)) {
      continue; // a cell without volume: the cosine is taken as larger than it is
    }
    Real dot = 0;
    for (size_t c = 0; c < k; ++c) {
      q[c] /= sqrt(norm);
      dot += v[c] * q[c];
    }
    for (size_t c = 0; c < k; ++c) {
      v[c] -= dot * q[c];
    }
    Q.insert(Q.end(), q.begin(), q.end());
  }
  Real cosine = 0;
  for (size_t c = 0; c < k; ++c) {
    cosine += v[c] * v[c];
  }
  cosine = (std::min)(sqrt(cosine), Real(1));
  // The dominance bound: with w = sum_j b_j (y_j - y_i), a point w + d r from corner i is dominated by it unless
  // d < max_c(-w_c / r_c), and dominates it unless d > min_c(-w_c / r_c). Each is convex or concave in b, and so
  // is bounded by its values at the corners; the least over i is taken.
  Real over = std::numeric_limits<Real>::infinity(), under = over;
  for (size_t i = 0; i < k; ++i) {
    Real most_over = 0, most_under = 0;
    for (size_t j = 0; j < k; ++j) {
      for (size_t c = 0; c < k && j != i; ++c) {
        const Real d = (Y[i * k + c] - Y[j * k + c]) / r[c];
        most_over = (std::max)(most_over, d);
        most_under = (std::max)(most_under, -d);
      }
    }
    over = (std::min)(over, most_over);
    under = (std::min)(under, most_under);
  }
  out.dominance = cosine * (std::max)(over, under);
  out.b.assign(k, Real(1) / static_cast<Real>(k));
  if (!convex && !concave) {
    return out;
  }
  out.flat = convex && concave;
  out.kind = convex && (!concave || sum >= 0) ? 0 : 1;
  for (size_t i = 0; i < k; ++i) {
    for (size_t j = 0; j < k && used(i); ++j) {
      A.push_back(((out.kind == 0 ? H[i * k + j] : -H[i * k + j]) + allowance[i]) / along[i]);
    }
  }
  if (A.empty()) {
    out.kind = 2;
    return out;
  }
  out.error = cosine * (std::max)(pareto_polish_game_value(A, k, out.b), Real(0));
  return out;
}

// Three objectives. The multipliers of a point are not unique where the Jacobian on its free variables has not
// the rank two: on the boundary of the front they fill a segment of the simplex, at a corner of it a piece of
// the simplex itself. The segment is along the combination of the two left singular vectors of the least
// singular values whose entries sum to zero; its ends, as far as the weak residual stays within the slack, are
// kept for the cells, each of which takes the tangent plane which faces it. The same search in six directions
// of the plane of the simplex tells whether the multipliers fill more than a segment.
template <class State, typename Real>
void pareto_polish_faces(State &s, std::vector<Real> const &Jh, std::vector<Real> const &lo,
                         std::vector<Real> const &hi, std::vector<Real> const &U, Real tolerance, size_t n_) {
  using Vec = std::vector<Real>;
  using std::abs;
  using std::sqrt;
  auto residual = [&](Vec const &l) {
    return pareto_polish_residual(
        [&](size_t k) { return l[0] * Jh[k] + l[1] * Jh[n_ + k] + l[2] * Jh[2 * n_ + k]; }, lo, hi);
  };
  Real largest = 0, size = 0;
  for (auto const &v : Jh) {
    size += v * v;
  }
  for (size_t i = 0; i < 3; ++i) {
    Vec unit(3, Real(0));
    unit[i] = 1;
    largest = (std::max)(largest, residual(unit));
  }
  // The residual which a multiplier may have: that at which a point counts as stationary. The noise is that of
  // the whole Jacobian, as at a corner of the front the residual of every multiplier is zero.
  const Real noise =
      8 * BOOST_MATH_NAMESPACE::tools::root_epsilon<Real>() * (std::max)(largest, sqrt(size));
  const Real limit = (std::max)(2 * s.weak, 200 * tolerance) + noise;
  const Real sum1 = U[3] + U[4] + U[5], sum2 = U[6] + U[7] + U[8];
  Real along[3], length = 0;
  for (size_t i = 0; i < 3; ++i) {
    along[i] = U[3 + i] * sum2 - U[6 + i] * sum1;
    length += along[i] * along[i];
  }
  length = sqrt(length);
  s.faces.clear();
  s.wide = false;
  s.slack = s.weak;
  for (int ray = 0; ray < 8; ++ray) {
    // The first two along the segment, one each way; the others towards and away from the corners of the simplex.
    Real d[3], reach = std::numeric_limits<Real>::infinity();
    for (size_t i = 0; i < 3; ++i) {
      const Real corner = (size_t(ray - 2) / 2 == i ? Real(2) : Real(-1)) / 3;
      d[i] = (ray % 2 ? Real(-1) : Real(1)) * (ray < 2 ? (length > 0 ? along[i] / length : Real(0)) : corner);
      reach = d[i] < 0 ? (std::min)(reach, s.lam[i] / -d[i]) : reach;
    }
    if (!(reach > 0) || !(reach < std::numeric_limits<Real>::infinity())) {
      continue;
    }
    Vec l(3);
    auto at = [&](Real t) {
      for (size_t i = 0; i < 3; ++i) {
        l[i] = (std::max)(s.lam[i] + t * d[i], Real(0));
      }
      return residual(l);
    };
    const Real t = pareto_polish_bisect(at, Real(0), reach, limit, 40), r = at(t);
    Real across = 0;
    for (size_t i = 0; i < 3 && length > 0; ++i) {
      across += d[i] * along[i] / length;
    }
    if (ray < 2 && t > Real(1) / 1000) {
      s.faces.insert(s.faces.end(), l.begin(), l.end());
      s.slack = (std::max)(s.slack, r);
    } else if (ray >= 2 && t * sqrt((std::max)(Real(2) / 3 - across * across, Real(0))) > Real(1) / 20) {
      s.wide = true;
    }
  }
  s.stationary = !s.inside && s.slack <= 200 * tolerance + noise;
}

// Three objectives: the triangles of a front, and what the searches into them have found. `tried` has the
// triangles by the numbers of their corners: the attempts made to divide each, or -1 for one which is presumed
// outside the front, the corners of which are in `beyond`, in the objectives. `sagitta` has the edges where the
// triangulation ends: the estimate of the error of each, once a point was sought between its ends, and `sought`
// how often. The chart: the direction r which the points are projected along, and two unit vectors which span the
// plane orthogonal to it.
template <typename Real> struct pareto_polish_triangulation {
  using Vec = std::vector<Real>;
  std::map<std::array<size_t, 3>, int> tried;
  std::map<std::pair<size_t, size_t>, Real> sagitta;
  std::map<std::pair<size_t, size_t>, int> sought;
  std::vector<Vec> beyond;
  Vec chart;
  template <class State, class Edge>
  static std::array<size_t, 3> corners(std::vector<State> const &states, Edge const &e) {
    std::array<size_t, 3> key{{states[e.i].id, states[e.j].id, states[e.k].id}};
    std::sort(key.begin(), key.end());
    return key;
  }
  // Three objectives: the cells of the points `order`, see the documentation. The points are mutually
  // nondominated, so no two of them differ by a multiple of a direction r with positive entries: they are a graph
  // over the plane orthogonal to r, and are triangulated there. r is the mean of their unit normals.
  template <class State, class Edge>
  void cells(std::vector<State> const &states, std::vector<size_t> const &order, Vec const &range, Real tol,
             Real longest, Real max_edge_length, std::vector<Edge> &edges) {
    using std::sqrt;
    using edge = Edge;
    using cell_class = decltype(Edge().kind);
    chart.assign(9, Real(0));
    for (size_t e : order) {
      const Real norm = states[e].lam.size() == 3 ? pareto_polish_norm(states[e].lam) : Real(0);
      for (size_t i = 0; i < 3 && norm > 0; ++i) {
        chart[i] += states[e].lam[i] / norm;
      }
    }
    const Real most = (std::max)(chart[0], (std::max)(chart[1], chart[2]));
    Real norm = 0;
    size_t least = 0;
    for (size_t i = 0; i < 3; ++i) {
      chart[i] = most > 0 ? (std::max)(chart[i], most / 10) : Real(1);
      norm += chart[i] * chart[i];
      least = chart[i] < chart[least] ? i : least;
    }
    Real norm1 = 0;
    for (size_t i = 0; i < 3; ++i) {
      chart[i] /= sqrt(norm);
    }
    for (size_t i = 0; i < 3; ++i) {
      chart[3 + i] = (i == least ? Real(1) : Real(0)) - chart[least] * chart[i];
      norm1 += chart[3 + i] * chart[3 + i];
    }
    for (size_t i = 0; i < 3; ++i) {
      chart[3 + i] /= sqrt(norm1);
    }
    for (size_t i = 0; i < 3; ++i) {
      const size_t j = (i + 1) % 3, k = (i + 2) % 3;
      chart[6 + i] = chart[j] * chart[3 + k] - chart[k] * chart[3 + j];
    }
    const Vec r(chart.begin(), chart.begin() + 3);
    pareto_polish_delaunay<Real> mesh;
    std::vector<Vec> Y(order.size(), Vec(3));
    for (size_t q = 0; q < order.size(); ++q) {
      Real u[2] = {0, 0};
      for (size_t i = 0; i < 3; ++i) {
        Y[q][i] = states[order[q]].f[i] / range[i];
        u[0] += Y[q][i] * chart[3 + i];
        u[1] += Y[q][i] * chart[6 + i];
      }
      mesh.insert(u[0], u[1]);
    }
    const auto triangles = mesh.triangles();
    std::map<std::pair<size_t, size_t>, size_t> owner; // the triangle of a directed edge
    std::vector<char> estimated;
    for (auto const &T : triangles) {
      edge e;
      e.i = order[T[0]];
      e.j = order[T[1]];
      e.k = order[T[2]];
      State const *corner[3] = {&states[e.i], &states[e.j], &states[e.k]};
      Vec P(9), N(9), allowance(3, Real(0)), centre(3, Real(0));
      std::vector<char> skip(3, 0);
      bool tangents = true, tight = true, estimate = false;
      size_t at = 0;
      for (size_t c = 0; c < 3; ++c) {
        for (size_t i = 0; i < 3; ++i) {
          P[3 * c + i] = Y[T[c]][i];
          N[3 * c + i] = r[i];
          centre[i] += P[3 * c + i] / 3;
        }
        tangents = tangents && corner[c]->lam.size() == 3;
        const Real length = pareto_polish_distance(Y[T[c]], Y[T[(c + 1) % 3]]);
        if (length > e.length) {
          e.length = length;
          at = c;
        }
      }
      for (size_t c = 0; c < 3 && tangents; ++c) {
        // Of the multipliers of the corner, those whose tangent plane is lowest at the centre of the cell: the
        // plane which faces it.
        State const &a = *corner[c];
        Real best = std::numeric_limits<Real>::infinity(), length = 1, along = 0;
        for (size_t v = 0; v <= a.faces.size() / 3; ++v) {
          Real const *l = v == 0 ? a.lam.data() : &a.faces[3 * (v - 1)];
          const Real size = sqrt(l[0] * l[0] + l[1] * l[1] + l[2] * l[2]);
          Real h = 0;
          for (size_t i = 0; i < 3; ++i) {
            h += l[i] * (centre[i] - P[3 * c + i]) / size;
          }
          if (h < best) {
            best = h;
            length = size;
            std::copy(l, l + 3, &N[3 * c]);
          }
        }
        for (size_t i = 0; i < 3; ++i) {
          along += N[3 * c + i] * r[i] / length;
        }
        const Real widest = (std::max)(pareto_polish_distance(a.x, corner[(c + 1) % 3]->x),
                                       pareto_polish_distance(a.x, corner[(c + 2) % 3]->x));
        allowance[c] = a.slack * widest / length;
        // The only multipliers of the corner are those of one objective alone: it is a stationary point of that
        // objective, which the minimizer of it is, and so is a point which is only weakly Pareto optimal, off the
        // front. Its tangent plane is not used, and the cell is no bound.
        skip[c] = a.faces.empty() && (std::max)(a.lam[0], (std::max)(a.lam[1], a.lam[2])) >= 1 - Real(1) / 1000000;
        tight = tight && a.stationary && !a.wide && !skip[c] && along >= Real(1) / 20;
        estimate = estimate || skip[c];
      }
      const auto cell = pareto_polish_simplex_error(P, N, allowance, r, 3, skip);
      for (size_t c = 0; c < 3 && tangents; ++c) {
        e.normals.insert(e.normals.end(), 3, Real(0));
        std::copy(&N[3 * c], &N[3 * c] + (skip[c] ? 0 : 3), e.normals.end() - 3);
      }
      e.kind = tangents ? static_cast<cell_class>(cell.kind) : cell_class::mixed;
      e.flat = tangents && cell.flat;
      e.bound = tangents && tight && cell.kind != 2;
      e.error = cell.error;
      e.box = cell.dominance;
      e.flat = e.flat && !estimate;
      estimated.push_back(estimate && cell.kind != 2); // the error of the envelope stays, as an estimate
      // A sliver, lower than a twentieth of its longest edge, is divided in the middle of that edge.
      Real cross = 0;
      for (size_t i = 0; i < 3; ++i) {
        const size_t j = (i + 1) % 3, k = (i + 2) % 3;
        const Real v = (P[3 + j] - P[j]) * (P[6 + k] - P[k]) - (P[3 + k] - P[k]) * (P[6 + j] - P[j]);
        cross += v * v;
        e.b[i] = cell.b[i];
      }
      if ((std::max)(e.b[0], (std::max)(e.b[1], e.b[2])) > Real(9) / 10) {
        // No point to divide it at near a corner: where a corner gives no tangent plane the envelope may be
        // farthest from the cell there.
        e.b[0] = e.b[1] = e.b[2] = Real(1) / 3;
      }
      if (!(sqrt(cross) >= e.length * e.length / 20)) {
        e.b[at] = e.b[(at + 1) % 3] = Real(1) / 2;
        e.b[(at + 2) % 3] = 0;
      }
      for (size_t c = 0; c < 3; ++c) {
        owner[std::make_pair(T[c], T[(c + 1) % 3])] = edges.size();
      }
      edges.push_back(e);
    }
    // Where a convex cell shares an edge with a concave one the curvature changes sign in one of the two:
    // neither is a bound.
    std::vector<char> clash(edges.size(), 0);
    for (auto const &o : owner) {
      auto const twin = owner.find(std::make_pair(o.first.second, o.first.first));
      if (twin != owner.end()) {
        edge const &a = edges[o.second], &b = edges[twin->second];
        clash[o.second] = clash[o.second] || (a.kind != cell_class::mixed && b.kind != cell_class::mixed && !a.flat &&
                                              !b.flat && a.kind != b.kind);
      }
    }
    for (size_t t = 0; t < clash.size(); ++t) {
      edge &e = edges[t];
      e.bound = e.bound && !clash[t];
      e.error = e.bound || (estimated[t] && !clash[t]) ? e.error : e.box;
      // Within a triangle which is presumed outside the front: outside as well.
      Real u[2] = {0, 0};
      for (size_t i = 0; i < 3; ++i) {
        const Real y = (states[e.i].f[i] + states[e.j].f[i] + states[e.k].f[i]) / (3 * range[i]);
        u[0] += y * chart[3 + i];
        u[1] += y * chart[6 + i];
      }
      bool within = false;
      for (size_t q = 0; q < beyond.size() && !within; ++q) {
        Real w[3][2];
        for (size_t c = 0; c < 3; ++c) {
          w[c][0] = w[c][1] = 0;
          for (size_t i = 0; i < 3; ++i) {
            w[c][0] += beyond[q][3 * c + i] / range[i] * chart[3 + i];
            w[c][1] += beyond[q][3 * c + i] / range[i] * chart[6 + i];
          }
        }
        auto side = [&](size_t a, size_t b) {
          return (w[b][0] - w[a][0]) * (u[1] - w[a][1]) - (w[b][1] - w[a][1]) * (u[0] - w[a][0]);
        };
        const Real s0 = side(0, 1), s1 = side(1, 2), s2 = side(2, 0);
        within = (s0 > 0 && s1 > 0 && s2 > 0) || (s0 < 0 && s1 < 0 && s2 < 0);
      }
      if (within) {
        tried[corners(states, e)] = -1;
      }
      auto const known = tried.find(corners(states, e));
      e.stage = known == tried.end() ? 0 : known->second;
      e.resolved = e.error <= tol && e.length <= longest;
      if (e.stage < 0) {
        e.kind = cell_class::boundary;
        e.bound = false;
        e.resolved = true;
      } else if (!e.resolved && e.stage >= 2) {
        // Within the tolerance and only longer than the quarter at which a triangle is divided for density: it
        // stays as it is, unless it is longer than max_edge_length allows.
        e.resolved = e.error <= tol && e.length <= max_edge_length;
        e.kind = e.resolved ? e.kind : cell_class::unresolved;
      }
    }
    // The edges where the triangles end.
    for (auto const &o : owner) {
      auto const twin = owner.find(std::make_pair(o.first.second, o.first.first));
      if (edges[o.second].kind == cell_class::boundary ||
          (twin != owner.end() && edges[twin->second].kind != cell_class::boundary)) {
        continue;
      }
      edge e;
      e.i = order[o.first.first];
      e.j = order[o.first.second];
      e.length = pareto_polish_distance(Y[o.first.first], Y[o.first.second]);
      const auto key = (std::minmax)(states[e.i].id, states[e.j].id);
      auto const known = sagitta.find(key);
      auto const before = sought.find(key);
      e.stage = before == sought.end() ? 0 : before->second;
      e.error = known != sagitta.end() ? known->second : e.length / 2;
      e.resolved = e.error <= tol;
      e.kind = e.resolved || e.stage < 2 ? cell_class::boundary : cell_class::unresolved;
      // The target of the second search: a quarter of the length of the edge beyond its middle, away from its
      // triangle, in the plane of the chart.
      const edge inner = edges[o.second];
      const size_t third =
          inner.i != e.i && inner.i != e.j ? inner.i : (inner.j != e.i && inner.j != e.j ? inner.j : inner.k);
      Real along[2] = {0, 0}, in[2] = {0, 0};
      e.apex.assign(3, Real(0));
      for (size_t i = 0; i < 3; ++i) {
        const Real a = states[e.i].f[i] / range[i], b = states[e.j].f[i] / range[i], c = states[third].f[i] / range[i];
        e.apex[i] = (a + b) / 2;
        for (size_t w = 0; w < 2; ++w) {
          along[w] += (b - a) * chart[3 + 3 * w + i];
          in[w] += (c - (a + b) / 2) * chart[3 + 3 * w + i];
        }
      }
      const Real length2 = along[0] * along[0] + along[1] * along[1];
      const Real dot = length2 > 0 ? (in[0] * along[0] + in[1] * along[1]) / length2 : Real(0);
      Real out[2] = {dot * along[0] - in[0], dot * along[1] - in[1]};
      const Real size = sqrt(out[0] * out[0] + out[1] * out[1]);
      for (size_t i = 0; i < 3 && size > 0; ++i) {
        e.apex[i] += sqrt(length2) / 4 * (out[0] * chart[3 + i] + out[1] * chart[6 + i]) / size;
      }
      edges.push_back(e);
    }
    // A point which the triangulation left out (its projection repeats another, to rounding) is in no triangle:
    // the edge to the nearest point which is in one is returned as unresolved, so that the point is not passed
    // over.
    for (size_t q = 0; q < order.size() && !triangles.empty(); ++q) {
      size_t nearest = order.size();
      for (size_t p = 0; p < order.size() && !mesh.inserted(q); ++p) {
        if (mesh.inserted(p) &&
            (nearest == order.size() || pareto_polish_distance(Y[q], Y[p]) < pareto_polish_distance(Y[q], Y[nearest]))) {
          nearest = p;
        }
      }
      if (nearest < order.size()) {
        edge e;
        e.i = order[q];
        e.j = order[nearest];
        e.length = e.error = pareto_polish_distance(Y[q], Y[nearest]);
        e.stage = 2;
        e.kind = cell_class::unresolved;
        edges.push_back(e);
      }
    }
  }
};

}} BOOST_MATH_NAMESPACE_END // namespace boost::math::optimization::detail
#endif
