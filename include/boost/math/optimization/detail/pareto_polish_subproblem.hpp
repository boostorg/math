/*
 * Copyright Nick Thompson, 2026
 * Use, modification and distribution are subject to the
 * Boost Software License, Version 1.0. (See accompanying file
 * LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef BOOST_MATH_OPTIMIZATION_DETAIL_PARETO_POLISH_SUBPROBLEM_HPP
#define BOOST_MATH_OPTIMIZATION_DETAIL_PARETO_POLISH_SUBPROBLEM_HPP
#include <boost/math/tools/config.hpp>
#ifndef BOOST_MATH_BUILD_MODULE
#include <algorithm>
#include <cmath>
#include <limits>
#include <utility>
#include <vector>
#endif

BOOST_MATH_NAMESPACE_BEGIN namespace optimization { namespace detail {

template <typename Real> Real pareto_polish_distance(std::vector<Real> const &a, std::vector<Real> const &b) {
  using std::sqrt;
  Real s = 0;
  for (size_t k = 0; k < a.size(); ++k) {
    s += (a[k] - b[k]) * (a[k] - b[k]);
  }
  return sqrt(s);
}

template <typename Real> Real pareto_polish_norm(std::vector<Real> const &v) {
  using std::sqrt;
  Real s = 0;
  for (auto const &a : v) {
    s += a * a;
  }
  return sqrt(s);
}

template <typename Real> bool pareto_polish_finite(std::vector<Real> const &v) {
  using std::isfinite;
  for (auto const &a : v) {
    if (!isfinite(a)) {
      return false;
    }
  }
  return true;
}

// The norm of the step -w(k) of each variable k cut to lo <= d <= hi: the weak residual of the multipliers lam
// with w = J'lam.
template <class W, typename Real>
Real pareto_polish_residual(W w, std::vector<Real> const &lo, std::vector<Real> const &hi) {
  using std::sqrt;
  Real r = 0;
  for (size_t k = 0; k < lo.size(); ++k) {
    const Real u = (std::min)((std::max)(-w(k), lo[k]), hi[k]);
    r += u * u;
  }
  return sqrt(r);
}

// The end of the segment from `good` to `bad` on which residual() stays within the limit, by bisection.
template <class F, typename Real> Real pareto_polish_bisect(F residual, Real good, Real bad, Real limit, int steps) {
  for (int it = 0; it < steps && residual(bad) > limit; ++it) {
    const Real mid = (bad + good) / 2;
    (residual(mid) > limit ? bad : good) = mid;
  }
  return residual(bad) > limit ? good : bad;
}

// Gaussian elimination with partial pivoting of the K by K + 1 augmented matrix a. False if singular.
template <typename Real> bool pareto_polish_solve(std::vector<Real> &a, size_t K, std::vector<Real> &x) {
  using std::abs;
  const size_t W = K + 1;
  Real scale = 0;
  for (size_t r = 0; r < K; ++r) {
    for (size_t c = 0; c < K; ++c) {
      scale = (std::max)(scale, abs(a[r * W + c]));
    }
  }
  for (size_t c = 0; c < K; ++c) {
    size_t p = c;
    for (size_t r = c + 1; r < K; ++r) {
      if (abs(a[r * W + c]) > abs(a[p * W + c])) {
        p = r;
      }
    }
    if (!(abs(a[p * W + c]) > 64 * std::numeric_limits<Real>::epsilon() * scale)) {
      return false;
    }
    for (size_t b = 0; b < W; ++b) {
      std::swap(a[p * W + b], a[c * W + b]);
    }
    for (size_t r = c + 1; r < K; ++r) {
      const Real w = a[r * W + c] / a[c * W + c];
      for (size_t b = c; b < W; ++b) {
        a[r * W + b] -= w * a[c * W + b];
      }
    }
  }
  x.resize(K);
  for (size_t r = K; r-- > 0;) {
    Real v = a[r * W + K];
    for (size_t b = r + 1; b < K; ++b) {
      v -= a[r * W + b] * x[b];
    }
    x[r] = v / a[r * W + r];
  }
  return true;
}

// Maximizes ct.lam - lam'H lam / 2 over the simplex, or over lam >= 1 (cone), H positive semidefinite, by
// enumeration of the faces: the maximizer is a stationary point of the face it is in the relative interior of.
// A face with a singular system is skipped, as then a maximizer lies on a smaller one.
template <typename Real>
void pareto_polish_face_qp(std::vector<Real> const &H, std::vector<Real> const &ct, bool cone,
                           std::vector<Real> &lam) {
  const size_t m = ct.size();
  auto value = [&](std::vector<Real> const &l) {
    Real v = 0;
    for (size_t a = 0; a < m; ++a) {
      v += ct[a] * l[a];
      for (size_t b = 0; b < m; ++b) {
        v -= l[a] * H[a * m + b] * l[b] / 2;
      }
    }
    return v;
  };
  Real best = -std::numeric_limits<Real>::infinity();
  std::vector<Real> trial(m), a, sol;
  std::vector<size_t> S;
  for (size_t mask = cone ? 0 : 1; mask < (size_t(1) << m); ++mask) {
    S.clear();
    for (size_t i = 0; i < m; ++i) {
      if (mask >> i & 1) {
        S.push_back(i);
      }
    }
    const size_t k = S.size();
    std::fill(trial.begin(), trial.end(), Real(cone ? 1 : 0));
    if (cone && k > 0) {
      // H_SS mu = ct_S - (H 1)_S, lam = 1 + mu
      a.assign(k * (k + 1), Real(0));
      for (size_t r = 0; r < k; ++r) {
        Real rhs = ct[S[r]];
        for (size_t b = 0; b < m; ++b) {
          rhs -= H[S[r] * m + b];
        }
        for (size_t c = 0; c < k; ++c) {
          a[r * (k + 1) + c] = H[S[r] * m + S[c]];
        }
        a[r * (k + 1) + k] = rhs;
      }
      if (!pareto_polish_solve(a, k, sol)) {
        continue;
      }
      bool ok = true;
      for (size_t r = 0; r < k; ++r) {
        ok = ok && sol[r] >= 0;
        trial[S[r]] = 1 + sol[r];
      }
      if (!ok) {
        continue;
      }
    } else if (!cone && k == 1) {
      trial[S[0]] = 1;
    } else if (!cone) {
      // H_SS lam + nu = ct_S, sum lam = 1; nu is scaled to keep the pivots comparable.
      const size_t K = k + 1;
      Real scale = 0;
      for (size_t r = 0; r < k; ++r) {
        scale = (std::max)(scale, H[S[r] * m + S[r]]);
      }
      if (!(scale > 0)) {
        continue;
      }
      a.assign(K * (K + 1), Real(0));
      for (size_t r = 0; r < k; ++r) {
        for (size_t c = 0; c < k; ++c) {
          a[r * (K + 1) + c] = H[S[r] * m + S[c]];
        }
        a[r * (K + 1) + k] = a[k * (K + 1) + r] = scale;
        a[r * (K + 1) + K] = ct[S[r]];
      }
      a[k * (K + 1) + K] = scale;
      if (!pareto_polish_solve(a, K, sol)) {
        continue;
      }
      bool ok = true;
      Real sum = 0;
      for (size_t r = 0; r < k; ++r) {
        ok = ok && sol[r] >= 0;
        sum += sol[r];
      }
      if (!ok || !(sum > 0)) {
        continue;
      }
      for (size_t r = 0; r < k; ++r) {
        trial[S[r]] = sol[r] / sum;
      }
    }
    const Real v = value(trial);
    if (v > best) {
      best = v;
      lam = trial;
    }
  }
}

// d = argmin over lo <= d <= hi (lo <= 0 <= hi, or lo = hi) of M(d) + |d|^2 / (2 tau), with g_i the rows of the
// m by n matrix J, and
//   M(d) = max_i (c_i + g_i.d),                             lam in the simplex, or
//   M(d) = sum_i g_i.d subject to g_i.d <= 0 for every i,   lam >= 1 (cone; c is not used).
// In both d = clip(-tau J'lam, lo, hi) at the multipliers lam, which maximize the dual
//   q(lam) = c.lam + sum_k min over lo_k <= u <= hi_k of (w_k u + u^2 / (2 tau)),   w = J'lam.
// q is concave, differentiable and piecewise quadratic: one quadratic for each set of clipped variables. Newton's
// method with an exact line search: the quadratic at lam (with a small proximal term, so that it is bounded) is
// maximized by pareto_polish_face_qp, and q is maximized on the ray through that point, where its derivative is
// piecewise linear. It stops where the segment has no ascent, which is
// where lam is optimal.
BOOST_MATH_TEST_EXPORT template <typename Real>
void pareto_polish_subproblem(std::vector<Real> const &J, size_t m, std::vector<Real> const &c,
                              std::vector<Real> const &lo, std::vector<Real> const &hi, Real tau, bool cone,
                              std::vector<Real> &lam, std::vector<Real> &d) {
  using std::abs;
  using std::sqrt;
  const size_t n = lo.size();
  const Real eps = std::numeric_limits<Real>::epsilon();
  lam.assign(m, Real(cone ? 1 : Real(1) / static_cast<Real>(m)));
  d.assign(n, Real(0));
  std::vector<Real> w(n), omega(n), H(m * m), ct(m), target(m), v(m);
  std::vector<std::pair<Real, size_t>> breaks;
  Real trace = 0;
  for (auto const &g : J) {
    trace += tau * g * g;
  }
  auto clip = [&]() {
    for (size_t k = 0; k < n; ++k) {
      w[k] = 0;
      for (size_t i = 0; i < m; ++i) {
        w[k] += lam[i] * J[i * n + k];
      }
      d[k] = (std::min)((std::max)(-tau * w[k], lo[k]), hi[k]);
    }
  };
  // The derivative of q along v at lam + t v, and the magnitude of its terms.
  auto slope = [&](Real t, Real &size) {
    Real s = 0;
    size = 0;
    for (size_t i = 0; i < m && !cone; ++i) {
      s += c[i] * v[i];
      size += abs(c[i] * v[i]);
    }
    for (size_t k = 0; k < n; ++k) {
      const Real u = (std::min)((std::max)(-tau * (w[k] + t * omega[k]), lo[k]), hi[k]);
      s += omega[k] * u;
      size += abs(omega[k] * u);
    }
    return s;
  };
  // The exact maximizer of q on the ray from lam along v, as far as the multipliers allow; lam moves there.
  // False if v is not a direction of ascent, to rounding.
  auto search = [&]() {
    for (size_t k = 0; k < n; ++k) {
      omega[k] = 0;
      for (size_t i = 0; i < m; ++i) {
        omega[k] += v[i] * J[i * n + k];
      }
    }
    // The gain along v to first order, against the size of q.
    Real size, q_size = 0;
    for (size_t k = 0; k < n; ++k) {
      q_size += abs(w[k] * d[k]) + d[k] * d[k] / (2 * tau);
    }
    for (size_t i = 0; i < m && !cone; ++i) {
      q_size += abs(c[i] * lam[i]);
    }
    Real s0 = slope(0, size);
    if (!(s0 > 16 * eps * (size + q_size))) {
      return false;
    }
    // The breakpoints of the derivative, in order; it is linear between them and decreasing.
    Real end = (std::numeric_limits<Real>::max)();
    for (size_t i = 0; i < m; ++i) {
      if (v[i] < 0) {
        end = (std::min)(end, (lam[i] - (cone ? 1 : 0)) / -v[i]);
      }
    }
    const bool bounded = end < (std::numeric_limits<Real>::max)();
    breaks.clear();
    for (size_t k = 0; k < n; ++k) {
      for (Real bound : {lo[k], hi[k]}) {
        const Real t = (-bound / tau - w[k]) / omega[k];
        if (lo[k] < hi[k] && t > 0 && t < end) {
          breaks.emplace_back(t, k);
        }
      }
    }
    std::sort(breaks.begin(), breaks.end());
    if (bounded) {
      breaks.emplace_back(end, n);
    }
    Real t0 = 0, step = -1;
    for (auto const &b : breaks) {
      const Real s1 = slope(b.first, size);
      if (s1 <= 0) {
        step = t0 + (b.first - t0) * s0 / (s0 - s1);
        break;
      }
      t0 = b.first;
      s0 = s1;
    }
    if (step < 0) {
      // Past the last breakpoint: the boundary, or on a ray without end the zero of the last linear piece.
      step = t0;
      if (!bounded) {
        const Real s1 = slope(2 * t0 + 1, size);
        if (s1 < s0) {
          step = t0 + (t0 + 1) * s0 / (s0 - s1);
        }
      }
    }
    for (size_t i = 0; i < m; ++i) {
      lam[i] += step * v[i];
      lam[i] = cone ? (std::max)(lam[i], Real(1)) : (std::max)(lam[i], Real(0));
    }
    clip();
    return step > 0;
  };
  clip();
  for (size_t iteration = 0; iteration < 100 + n; ++iteration) {
    const Real rho = 1024 * eps * trace + (std::numeric_limits<Real>::min)();
    std::fill(H.begin(), H.end(), Real(0));
    for (size_t i = 0; i < m; ++i) {
      ct[i] = (cone ? Real(0) : c[i]) + rho * lam[i];
      H[i * m + i] = rho;
    }
    for (size_t k = 0; k < n; ++k) {
      const bool clipped = lo[k] == hi[k] || -tau * w[k] < lo[k] || -tau * w[k] > hi[k];
      for (size_t a = 0; a < m; ++a) {
        if (clipped) {
          ct[a] += J[a * n + k] * d[k];
        } else {
          for (size_t b = 0; b < m; ++b) {
            H[a * m + b] += tau * J[a * n + k] * J[b * n + k];
          }
        }
      }
    }
    target = lam;
    pareto_polish_face_qp(H, ct, cone, target);
    Real moved = 0, largest = 1;
    for (size_t i = 0; i < m; ++i) {
      v[i] = target[i] - lam[i];
      moved = (std::max)(moved, abs(v[i]));
      largest = (std::max)(largest, lam[i]);
    }
    if (!(moved > 64 * eps * largest) || !search()) {
      break;
    }
  }
  // The directions of Newton's method are long where the quadratic is nearly flat, and the derivative along them
  // is then lost to cancellation. The coordinates do not have this problem: they finish the multipliers of the
  // cone, so that g_i.d <= 0 holds to rounding.
  for (size_t sweep = 0; cone && sweep < 50; ++sweep) {
    bool any = false;
    for (size_t i = 0; i < m; ++i) {
      for (Real sign : {Real(1), Real(-1)}) {
        if (sign > 0 || lam[i] > 1) {
          std::fill(v.begin(), v.end(), Real(0));
          v[i] = sign * lam[i];
          const Real before = lam[i];
          any = (search() && abs(lam[i] - before) > 64 * eps * before) || any;
        }
      }
    }
    if (!any) {
      break;
    }
  }
}

}} BOOST_MATH_NAMESPACE_END // namespace boost::math::optimization::detail
#endif
