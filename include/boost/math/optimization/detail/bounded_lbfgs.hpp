/*
 * Copyright Nick Thompson, 2026
 * Use, modification and distribution are subject to the
 * Boost Software License, Version 1.0. (See accompanying file
 * LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef BOOST_MATH_OPTIMIZATION_DETAIL_BOUNDED_LBFGS_HPP
#define BOOST_MATH_OPTIMIZATION_DETAIL_BOUNDED_LBFGS_HPP

#include <boost/math/tools/config.hpp>
#ifndef BOOST_MATH_BUILD_MODULE
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <utility>
#include <vector>
#endif

BOOST_MATH_NAMESPACE_BEGIN namespace optimization::detail {

BOOST_MATH_TEST_EXPORT struct bounded_lbfgs_parameters {
  size_t max_iterations = 200;
  // The number of correction pairs kept.
  size_t memory = 10;
  // Stop when the infinity norm of the projected gradient P(x - g) - x is at most this,
  double gradient_tolerance = 1e-5;
  // or when an iteration decreases f by at most this times max(|f|, 1) (the defaults are those of SciPy's L-BFGS-B).
  double function_tolerance = 2.2e-9;
  bool maximize = false;
  // The line search of L-BFGS-B in place of backtracking: that of More and Thuente, to the strong Wolfe conditions,
  // along the step cut short at the first bound it meets. While no bound is met, the iterates are then those of
  // L-BFGS-B (Zhu, Byrd, Lu and Nocedal, Algorithm 778, ACM Trans. Math. Software 23 (1997)), which matters where the
  // function has several local minima and the one that SciPy finds is wanted.
  bool strong_wolfe = false;
};

BOOST_MATH_TEST_EXPORT struct bounded_lbfgs_result {
  double value;              // f at the returned x; not finite if f is not finite at the starting point
  double projected_gradient; // the infinity norm of P(x - g) - x there (of P(x + g) - x when maximizing)
  size_t iterations;
  size_t evaluations;
};

// The line search of More, J. J., Thuente, D. J. (1994), Line search algorithms with guaranteed sufficient decrease,
// ACM Trans. Math. Software 20, 286-307, as dcsrch and dcstep of MINPACK-2 have it, with the constants of L-BFGS-B:
// finds a step with f(stp) <= f(0) + 1e-3 stp f'(0) and |f'(stp)| <= 0.9 |f'(0)| in [0, stpmax]. Construct it with
// f(0), f'(0) < 0 and the first step to try; then, while next(f(stp), f'(stp)) is true, evaluate at the new stp.
// When it is false, stp is the step to take: it has been evaluated, and satisfies the conditions unless they cannot
// be satisfied (the step is stpmax, or rounding errors prevent progress).
struct more_thuente {
  double stp;
  more_thuente(double f0, double g0, double first, double stpmax)
      : stp(first), finit(f0), ginit(g0), gtest(1e-3 * g0), stpmax(stpmax), width(stpmax), width1(2 * stpmax), fx_(f0),
        gx_(g0), fy_(f0), gy_(g0), stmax(5 * first) {}

  bool next(double f, double g) {
    using std::abs;
    double const ftest = finit + stp * gtest;
    if (stage == 1 && f <= ftest && g >= 0) {
      stage = 2;
    }
    bool const stuck = brackt && (stp <= stmin || stp >= stmax || stmax - stmin <= 0.1 * stmax);
    bool const at_max = stp == stpmax && f <= ftest && g <= gtest;
    bool const at_min = stp == 0 && (f > ftest || g >= gtest);
    if (stuck || at_max || at_min || (f <= ftest && abs(g) <= 0.9 * -ginit)) {
      return false;
    }
    // While the sufficient decrease has not been reached with a nonnegative slope, the interval is updated for the
    // function less the line of sufficient decrease.
    if (stage == 1 && f <= fx_ && f > ftest) {
      double fm = f - stp * gtest, fxm = fx_ - stx * gtest, fym = fy_ - sty * gtest;
      double gm = g - gtest, gxm = gx_ - gtest, gym = gy_ - gtest;
      step(fxm, gxm, fym, gym, fm, gm);
      fx_ = fxm + stx * gtest;
      fy_ = fym + sty * gtest;
      gx_ = gxm + gtest;
      gy_ = gym + gtest;
    }
    else {
      step(fx_, gx_, fy_, gy_, f, g);
    }
    if (brackt) {
      // A bisection if the interval has not shrunk enough:
      if (abs(sty - stx) >= 0.66 * width1) {
        stp = stx + 0.5 * (sty - stx);
      }
      width1 = width;
      width = abs(sty - stx);
      stmin = (std::min)(stx, sty);
      stmax = (std::max)(stx, sty);
    }
    else {
      stmin = stp + 1.1 * (stp - stx);
      stmax = stp + 4 * (stp - stx);
    }
    stp = (std::min)((std::max)(stp, 0.0), stpmax);
    if (brackt && (stp <= stmin || stp >= stmax || stmax - stmin <= 0.1 * stmax)) {
      stp = stx;
    }
    return true;
  }

private:
  // dcstep: the next trial step from the best step so far (stx), the other end of the interval (sty) and the trial
  // (stp), by safeguarded cubic and quadratic interpolation; updates the interval.
  void step(double &fx, double &dx, double &fy, double &dy, double fp, double dp) {
    using std::abs;
    using std::sqrt;
    double const sgnd = dp * (dx / abs(dx));
    double stpf;
    auto cubic = [](double theta, double a, double b) {
      double const s = (std::max)({abs(theta), abs(a), abs(b)});
      return s * sqrt((std::max)(0.0, (theta / s) * (theta / s) - (a / s) * (b / s)));
    };
    if (fp > fx) {
      // A higher value: the minimum is bracketed.
      double const theta = 3 * (fx - fp) / (stp - stx) + dx + dp;
      double gamma = cubic(theta, dx, dp);
      gamma = stp < stx ? -gamma : gamma;
      double const r = ((gamma - dx) + theta) / (((gamma - dx) + gamma) + dp);
      double const stpc = stx + r * (stp - stx);
      double const stpq = stx + ((dx / ((fx - fp) / (stp - stx) + dx)) / 2) * (stp - stx);
      stpf = abs(stpc - stx) < abs(stpq - stx) ? stpc : stpc + (stpq - stpc) / 2;
      brackt = true;
    }
    else if (sgnd < 0) {
      // A lower value and slopes of opposite sign: bracketed too.
      double const theta = 3 * (fx - fp) / (stp - stx) + dx + dp;
      double gamma = cubic(theta, dx, dp);
      gamma = stp > stx ? -gamma : gamma;
      double const r = ((gamma - dp) + theta) / (((gamma - dp) + gamma) + dx);
      double const stpc = stp + r * (stx - stp);
      double const stpq = stp + (dp / (dp - dx)) * (stx - stp);
      stpf = abs(stpc - stp) > abs(stpq - stp) ? stpc : stpq;
      brackt = true;
    }
    else if (abs(dp) < abs(dx)) {
      // A lower value, slopes of the same sign, and the slope decreases in magnitude.
      double const theta = 3 * (fx - fp) / (stp - stx) + dx + dp;
      double gamma = cubic(theta, dx, dp);
      gamma = stp > stx ? -gamma : gamma;
      double const r = ((gamma - dp) + theta) / ((gamma + (dx - dp)) + gamma);
      double const stpc = r < 0 && gamma != 0 ? stp + r * (stx - stp) : (stp > stx ? stmax : stmin);
      double const stpq = stp + (dp / (dp - dx)) * (stx - stp);
      if (brackt) {
        stpf = abs(stpc - stp) < abs(stpq - stp) ? stpc : stpq;
        stpf = stp > stx ? (std::min)(stp + 0.66 * (sty - stp), stpf) : (std::max)(stp + 0.66 * (sty - stp), stpf);
      }
      else {
        stpf = abs(stpc - stp) > abs(stpq - stp) ? stpc : stpq;
        stpf = (std::max)(stmin, (std::min)(stmax, stpf));
      }
    }
    else if (brackt) {
      // A lower value, slopes of the same sign, and the slope does not decrease.
      double const theta = 3 * (fp - fy) / (sty - stp) + dy + dp;
      double gamma = cubic(theta, dy, dp);
      gamma = stp > sty ? -gamma : gamma;
      double const r = ((gamma - dp) + theta) / (((gamma - dp) + gamma) + dy);
      stpf = stp + r * (sty - stp);
    }
    else {
      stpf = stp > stx ? stmax : stmin;
    }
    if (fp > fx) {
      sty = stp;
      fy = fp;
      dy = dp;
    }
    else {
      if (sgnd < 0) {
        sty = stx;
        fy = fx;
        dy = dx;
      }
      stx = stp;
      fx = fp;
      dx = dp;
    }
    stp = stpf;
  }

  double finit, ginit, gtest, stpmax, width, width1;
  double stx = 0, fx_, gx_, sty = 0, fy_, gy_, stmin = 0, stmax;
  bool brackt = false;
  int stage = 1;
};

// Minimizes (or maximizes) a smooth function over the box lower <= x <= upper, whose bounds may be infinite, by a
// projected quasi-Newton method: the two-metric projection of Bertsekas, "Projected Newton methods for optimization
// problems with simple constraints", SIAM J. Control Optim. 20 (1982), with the L-BFGS matrix of Nocedal, "Updating
// quasi-Newton matrices with limited storage", Math. Comp. 35 (1980), in place of the Hessian, as in Schmidt, Kim and
// Sra, "Projected Newton-type methods in machine learning" (2011), Section 11.3.
//
// f(x, g) returns the value at x (always a point of the box) and fills g, already of the size of x, with the gradient.
// x is the starting point, projected onto the box if it is outside, and on return the best point found, which need
// not be where f was evaluated last. A value or gradient that is not finite rejects a trial step, so f may be
// undefined on parts of the box, though not at the starting point.
//
// In each iteration the variables on a bound whose gradient points out of the box are held there; the others move
// along the L-BFGS direction of the face they span, and a backtracking line search follows the projection of that
// line onto the box, so that variables reaching a bound stay on it while the rest of the step is taken in full.
BOOST_MATH_TEST_EXPORT template <class F>
bounded_lbfgs_result bounded_lbfgs(F &&f, std::vector<double> &x, std::vector<double> const &lower,
                                   std::vector<double> const &upper, bounded_lbfgs_parameters const &params = {}) {
  using std::abs;
  using std::isfinite;
  using std::sqrt;
  size_t const n = x.size();
  double const inf = std::numeric_limits<double>::infinity();
  double const sign = params.maximize ? -1 : 1;
  bounded_lbfgs_result result{inf, inf, 0, 0};
  auto evaluate = [&](std::vector<double> const &p, std::vector<double> &g) {
    double v = sign * f(p, g);
    ++result.evaluations;
    for (size_t i = 0; i < n; ++i) {
      g[i] *= sign;
      if (!isfinite(g[i])) {
        v = inf;
      }
    }
    return isfinite(v) ? v : inf;
  };
  auto project = [&](double v, size_t i) { return (std::min)((std::max)(v, lower[i]), upper[i]); };
  // The dot product over the free variables.
  std::vector<double> free(n);
  auto dot = [&](std::vector<double> const &a, std::vector<double> const &b) {
    double s = 0;
    for (size_t i = 0; i < n; ++i) {
      s += free[i] * a[i] * b[i];
    }
    return s;
  };
  for (size_t i = 0; i < n; ++i) {
    x[i] = project(x[i], i);
  }
  std::vector<double> g(n), xt(n), gt(n), d(n), a(params.memory), sy(params.memory);
  // The pairs s = x_{k+1} - x_k and y = g_{k+1} - g_k, oldest first.
  std::vector<std::vector<double>> S, Y;
  double fx = evaluate(x, g);
  double decrease = inf;
  while (fx < inf) {
    result.projected_gradient = 0;
    for (size_t i = 0; i < n; ++i) {
      result.projected_gradient = (std::max)(result.projected_gradient, abs(project(x[i] - g[i], i) - x[i]));
      free[i] = (x[i] <= lower[i] && g[i] > 0) || (x[i] >= upper[i] && g[i] < 0) ? 0.0 : 1.0;
    }
    if (result.projected_gradient <= params.gradient_tolerance || result.iterations >= params.max_iterations ||
        decrease <= params.function_tolerance * (std::max)(abs(fx), 1.0)) {
      break;
    }
    // d = -H g on the face of the free variables, by the two-loop recursion over the pairs that have positive
    // curvature there:
    for (size_t i = 0; i < n; ++i) {
      d[i] = -free[i] * g[i];
    }
    double gamma = 0;
    for (size_t k = S.size(); k-- > 0;) {
      double const yy = dot(Y[k], Y[k]);
      sy[k] = dot(S[k], Y[k]);
      if (!(sy[k] > 1e-10 * yy)) {
        sy[k] = 0;
        continue;
      }
      a[k] = dot(S[k], d) / sy[k];
      for (size_t i = 0; i < n; ++i) {
        d[i] -= a[k] * Y[k][i];
      }
      if (gamma == 0) {
        gamma = sy[k] / yy;
      }
    }
    for (size_t i = 0; i < n; ++i) {
      d[i] *= gamma > 0 ? gamma : 1.0;
    }
    for (size_t k = 0; k < S.size(); ++k) {
      if (sy[k] > 0) {
        double const b = dot(Y[k], d) / sy[k];
        for (size_t i = 0; i < n; ++i) {
          d[i] += (a[k] - b) * S[k][i];
        }
      }
    }
    for (size_t i = 0; i < n; ++i) {
      d[i] *= free[i];
    }
    // The projected line search, first along d, then, if no step along it is accepted or it does not descend, along
    // the projected gradient path, starting with a step of unit length:
    double ft = inf;
    bool accepted = false;
    // (after the line search of More and Thuente, if it is asked for and fails, as it does where f is not finite).
    for (int attempt = gamma > 0 && dot(g, d) < 0 ? 0 : 1; attempt < (params.strong_wolfe ? 3 : 2) && !accepted;
         ++attempt) {
      double t = 1;
      if (attempt >= 1) {
        S.clear();
        Y.clear();
        for (size_t i = 0; i < n; ++i) {
          d[i] = -free[i] * g[i];
        }
        t = 1 / sqrt(dot(d, d));
      }
      if (params.strong_wolfe && attempt < 2) {
        // As L-BFGS-B where its generalized Cauchy point meets no bound: the step is projected onto the box (cut
        // short, if that does not descend), and searched to at most the first bound it meets (lnsrlb), from a
        // first step of length at most one, then from unit steps, in at most 20 evaluations. Where every variable has
        // both bounds, the first step of L-BFGS-B is 1 instead, so the iterates are its own only if one does not.
        auto feasible = [&]() {
          double longest = 1e10;
          for (size_t i = 0; i < n; ++i) {
            // A variable on a bound which the step would leave the box through stays:
            d[i] = (d[i] < 0 && x[i] <= lower[i]) || (d[i] > 0 && x[i] >= upper[i]) ? 0.0 : d[i];
            longest = d[i] < 0 ? (std::min)(longest, (lower[i] - x[i]) / d[i]) : longest;
            longest = d[i] > 0 ? (std::min)(longest, (upper[i] - x[i]) / d[i]) : longest;
          }
          return longest;
        };
        double slope = 0;
        for (size_t i = 0; i < n; ++i) {
          slope += g[i] * (project(x[i] + d[i], i) - x[i]);
        }
        double const cut = slope < 0 ? 1.0 : (std::min)(1.0, feasible());
        slope = 0;
        for (size_t i = 0; i < n; ++i) {
          d[i] = cut < 1 ? cut * d[i] : project(x[i] + d[i], i) - x[i];
          slope += g[i] * d[i];
        }
        if (!(slope < -1e-12 * sqrt(dot(g, g) * dot(d, d)))) {
          continue;
        }
        double const stpmax = result.iterations == 0 ? (std::min)(feasible(), 1.0) : feasible();
        t = result.iterations == 0 ? 1 / sqrt(dot(d, d)) : 1.0;
        more_thuente search(fx, slope, (std::min)(t, stpmax), stpmax);
        for (int evaluations = 0; evaluations < 20 && !accepted; ++evaluations) {
          for (size_t i = 0; i < n; ++i) {
            xt[i] = project(x[i] + search.stp * d[i], i);
          }
          ft = evaluate(xt, gt);
          if (!(ft < inf)) {
            break;
          }
          double slope_t = 0;
          for (size_t i = 0; i < n; ++i) {
            slope_t += gt[i] * d[i];
          }
          accepted = !search.next(ft, slope_t);
        }
        accepted = accepted && ft < fx;
        continue;
      }
      for (int halvings = 0; halvings < 50 && !accepted; ++halvings, t /= 2) {
        double slope = 0;
        for (size_t i = 0; i < n; ++i) {
          xt[i] = project(x[i] + t * d[i], i);
          slope += g[i] * (xt[i] - x[i]);
        }
        if (slope < 0) {
          ft = evaluate(xt, gt);
          // The sufficient decrease condition of Armijo along the projection arc (Bertsekas, equation (27)):
          accepted = ft <= fx + 1e-4 * slope;
        }
      }
    }
    if (!accepted) {
      break;
    }
    if (S.size() == params.memory && !S.empty()) {
      std::rotate(S.begin(), S.begin() + 1, S.end());
      std::rotate(Y.begin(), Y.begin() + 1, Y.end());
    }
    else if (params.memory > 0) {
      S.emplace_back(n);
      Y.emplace_back(n);
    }
    for (size_t i = 0; i < n && !S.empty(); ++i) {
      S.back()[i] = xt[i] - x[i];
      Y.back()[i] = gt[i] - g[i];
    }
    decrease = fx - ft;
    x.swap(xt);
    g.swap(gt);
    fx = ft;
    ++result.iterations;
  }
  result.value = sign * fx;
  return result;
}

} BOOST_MATH_NAMESPACE_END // namespace boost::math::optimization::detail
#endif
