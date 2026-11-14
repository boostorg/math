/*
 * Copyright Nick Thompson, 2026
 * Use, modification and distribution are subject to the
 * Boost Software License, Version 1.0. (See accompanying file
 * LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef BOOST_MATH_OPTIMIZATION_DETAIL_HYPERVOLUME_HPP
#define BOOST_MATH_OPTIMIZATION_DETAIL_HYPERVOLUME_HPP

#include <boost/math/tools/config.hpp>
#ifndef BOOST_MATH_BUILD_MODULE
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <numeric>
#include <type_traits>
#include <vector>
#endif

BOOST_MATH_NAMESPACE_BEGIN namespace optimization::detail {

// Lacour, R., Klamroth, K., Fonseca, C. M. (2017). A box decomposition algorithm to compute the hypervolume indicator.
// Computers & Operations Research, 79, 347-360.
// All objectives are minimized, and the reference point bounds the region of interest from above: only points which
// are less than the reference in every objective count. A point with a NaN objective is therefore ignored, as are
// dominated points and duplicates. Objectives may be equal to one another. They may not be -infinity.

// The local upper bounds of a set of points: the corners u of the maximal boxes {y < u} which hold no point. Their
// union is the region no point dominates. Bound i is u[i*m + j], and coordinate j of the point which defines it in
// dimension k is z[(i*m + k)*m + j].
template <class Real> struct local_upper_bounds {
  size_t m;
  std::vector<Real> u, z, new_u, new_z;

  // A single bound, the reference, defined by dummy points which are equal to `least` but for one coordinate.
  local_upper_bounds(std::vector<Real> const &reference, std::vector<Real> const &least)
      : m(reference.size()), u(reference) {
    for (size_t k = 0; k < m; ++k) {
      z.insert(z.end(), least.begin(), least.end());
      z[k * m + k] = reference[k];
    }
  }

  size_t size() const { return u.size() / m; }

  // Algorithm 1, but with p[j] >= rather than > in the test for a new bound, so that points may share coordinates;
  // bounds whose boxes are empty may then appear. The point must not dominate one that was inserted before it.
  void insert(Real const *p) {
    new_u.clear();
    new_z.clear();
    for (size_t i = 0; i < size();) {
      Real const *ui = u.data() + i * m;
      size_t k = 0;
      while (k < m && p[k] < ui[k]) {
        ++k;
      }
      if (k < m) {
        ++i;
        continue;
      }
      // p is in the box of bound i, which is replaced by up to m bounds, each with p in place of a defining point:
      Real const *zi = z.data() + i * m * m;
      for (size_t j = 0; j < m; ++j) {
        bool keep = true;
        for (k = 0; k < m && keep; ++k) {
          keep = k == j || p[j] >= zi[k * m + j];
        }
        if (keep) {
          new_u.insert(new_u.end(), ui, ui + m);
          new_u[new_u.size() - m + j] = p[j];
          new_z.insert(new_z.end(), zi, zi + m * m);
          std::copy(p, p + m, new_z.data() + new_z.size() - (m - j) * m);
        }
      }
      if (i + 1 < size()) {
        std::copy(u.data() + u.size() - m, u.data() + u.size(), u.data() + i * m);
        std::copy(z.data() + z.size() - m * m, z.data() + z.size(), z.data() + i * m * m);
      }
      u.resize(u.size() - m);
      z.resize(z.size() - m * m);
    }
    u.insert(u.end(), new_u.begin(), new_u.end());
    z.insert(z.end(), new_z.begin(), new_z.end());
  }

  // Inserts points[i*m + j], in lexicographic order, so that a point comes after those which dominate it, and changes
  // nothing. Any such order will do, but of those we tried, the last objective as the most significant gives the
  // fewest boxes: about half as many as the first objective as the most significant, for four to six objectives.
  void insert(std::vector<Real> const &points) {
    std::vector<size_t> order(points.size() / m);
    std::iota(order.begin(), order.end(), size_t(0));
    Real const *p = points.data();
    std::sort(order.begin(), order.end(), [&](size_t a, size_t b) {
      size_t j = m - 1;
      while (j > 0 && p[a * m + j] == p[b * m + j]) {
        --j;
      }
      return p[a * m + j] < p[b * m + j];
    });
    for (size_t i : order) {
      insert(p + i * m);
    }
  }

  // Equation 2: the region which the points dominate, below the reference, is the disjoint union of one box for each
  // bound. Calls f(lower, upper) with the corners of each box that is not empty.
  template <class F> void for_each_box(Real reference0, F &&f) const {
    std::vector<Real> lower(m), upper(m);
    for (size_t i = 0; i < size(); ++i) {
      Real const *zi = z.data() + i * m * m;
      std::copy(u.data() + i * m, u.data() + (i + 1) * m, upper.begin());
      upper[0] = reference0;
      lower[0] = zi[0];
      bool empty = !(lower[0] < upper[0]);
      for (size_t j = 1; j < m && !empty; ++j) {
        lower[j] = zi[j];
        for (size_t k = 1; k < j; ++k) {
          lower[j] = (std::max)(lower[j], zi[k * m + j]);
        }
        empty = !(lower[j] < upper[j]);
      }
      if (!empty) {
        f(lower.data(), upper.data());
      }
    }
  }
};

// The points which are less than the reference in every objective, one after the other.
template <class ObjectiveContainer>
auto below_reference(std::vector<ObjectiveContainer> const &points, ObjectiveContainer const &reference) {
  using Real = typename ObjectiveContainer::value_type;
  const size_t m = reference.size();
  std::vector<Real> flat;
  for (auto const &p : points) {
    size_t k = 0;
    while (k < m && p[k] < reference[k]) {
      flat.push_back(p[k++]);
    }
    if (k < m) {
      flat.resize(flat.size() - k);
    }
  }
  return flat;
}

// The hypervolume indicator: the volume of the region which is dominated by one of the points and is below the
// reference.
BOOST_MATH_TEST_EXPORT template <class ObjectiveContainer>
auto hypervolume(std::vector<ObjectiveContainer> const &points, ObjectiveContainer const &reference) {
  using Real = typename ObjectiveContainer::value_type;
  const size_t m = reference.size();
  local_upper_bounds<Real> bounds(std::vector<Real>(reference.begin(), reference.end()),
                                  std::vector<Real>(m, -std::numeric_limits<Real>::infinity()));
  bounds.insert(below_reference(points, reference));
  Real volume = 0;
  bounds.for_each_box(reference[0], [&](Real const *lower, Real const *upper) {
    Real box = 1;
    for (size_t j = 0; j < m; ++j) {
      box *= upper[j] - lower[j];
    }
    volume += box;
  });
  return volume;
}

// Disjoint boxes: box i is the product over j of the intervals from lower[i*m + j], which may be -infinity, to
// upper[i*m + j], which is at most the reference.
BOOST_MATH_TEST_EXPORT template <class Real> struct box_decomposition {
  size_t objectives = 0;
  std::vector<Real> lower, upper;
  size_t size() const { return objectives ? lower.size() / objectives : 0; }
};

// The region below the reference which no point dominates, as disjoint boxes. As in Yang, K., Emmerich, M., Deutz,
// A., Back, T. (2019), Efficient computation of expected hypervolume improvement using box decomposition algorithms,
// Journal of Global Optimization, 75, 3-34, and BoTorch's FastNondominatedPartitioning: this region is the union of
// the boxes below the local upper bounds of the points, so that it is, once the objectives are negated, the region
// which the bounds dominate.
BOOST_MATH_TEST_EXPORT template <class ObjectiveContainer>
auto nondominated_boxes(std::vector<ObjectiveContainer> const &points, ObjectiveContainer const &reference) {
  using Real = typename ObjectiveContainer::value_type;
  const size_t m = reference.size();
  const Real infinity = std::numeric_limits<Real>::infinity();
  std::vector<Real> r(reference.begin(), reference.end());
  local_upper_bounds<Real> bounds(r, std::vector<Real>(m, -infinity));
  bounds.insert(below_reference(points, reference));
  for (auto &x : bounds.u) {
    x = -x;
  }
  for (auto &x : r) {
    x = -x;
  }
  local_upper_bounds<Real> negated(std::vector<Real>(m, infinity), r);
  negated.insert(bounds.u);
  box_decomposition<Real> boxes;
  boxes.objectives = m;
  negated.for_each_box(infinity, [&](Real const *lower, Real const *upper) {
    for (size_t j = 0; j < m; ++j) {
      boxes.lower.push_back(-upper[j]);
      boxes.upper.push_back(-lower[j]);
    }
  });
  return boxes;
}

// The hypervolume that y adds to that of the points whose nondominated region the boxes are: the volume of the part
// of the boxes which y dominates. Zero if y has a NaN.
BOOST_MATH_TEST_EXPORT template <class ObjectiveContainer>
auto hypervolume_improvement(ObjectiveContainer const &y,
                             box_decomposition<typename ObjectiveContainer::value_type> const &boxes) {
  using Real = typename ObjectiveContainer::value_type;
  const size_t m = boxes.objectives;
  Real improvement = 0;
  for (size_t i = 0; i < boxes.size(); ++i) {
    Real volume = 1;
    for (size_t j = 0; j < m && volume > 0; ++j) {
      const Real side = boxes.upper[i * m + j] - (y[j] > boxes.lower[i * m + j] ? y[j] : boxes.lower[i * m + j]);
      // A NaN in y must not pass for -infinity:
      volume = side > 0 && y[j] == y[j] ? volume * side : 0;
    }
    improvement += volume;
  }
  return improvement;
}

// The smooth approximations of Ament, S., Daulton, S., Eriksson, D., Balandat, M., Bakshy, E. (2023), Unexpected
// improvements to expected improvement for Bayesian optimization, NeurIPS 36, with the constants of BoTorch, whose
// utils/safe_math.py and qLogExpectedHypervolumeImprovement are the reference implementation.

// log(tau f(x/tau)) and its derivative, where f(t) = log(1 + e^t) + (1/10)/(1 + t^2) is a smooth max(t, 0) which
// is positive, increasing and convex. As x -> -infinity it decays only as 1/(10 t^2), so that the derivative of its
// logarithm, -2/x, still says which way to go.
BOOST_MATH_TEST_EXPORT template <class Real> Real log_fatplus(Real x, Real tau, Real &derivative) {
  using std::exp;
  using std::log;
  using std::log1p;
  const Real t = x / tau;
  if (t < 0 && t * t * std::numeric_limits<Real>::epsilon() > 2) {
    // f(t) rounds to 1/(10 t^2), which may underflow:
    derivative = -2 / x;
    return log(tau / 10) - 2 * log(-t);
  }
  // e^-|t|, which is too small to count beyond 50:
  const Real e = t < -50 || t > 50 ? Real(0) : exp(t > 0 ? -t : t);
  const Real c = 1 / (1 + t * t);
  const Real f = (t > 0 ? t : Real(0)) + log1p(e) + c / 10;
  derivative = ((t > 0 ? Real(1) : e) / (1 + e) - t * c * c / 5) / (tau * f);
  return log(tau * f);
}

// A smooth min(a, b), min(a, b) - tau log(1 + p(|a - b|/tau)) with p(d) = 2/(2 + 2d + d^2), and its derivative with
// respect to a, which decays only as 2/d^3 where a is the greater.
BOOST_MATH_TEST_EXPORT template <class Real> Real fatmin(Real a, Real b, Real tau, Real &derivative) {
  using std::log1p;
  const Real d = (a < b ? b - a : a - b) / tau;
  const Real p = 2 / (2 + 2 * d + d * d);
  derivative = p * p * (1 + d) / (1 + p);
  derivative = a < b ? 1 - derivative : derivative;
  return (a < b ? a : b) - tau * log1p(p);
}

// log(sum_i exp(x_i)), one finite term at a time, without overflow or underflow. add(x, rescale) returns the weight
// of x, once the weights it returned before, and any sums of them, have been multiplied by rescale; the weights then
// add up to `sum`, so that sum_i (weight_i/sum) dx_i is the derivative of value().
BOOST_MATH_TEST_EXPORT template <class Real> struct log_sum_exp {
  Real largest = -std::numeric_limits<Real>::infinity(), sum = 0;
  Real add(Real x, Real &rescale) {
    using std::exp;
    rescale = x > largest ? exp(largest - x) : Real(1);
    const Real weight = x > largest ? Real(1) : exp(x - largest);
    largest = x > largest ? x : largest;
    sum = sum * rescale + weight;
    return weight;
  }
  Real value() const {
    using std::log;
    return largest + log(sum);
  }
};

// A smooth log(hypervolume_improvement(y, boxes)), which is finite and has a gradient, written to gradient[0..m)
// unless it is null, of which no component is positive, even where y improves on nothing. It is the q = 1 case of
// BoTorch's qLogExpectedHypervolumeImprovement (fat = True) for a single sample, with the objectives negated:
//   log sum_boxes prod_k exp(fatmin(log_fatplus(u_k - y_k; tau_relu), log(u_k - l_k); tau_max)),
// since u_k - max(l_k, y_k) = min(u_k - y_k, u_k - l_k). As there, a lower corner at -infinity is moved to -10^10, or
// -10^8 in single precision. The relative error of the improvement is less than m tau_max log(2) where y improves
// on the front by much more than tau_relu.
BOOST_MATH_TEST_EXPORT template <class Real>
Real log_hypervolume_improvement(Real const *y, box_decomposition<Real> const &boxes,
                                 std::common_type_t<Real> *gradient = nullptr,
                                 std::common_type_t<Real> tau_relu = Real(1e-6),
                                 std::common_type_t<Real> tau_max = Real(1e-2)) {
  using std::log;
  const size_t m = boxes.objectives;
  const Real least = std::numeric_limits<Real>::digits > 24 ? Real(-1e10) : Real(-1e8);
  // The derivatives of the logarithm of the volume of one box:
  Real few[8];
  std::vector<Real> many(m > 8 ? m : 0);
  Real *slope = m > 8 ? many.data() : few;
  log_sum_exp<Real> total;
  if (gradient) {
    std::fill(gradient, gradient + m, Real(0));
  }
  for (size_t i = 0; i < boxes.size(); ++i) {
    Real const *lower = boxes.lower.data() + i * m;
    Real const *upper = boxes.upper.data() + i * m;
    Real log_volume = 0;
    for (size_t k = 0; k < m; ++k) {
      Real d_side, d_min;
      const Real log_side = log_fatplus(upper[k] - y[k], tau_relu, d_side);
      log_volume += fatmin(log_side, log(upper[k] - (std::max)(lower[k], least)), tau_max, d_min);
      slope[k] = -d_side * d_min;
    }
    Real rescale;
    const Real weight = total.add(log_volume, rescale);
    for (size_t k = 0; gradient && k < m; ++k) {
      gradient[k] = gradient[k] * rescale + weight * slope[k];
    }
  }
  for (size_t k = 0; gradient && k < m; ++k) {
    gradient[k] /= total.sum;
  }
  return total.value();
}

} BOOST_MATH_NAMESPACE_END // namespace boost::math::optimization::detail
#endif
