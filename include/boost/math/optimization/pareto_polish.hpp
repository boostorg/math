/*
 * Copyright Nick Thompson, 2026
 * Use, modification and distribution are subject to the
 * Boost Software License, Version 1.0. (See accompanying file
 * LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef BOOST_MATH_OPTIMIZATION_PARETO_POLISH_HPP
#define BOOST_MATH_OPTIMIZATION_PARETO_POLISH_HPP
#include <boost/math/tools/config.hpp>
#include <boost/math/tools/precision.hpp>
#include <boost/math/optimization/detail/common.hpp>
#include <boost/math/optimization/detail/parallel_for.hpp>
#include <boost/math/optimization/detail/pareto_polish_cells.hpp>
#include <boost/math/optimization/detail/pareto_polish_corrector.hpp>
#include <boost/math/optimization/detail/pareto_polish_subproblem.hpp>
#include <boost/math/optimization/detail/rdiff_optimization_policies.hpp>
#ifndef BOOST_MATH_BUILD_MODULE
#include <algorithm>
#include <array>
#include <cmath>
#include <iterator>
#include <limits>
#include <map>
#include <mutex>
#include <sstream>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>
#endif

BOOST_MATH_NAMESPACE_BEGIN namespace optimization {

// See the documentation for each of these.
BOOST_MATH_EXPORT template <typename ArgumentContainer> struct pareto_polish_parameters {
  using Real = typename ArgumentContainer::value_type;
  ArgumentContainer lower_bounds;
  ArgumentContainer upper_bounds;
  Real tolerance = BOOST_MATH_NAMESPACE::tools::root_epsilon<Real>(); // of the steps of the polish of a point
  size_t max_iterations = 100; // accepted steps per point; zero: the points are measured as they are
  // New points are added until the error of every cell is at most this, in objectives divided by their ranges over
  // the points returned. Zero: no new points.
  Real refinement_tolerance = Real(1) / 100;
  size_t max_new_points = 1000; // a guard on the cost of the refinement, not its criterion
  // A cell longer than this is divided whatever its error. At least a quarter of refinement_tolerance, or twice it
  // for three or more objectives: a shorter cell cannot be divided.
  Real max_edge_length = std::numeric_limits<Real>::infinity();
  // An end of the front which does not move is sought again from this fraction of the box inside the bounds.
  Real end_perturbation = Real(1) / 100;
  // With a policy for Hessian-vector products: the corrector of second order is used at a point while no more
  // than this many of its variables are free, and least_reduced_curvature is computed for the points returned.
  size_t corrector_dimension_limit = 100;
  bool curvature_diagnostic = false;
  unsigned threads = 1; // the calling thread is one of them
};

BOOST_MATH_EXPORT enum class pareto_polish_termination {
  stationarity_tolerance, // the polish stopped by `tolerance`, and weak_kkt_residual is at most a hundred times it
  no_acceptable_step,     // short of that: the line search found no step that keeps the guarantee
  iteration_cap,
  non_finite // an objective or a derivative is not finite at the point
};

BOOST_MATH_EXPORT BOOST_MATH_INLINE_CONSTEXPR size_t pareto_polish_created = static_cast<size_t>(-1);

BOOST_MATH_EXPORT enum class pareto_polish_cell_class {
  convex,     // each end lies on or above the tangent of the other: chord above the front, tangents below
  concave,    // each lies on or below it: the chord and the tangents change places
  mixed,      // neither, or an end has no tangent
  gap,        // presumed: a local search from each end into the cell stopped at a stationary point
  unresolved, // above the tolerance, and the attempts to divide it gave no point between its ends
  boundary    // three objectives: an edge where the triangles end, or a triangle presumed outside the front
};

// A cell of the front: the indices of its corners in the points of the result. For two objectives the cells are
// the segments between neighbours; for three, the triangles of a triangulation of the points, and the edges where
// it ends; for more they are, for now, the edges between Gabriel neighbours.
BOOST_MATH_EXPORT template <typename Real> struct pareto_polish_cell {
  std::vector<size_t> vertices;
  Real error;           // what the refinement compares with its tolerance, in objectives divided by their ranges
  Real dominance_bound; // a bound which needs no convexity, for two and three objectives; NaN otherwise
  pareto_polish_cell_class kind;
  // The error bounds the distance of the front from the cell if the front is connected and convex (or concave)
  // there, and under the assumptions on the tangents which the documentation lists.
  bool bound_under_convexity;
  // The multipliers at each vertex in turn, m numbers each, of the tangent there which the cell was compared with
  // (where a point has more than one, the one which faces the cell); zeros for a vertex whose tangent was not
  // used, and empty for a cell without tangents.
  std::vector<Real> multipliers;
};

BOOST_MATH_EXPORT enum class pareto_polish_status {
  discovered_cells_resolved, // every cell which is not a presumed gap is within the refinement tolerance. This
                             // says nothing of parts of the front which no point of the result is on.
  cells_unresolved,          // some cells are unresolved: see unresolved_cells
  budget_exhausted,          // max_new_points was reached first: see achieved_tolerance
  no_cells,                  // there are not two points with finite objectives to make a cell of
  not_refined                // the refinement tolerance is zero: the cells are measured, and not divided
};

// What happened, in order, for pareto_polish to record: the points by number, those of the front first.
BOOST_MATH_EXPORT enum class pareto_polish_event_kind {
  start,             // a point of the front
  achievement_step,  // an accepted step of the first kind
  nonworsening_step, // an accepted step of the second kind
  end_extension,     // a new point starts from point `first`, to decrease one objective
  refinement,        // a new point starts between the points `first` and `second` (and `third`)
  removed,           // the point is not returned: it is dominated, or too close to another
  round,             // a round of the refinement begins: the cells as they are follow
  cell,              // the cell of the points `first`, `second` (and `third`), with its class and error
  corrector_step,    // an accepted step of the first kind by the corrector of second order
  corrector_fallback // the corrector gave no acceptable step, for the reason `first` (a pareto_polish_fallback)
};

// Why the corrector of second order left an iteration to the first order step.
BOOST_MATH_EXPORT enum class pareto_polish_fallback : size_t {
  products_not_finite, // a Hessian-vector product is not finite, or not of the size asked for
  no_step,             // its system is singular, or its step leaves the bounds
  no_predicted_decrease,
  rejected // the objectives at its trial points break the guarantee or do not decrease enough
};

BOOST_MATH_EXPORT template <typename Real> struct pareto_polish_event {
  pareto_polish_event_kind kind;
  size_t point, first, second;
  std::vector<Real> x, f;
  size_t objective_evaluations, jacobian_evaluations; // of this point so far
  pareto_polish_cell_class cell_class = pareto_polish_cell_class::convex;
  Real error = 0;
  size_t hvp_evaluations = 0;           // of this point so far
  size_t third = pareto_polish_created; // the third corner of a triangle
  std::vector<Real> multipliers;        // of a cell: see pareto_polish_cell
};

BOOST_MATH_EXPORT template <typename ArgumentContainer, typename ObjectiveContainer> struct pareto_polish_point {
  using Real = typename ArgumentContainer::value_type;
  ArgumentContainer x;
  ObjectiveContainer f;
  // First order quantities of the objectives divided by the ranges of the result; neither certifies Pareto
  // optimality. NaN when there is no finite Jacobian.
  Real weak_kkt_residual = std::numeric_limits<Real>::quiet_NaN();
  Real linearized_nonworsening_gain = std::numeric_limits<Real>::quiet_NaN();
  // The multipliers of weak_kkt_residual, which sum to one: sum_i multipliers[i] (y[i] - f[i]) / ranges[i] = 0 is
  // the tangent of the front at f, in the objectives y. Empty when there is no finite Jacobian.
  std::vector<Real> multipliers;
  pareto_polish_termination termination = pareto_polish_termination::non_finite;
  // The index of the point of the input front which this one came from and weakly dominates, or
  // pareto_polish_created for a new point.
  size_t origin = 0;
  size_t objective_evaluations = 0, jacobian_evaluations = 0, hvp_evaluations = 0;
  size_t corrector_steps = 0;     // accepted steps of the corrector of second order
  size_t corrector_fallbacks = 0; // iterations which it left to the first order step
  // With curvature_diagnostic: the least eigenvalue of the Hessian of the weighted sum of the objectives on the
  // directions which change no variable on a bound and no objective to first order. A diagnostic, no certificate.
  Real least_reduced_curvature = std::numeric_limits<Real>::quiet_NaN();
  // The Jacobian is not finite at x, which is on the bounds, and that of a point just inside them stands for it:
  // the cells at this point are not bounds.
  bool jacobian_from_inside = false;
};

BOOST_MATH_EXPORT template <typename ArgumentContainer, typename ObjectiveContainer> struct pareto_polish_result {
  using Real = typename ArgumentContainer::value_type;
  using point_type = pareto_polish_point<ArgumentContainer, ObjectiveContainer>;
  std::vector<point_type> points;
  std::vector<pareto_polish_cell<Real>> cells;
  pareto_polish_status status = pareto_polish_status::not_refined;
  // The largest error of a cell which is neither a presumed gap nor of the kind boundary: what the refinement
  // reached, whatever the status.
  Real achieved_tolerance = 0;
  size_t gap_cells = 0, unresolved_cells = 0, boundary_cells = 0;
  size_t new_points = 0;
  // The points of the front which were polished to the stationarity tolerance, and those of them which are not
  // returned because another point dominates them: many of these tell of a front which is not mature.
  size_t polished_inputs = 0, dominated_inputs = 0;
  // The ranges of the objectives over the points returned with finite objectives (one where an objective is
  // constant), which the errors and residuals are divided by.
  std::vector<Real> ranges;
  // Three objectives: the unit vector r, in objectives divided by their ranges, along which the points were
  // projected to be triangulated, and along which the front over a triangle is compared with it.
  std::vector<Real> chart_direction;
  // The cost, with the points which are not returned.
  size_t objective_evaluations = 0, jacobian_evaluations = 0, hvp_evaluations = 0;
  size_t size() const { return points.size(); }
  point_type const &operator[](size_t i) const { return points[i]; }
  typename std::vector<point_type>::const_iterator begin() const { return points.begin(); }
  typename std::vector<point_type>::const_iterator end() const { return points.end(); }
};

// The default Jacobian: reverse mode autodiff of the objectives, which are then called with a
// std::vector<rvar<Real, 1>>. A policy is called as policy(objectives, x, values, J) with x a std::vector<Real>, to
// fill the values and the row-major Jacobian. The arguments are recorded on the tape of the calling thread after
// what is there, which is kept; pareto_polish takes them off again with release(), and a copy records its own.
BOOST_MATH_EXPORT template <typename Real> struct pareto_polish_reverse_mode_policy {
  template <class Objective>
  void operator()(Objective &&objectives, std::vector<Real> const &x, std::vector<Real> &values, std::vector<Real> &J) {
    reverse_mode_jacobian_evaluation_policy<Real>{}(objectives, args_.at(x), values, J);
  }
  void release() { args_.release(); }

private:
  detail::pareto_polish_arguments<Real, 1> args_;
};

namespace detail {

// The policies of a thread.
template <class JacobianPolicy, class HvpPolicy> struct pareto_polish_policies {
  JacobianPolicy jacobian;
  HvpPolicy hvp;
};

template <typename Real> struct pareto_polish_state {
  std::vector<Real> x, f, a, J, lam; // lam: the multipliers of the weak KKT residual; empty without a Jacobian
  // A new point is corrected along the ray from a in the direction dir (times the ranges), with nothing asked of
  // f <= a.
  std::vector<Real> dir;
  // For two objectives, the multipliers (s, 1 - s) for s in [lam_lo, lam_hi] have a weak residual of at most
  // slack; stationary: that is small, see certify.
  Real lam_lo = 0, lam_hi = 1, slack = std::numeric_limits<Real>::infinity();
  bool stationary = false;
  // inside: J is that of a point just inside the bounds. no_J: none is finite at x.
  bool have_J = false, no_J = false, inside = false;
  bool certified = false; // by certify, in the ranges in use
  // As the last polish left them: certify may take the termination back, for the ranges in use.
  pareto_polish_termination polished = pareto_polish_termination::non_finite;
  bool reached = false;
  bool pending = true;
  bool guard = true;
  bool leader = false; // its model of second order is for other points as well
  // A search from here into the cell on that side stopped at a stationary point: see the refinement.
  bool gap_left = false, gap_right = false;
  bool converged = false; // the last polish ended where its own scalarization is stationary
  Real radius = std::numeric_limits<Real>::infinity(); // no step changes a variable by more than this
  size_t origin = pareto_polish_created, id = 0;
  Real weak = std::numeric_limits<Real>::quiet_NaN();
  Real gain = std::numeric_limits<Real>::quiet_NaN();
  pareto_polish_termination termination = pareto_polish_termination::non_finite;
  size_t evaluations = 0, jacobians = 0, hvps = 0, corrector_steps = 0, fallbacks = 0;
  Real curvature = std::numeric_limits<Real>::quiet_NaN();
  // Three objectives: multipliers other than lam with a weak residual of at most slack, m numbers each, and
  // whether there are such in more than one direction from lam, as at a corner of the front.
  std::vector<Real> faces;
  bool wide = false;
  pareto_polish_corrector<Real> model; // of the corrector of second order, kept for the points made from this one
};

template <class Policy> auto pareto_polish_release(Policy &policy, int) -> decltype(policy.release(), void()) {
  policy.release();
}
template <class Policy> void pareto_polish_release(Policy &, long) {}

// Runs tasks on the calling thread and on threads - 1 more, each of those with its own copy of the policies. What
// the policies of the calling thread recorded on its reverse mode tapes is taken off them at the end.
template <class JacobianPolicy, class HvpPolicy> class pareto_polish_workers {
public:
  using policies_type = pareto_polish_policies<JacobianPolicy, HvpPolicy>;
  pareto_polish_workers(JacobianPolicy const &jacobian, HvpPolicy const &hvp, size_t threads)
      : policies_{jacobian, hvp}, mine_(policies_), threads_(threads) {}
  ~pareto_polish_workers() {
    pareto_polish_release(mine_.hvp, 0);
    pareto_polish_release(mine_.jacobian, 0);
  }
  // task(j, policies) for j = 0, ..., count - 1.
  template <class Task> void run(size_t count, Task task) {
    std::vector<policies_type> others((std::min)(threads_, (std::max)(count, size_t(1))) - 1, policies_);
    parallel_for(count, threads_, [&](size_t j, size_t thread) { task(j, thread == 0 ? mine_ : others[thread - 1]); });
  }

private:
  const policies_type policies_;
  policies_type mine_;
  size_t threads_;
};

template <class Func, class ArgumentContainer> class pareto_polisher {
public:
  using Real = typename ArgumentContainer::value_type;
  using Vec = std::vector<Real>;

  pareto_polisher(Func const &objectives, pareto_polish_parameters<ArgumentContainer> const &params, Vec range,
                  std::vector<pareto_polish_event<Real>> *trace)
      : objectives_(objectives), params_(params), range_(std::move(range)), n_(params.lower_bounds.size()),
        m_(range_.size()), trace_(trace) {}

  // The ranges of the objectives: the scale of the rays and of the results.
  void set_range(Vec range) { range_ = std::move(range); }

  void record(pareto_polish_state<Real> const &s, pareto_polish_event_kind kind, size_t first, size_t second,
              pareto_polish_cell_class cell_class = pareto_polish_cell_class::convex, Real error = 0,
              size_t third = pareto_polish_created, std::vector<Real> const &multipliers = {}) const {
    if (trace_) {
      std::lock_guard<std::mutex> lock(mutex_);
      trace_->push_back({kind, s.id, first, second, s.x, s.f, s.evaluations, s.jacobians, cell_class, error, s.hvps,
                         third, multipliers});
    }
  }

  // The objectives at x; false if they are not all finite.
  bool evaluate(Vec const &x, Vec &f) const {
    ArgumentContainer arg = params_.lower_bounds;
    for (size_t k = 0; k < n_; ++k) {
      arg[k] = x[k];
    }
    auto y = objectives_(arg);
    f.assign(std::begin(y), std::end(y));
    return f.size() == m_ && pareto_polish_finite(f);
  }

  // From s.x with objectives s.f: on return s.f <= s.a (the objectives at the start, unless given), by evaluation,
  // if s.guard is set. Otherwise s.a is any point of the space of the objectives.
  template <class Policy, class Hvp> void polish(pareto_polish_state<Real> &s, Policy &policy, Hvp &hvp) const {
    using std::abs;
    using event = pareto_polish_event_kind;
    s.termination = pareto_polish_termination::non_finite;
    s.converged = s.certified = false;
    if (s.a.empty()) {
      s.a = s.f;
    }
    const Real sigma = Real(1) / 10000;
    Real tau[2] = {1, 1};
    bool have_prev[2] = {false, false}, stage1_converged = false, stage1_done = false;
    Vec prev_x[2], prev_lam[2], prev_w[2], lo(n_), hi(n_), Jh, c(m_), lam, d, last, xn(n_), fn, scale = range_;
    for (size_t i = 0; i < m_ && s.dir.size() == m_; ++i) {
      scale[i] *= s.dir[i];
    }
    int stage = 1;
    size_t iterations = 0, consecutive = 0, shrinks = 0, stable = 0, first_order = 0, failures = 0;
    bool second = !std::is_same<Hvp, pareto_polish_no_hvp>::value;
    pareto_polish_corrector<Real> &newton = s.model;
    newton.resize(second ? n_ : 0, second ? m_ : 0);
    bool inherited = !newton.empty();
    Vec last_x, last_J, weights(m_);
    auto multiply = [&](Vec const &l) {
      Vec w(n_, Real(0));
      for (size_t i = 0; i < m_; ++i) {
        for (size_t k = 0; k < n_; ++k) {
          w[k] += l[i] * Jh[i * n_ + k];
        }
      }
      return w;
    };
    // max_i (c_i + g_i.v), g_i the rows of Jh: what stage 1 decreases, to first order, after the step v.
    auto linear = [&](Vec const &v) {
      Real most = -std::numeric_limits<Real>::infinity();
      for (size_t i = 0; i < m_; ++i) {
        Real gd = 0;
        for (size_t k = 0; k < n_; ++k) {
          gd += Jh[i * n_ + k] * v[k];
        }
        most = (std::max)(most, c[i] + gd);
      }
      return most;
    };
    if (params_.max_iterations == 0) {
      // No step: the point is certified as it is.
      s.termination = jacobian(s, policy) ? pareto_polish_termination::iteration_cap : s.termination;
      return;
    }
    for (;;) {
      if (!jacobian(s, policy)) {
        s.termination = pareto_polish_termination::non_finite;
        return;
      }
      Jh = s.J;
      if (second && !last_x.empty() && last_x != s.x) {
        // The secant updates of the model, from the change of each gradient over the last step.
        Vec step(n_), change(n_);
        for (size_t i = 0; i < m_; ++i) {
          for (size_t k = 0; k < n_; ++k) {
            step[k] = s.x[k] - last_x[k];
            change[k] = s.J[i * n_ + k] - last_J[i * n_ + k];
          }
          newton.secant(i, step, change);
        }
      }
      if (second) {
        last_x = s.x;
        last_J = s.J;
      }
      Real phi = -std::numeric_limits<Real>::infinity();
      for (size_t i = 0; i < m_; ++i) {
        for (size_t k = 0; k < n_; ++k) {
          Jh[i * n_ + k] /= scale[i];
        }
        c[i] = (s.f[i] - s.a[i]) / scale[i];
        phi = (std::max)(phi, c[i]);
      }
      for (size_t k = 0; k < n_; ++k) {
        lo[k] = (std::max)(params_.lower_bounds[k] - s.x[k], -s.radius);
        hi[k] = (std::min)(params_.upper_bounds[k] - s.x[k], s.radius);
      }
      const int z = stage - 1;
      const bool cone = stage == 2;
      if (have_prev[z]) {
        // Barzilai-Borwein: tau = s.s / s.y, y the change of J'lam at the multipliers of the last step.
        const Vec w = multiply(prev_lam[z]);
        Real sy = 0, ss = 0;
        for (size_t k = 0; k < n_; ++k) {
          const Real sk = s.x[k] - prev_x[z][k];
          sy += sk * (w[k] - prev_w[z][k]);
          ss += sk * sk;
        }
        if (sy > ss * std::numeric_limits<Real>::epsilon()) {
          tau[z] = (std::min)((std::max)(ss / sy, Real(1) / 100000000), Real(100000000));
        }
        have_prev[z] = false;
      }
      // Stage 1 decreases max_i (f_i - a_i)/r_i, the scalarization of Pascoletti and Serafini with reference point
      // a and direction r (the achievement function of Wierzbicki), which is at most zero where f <= a. Stage 2
      // decreases the sum of the objectives with none of them increasing.
      auto acceptable = [&]() {
        Real model = -std::numeric_limits<Real>::infinity(), phi_new = model, gain = 0, predicted = 0;
        for (size_t i = 0; i < m_; ++i) {
          if (cone ? !(fn[i] <= s.f[i]) : (s.guard && !(fn[i] <= s.a[i]))) {
            return false;
          }
          Real gd = 0;
          for (size_t k = 0; k < n_; ++k) {
            gd += Jh[i * n_ + k] * d[k];
          }
          model = (std::max)(model, c[i] + gd);
          phi_new = (std::max)(phi_new, (fn[i] - s.a[i]) / scale[i]);
          gain += (s.f[i] - fn[i]) / scale[i];
          predicted -= gd;
        }
        return cone ? (gain > 0 && gain >= sigma * predicted) : (model < phi && phi_new <= phi + sigma * (model - phi));
      };
      // The trial point s.x + d, which is taken if it is acceptable.
      auto take = [&](event what) {
        for (size_t k = 0; k < n_; ++k) {
          xn[k] = (std::min)((std::max)(s.x[k] + d[k], params_.lower_bounds[k]), params_.upper_bounds[k]);
        }
        ++s.evaluations;
        if (!evaluate(xn, fn) || !acceptable()) {
          return false;
        }
        prev_x[z] = s.x;
        prev_lam[z] = lam;
        prev_w[z] = multiply(lam);
        have_prev[z] = true;
        have_prev[1 - z] = false;
        s.x = xn;
        s.f = fn;
        s.have_J = false;
        record(s, what, s.id, s.id);
        return true;
      };
      Vec lam0, d0;
      pareto_polish_subproblem(Jh, m_, c, lo, hi, tau[z], cone, lam0, d0);
      // Stationary: the step is short, and so is that on the box itself, without the radius of the search.
      const Real limit = (cone ? 100 : 1) * params_.tolerance * tau[z];
      int result = 0; // 1: a step, 0: none found, -1: stationary
      if (pareto_polish_norm(d0) <= limit) {
        result = -1;
        if (s.radius < (std::numeric_limits<Real>::max)()) {
          Vec lo2(n_), hi2(n_);
          for (size_t k = 0; k < n_; ++k) {
            lo2[k] = params_.lower_bounds[k] - s.x[k];
            hi2[k] = params_.upper_bounds[k] - s.x[k];
          }
          pareto_polish_subproblem(Jh, m_, c, lo2, hi2, tau[z], cone, lam, d);
          result = pareto_polish_norm(d) <= limit ? -1 : 0;
        }
      }
      if (second && !cone && result == 0 && !s.inside) {
        // The corrector (see the documentation): the step of Newton's method on the KKT conditions of
        //   minimize t  subject to  (f_i(x) - a_i) / r_i <= t,  the bounds,
        // with W = sum_i lam_i H_i / r_i on the free variables, at the multipliers of the first order step. The
        // products for the model are asked for once the free variables have been the same for two iterations
        // with no new bound reached, and the first order steps are slow; a model from another point is used at once.
        bool quiet = true;
        for (size_t k = 0; k < n_; ++k) {
          quiet = quiet && !(lo[k] < hi[k] && ((d0[k] <= lo[k] && lo[k] < 0) || (d0[k] >= hi[k] && hi[k] > 0)));
        }
        const bool same = newton.free_variables(multiply(lam0), lo, hi);
        stable = quiet ? (same ? stable + 1 : 1) : 0;
        const size_t nf = newton.free().size();
        // The objectives whose Hessians the step needs: those with multipliers which are not negligible.
        const Real largest = *std::max_element(lam0.begin(), lam0.end());
        size_t unknown = 0;
        for (size_t i = 0; i < m_; ++i) {
          weights[i] = lam0[i] > largest / 1000000 ? lam0[i] / scale[i] : Real(0);
          unknown += weights[i] != 0 ? newton.missing(i).size() : 0;
        }
        inherited = inherited && unknown == 0;
        // What the first order step would gain by its model: the step of the corrector is to promise a tenth of
        // that at least, so that it is not used where it is the shorter of the two.
        const Real cauchy = linear(d0);
        const size_t slow = s.leader ? size_t(8) : (std::max)(size_t(8), m_ * nf);
        if ((inherited ? stable >= 1 : (stable >= 2 && first_order >= slow)) && nf > 0 &&
            nf <= params_.corrector_dimension_limit) {
          using reason = pareto_polish_fallback;
          reason why = reason::products_not_finite;
          for (size_t i = 0; i < m_ && second && unknown > 0; ++i) {
            const std::vector<size_t> missing = newton.missing(i);
            if (weights[i] == 0 || missing.empty()) {
              continue;
            }
            Vec unit(m_, Real(0)), directions(missing.size() * n_, Real(0)), columns;
            unit[i] = 1;
            for (size_t j = 0; j < missing.size(); ++j) {
              directions[j * n_ + missing[j]] = 1;
            }
            second = products(s, policy, hvp, unit, directions, columns);
            for (size_t j = 0; j < missing.size() && second; ++j) {
              newton.set_column(i, missing[j], &columns[j * n_]);
            }
          }
          for (int attempt = 0; attempt < 2 && second && result == 0; ++attempt) {
            Real shift = 0;
            if (attempt == 1) {
              Vec W, G;
              size_t rank;
              newton.reduced(Jh, m_, weights, Real(0), W, G);
              const Real least = pareto_polish_reduced_curvature(W, Vec(), nf, Real(0), rank);
              shift = 1 / tau[0] + (least < 0 ? -2 * least : Real(0));
            }
            why = reason::no_step;
            if (!newton.step(Jh, m_, c, lo, hi, lam0, weights, shift, d, lam)) {
              continue;
            }
            const Real model = linear(d);
            why = reason::no_predicted_decrease;
            if (!(model < phi && phi - model >= (phi - cauchy) / 10)) {
              continue;
            }
            why = reason::rejected;
            if (take(event::corrector_step)) {
              ++s.corrector_steps;
              result = 1;
            }
          }
          if (result == 0) {
            ++s.fallbacks;
            record(s, event::corrector_fallback, static_cast<size_t>(why), s.id);
            failures += why == reason::rejected;
            if (inherited && failures >= 3) {
              // The model of the neighbours does not serve here: one of this point is made, once it is due.
              newton.forget();
              inherited = false;
              failures = 0;
            }
            second = second && failures < 3;
          }
        }
      }
      bool first_tried = false;
      for (int pass = 0; pass < 2 && result == 0; ++pass) {
        // Pass 0 (stage 1): the variables which the first trial step puts on a bound stay there while the rest of
        // the step is shortened. Pass 1: all of the step is shortened.
        Vec lo2 = lo, hi2 = hi;
        if (pass == 0) {
          bool any = false;
          for (size_t k = 0; k < n_ && !cone; ++k) {
            if (lo[k] < hi[k] && (d0[k] <= lo[k] || d0[k] >= hi[k])) {
              lo2[k] = hi2[k] = d0[k];
              any = true;
            }
          }
          if (!any) {
            continue;
          }
        }
        Real t = tau[z];
        int cut = 0;
        last.clear();
        for (int halving = 0; halving <= (cone ? 6 : 25); ++halving, t /= 2) {
          if (halving == 0) {
            lam = lam0;
            d = d0;
          } else {
            pareto_polish_subproblem(Jh, m_, c, lo2, hi2, t, cone, lam, d);
          }
          if (pareto_polish_norm(d) == 0) {
            break;
          }
          if (halving == 0 && first_tried) {
            last = d; // the first pass tried this step
            continue;
          }
          first_tried = first_tried || halving == 0;
          if (d == last) {
            // The bounds cut this step to the last one, which was tried. A shorter one may differ: twelve more
            // lengths are tried from the first which does (none that was accepted needed more than ten).
            cut = cut == 0 ? 1 : cut;
            continue;
          }
          if (cut > 0 && ++cut > 13) {
            break;
          }
          last = d;
          if (take(cone ? event::nonworsening_step : event::achievement_step)) {
            tau[z] = 2 * t;
            first_order += !cone;
            result = 1;
            break;
          }
        }
      }
      if (result == 0 && cone) {
        // Where one objective is stationary on its own, as on a face of the box on which it is least, the step
        // of the second kind may move a variable off that face, or lie on the boundary of its cone with no first
        // order decrease of another objective, and then no length of it is acceptable. The step which decreases
        // the other objectives alike to first order is tried: it leaves a point which is only weakly optimal.
        Vec others(m_, Real(0)), sizes(m_, Real(0));
        Real most = 0;
        for (size_t i = 0; i < m_; ++i) {
          for (size_t k = 0; k < n_; ++k) {
            sizes[i] += lo[k] < 0 && hi[k] > 0 ? Jh[i * n_ + k] * Jh[i * n_ + k] : Real(0);
          }
          most = (std::max)(most, sizes[i]);
        }
        bool any = false;
        for (size_t i = 0; i < m_; ++i) {
          if (!(sizes[i] > most * std::numeric_limits<Real>::epsilon())) {
            others[i] = -1 / std::numeric_limits<Real>::epsilon();
            any = true;
          }
        }
        // The variables on a bound stay there: they are what holds the objective which is left out.
        Vec lo2 = lo, hi2 = hi;
        for (size_t k = 0; k < n_; ++k) {
          if (!(lo[k] < 0 && hi[k] > 0)) {
            lo2[k] = hi2[k] = 0;
          }
        }
        Real t = tau[z];
        for (int halving = 0; halving <= 6 && any && most > 0; ++halving, t /= 2) {
          pareto_polish_subproblem(Jh, m_, others, lo2, hi2, t, false, lam, d);
          if (!(pareto_polish_norm(d) > limit)) {
            break;
          }
          if (take(event::nonworsening_step)) {
            tau[z] = 2 * t;
            result = 1;
            break;
          }
        }
      }
      if (result == 0 && stage == 1 && s.radius < (std::numeric_limits<Real>::max)() && ++shrinks <= 20) {
        s.radius /= 4; // the steps of the line search do not get shorter than this
        continue;
      }
      if (result == 1) {
        stage1_done = false;
        if (++iterations >= params_.max_iterations) {
          s.termination = pareto_polish_termination::iteration_cap;
          return;
        }
        // Stage 2 goes on while its steps are accepted, and is tried after every fourth step of stage 1: where
        // the ray of stage 1 leaves the image of the box, stage 1 only creeps toward the bound.
        if (stage == 1 && ++consecutive >= 4) {
          stage = 2;
          consecutive = 0;
        }
      } else if (stage == 1) {
        stage1_done = true;
        stage1_converged = result == -1;
        stage = 2;
      } else if (stage1_done) {
        s.converged = stage1_converged;
        s.termination = (stage1_converged || result == -1) ? pareto_polish_termination::stationarity_tolerance
                                                           : pareto_polish_termination::no_acceptable_step;
        return;
      } else {
        stage = 1;
        consecutive = 0;
      }
    }
  }

  // The two first order quantities of the result at s.x, and for two objectives the tangents of the front.
  template <class Policy> void certify(pareto_polish_state<Real> &s, Policy &policy) const {
    using std::sqrt;
    s.lam.clear();
    s.stationary = false;
    // How the polish ended is kept, so that the termination is that of the ranges in use, whatever it was in
    // those of an earlier round.
    if (!s.certified) {
      s.polished = s.termination;
      s.reached = s.converged;
    }
    s.termination = s.polished;
    s.converged = s.reached;
    s.certified = true;
    if (!jacobian(s, policy)) {
      return;
    }
    Vec Jh = s.J, lo(n_), hi(n_), lam, d;
    for (size_t i = 0; i < m_; ++i) {
      for (size_t k = 0; k < n_; ++k) {
        Jh[i * n_ + k] /= range_[i];
      }
    }
    for (size_t k = 0; k < n_; ++k) {
      lo[k] = s.x[k] <= params_.lower_bounds[k] ? Real(0) : -std::numeric_limits<Real>::infinity();
      hi[k] = s.x[k] >= params_.upper_bounds[k] ? Real(0) : std::numeric_limits<Real>::infinity();
    }
    const Vec zero(m_, Real(0));
    pareto_polish_subproblem(Jh, m_, zero, lo, hi, Real(1), false, s.lam, d);
    s.weak = pareto_polish_norm(d);
    pareto_polish_subproblem(Jh, m_, zero, lo, hi, Real(1), true, lam, d);
    s.gain = pareto_polish_norm(d);
    if (s.termination == pareto_polish_termination::stationarity_tolerance && !(s.weak <= 100 * params_.tolerance)) {
      // The step of the polish was short for a bound near by, or in other ranges: the point is not stationary.
      s.termination = pareto_polish_termination::no_acceptable_step;
      s.converged = false;
    }
    if (m_ == 3) {
      std::vector<size_t> F;
      for (size_t k = 0; k < n_; ++k) {
        if (lo[k] < 0 && hi[k] > 0) {
          F.push_back(k);
        }
      }
      Vec JF(m_ * F.size()), U, sigma, V;
      for (size_t i = 0; i < m_; ++i) {
        for (size_t k = 0; k < F.size(); ++k) {
          JF[i * F.size() + k] = Jh[i * n_ + F[k]];
        }
      }
      pareto_polish_svd(std::move(JF), m_, F.size(), U, sigma, V);
      pareto_polish_faces(s, Jh, lo, hi, U, params_.tolerance, n_);
    }
    if (m_ != 2) {
      return;
    }
    // The multipliers are not unique at a kink of the front, at a corner of the box and at the ends: those with
    // a residual near that of the best form an interval, as the residual is convex in s. Each segment of the
    // front then takes the tangent which faces it, the one sided derivative of the Sandwich algorithm. The
    // residual of every multiplier of the interval is at most the slack, which the cells are charged for.
    auto residual = [&](Real t) {
      return pareto_polish_residual([&](size_t k) { return t * Jh[k] + (1 - t) * Jh[n_ + k]; }, lo, hi);
    };
    const Real largest = (std::max)(residual(Real(1)), residual(Real(0)));
    const Real noise = 8 * BOOST_MATH_NAMESPACE::tools::root_epsilon<Real>() * largest;
    const Real limit = 2 * s.weak + noise;
    s.lam_lo = pareto_polish_bisect(residual, s.lam[0], Real(0), limit, 60);
    s.lam_hi = pareto_polish_bisect(residual, s.lam[0], Real(1), limit, 60);
    s.slack = (std::max)(residual(s.lam_lo), residual(s.lam_hi));
    // Stationary: the weak residual is within the hundred tolerances at which the second stage of the polish stops.
    s.stationary = !s.inside && s.slack <= 200 * params_.tolerance + noise;
  }

  // The least curvature of the weighted sum of the objectives at the multipliers of s, a stationary point, on the
  // directions which change no variable on a bound and no objective to first order. It is left as NaN if the
  // Jacobian on the other variables is not of rank m - 1: of less, the multipliers are not unique; of more, the
  // point is not stationary.
  template <class Policy, class Hvp> void curvature(pareto_polish_state<Real> &s, Policy &policy, Hvp &hvp) const {
    using std::sqrt;
    if (s.termination != pareto_polish_termination::stationarity_tolerance || s.inside || !s.have_J ||
        s.lam.size() != m_) {
      return;
    }
    Vec Jh = s.J, lo(n_), hi(n_), weights(m_), W, G, columns;
    Real largest = 0;
    for (size_t i = 0; i < m_; ++i) {
      Real norm2 = 0;
      for (size_t k = 0; k < n_; ++k) {
        Jh[i * n_ + k] /= range_[i];
        norm2 += Jh[i * n_ + k] * Jh[i * n_ + k];
      }
      largest = (std::max)(largest, sqrt(norm2));
      weights[i] = s.lam[i] / range_[i];
    }
    for (size_t k = 0; k < n_; ++k) {
      lo[k] = s.x[k] <= params_.lower_bounds[k] ? Real(0) : Real(-1);
      hi[k] = s.x[k] >= params_.upper_bounds[k] ? Real(0) : Real(1);
    }
    pareto_polish_corrector<Real> model;
    const Vec one(1, Real(1));
    model.resize(n_, 1);
    model.free_variables(lo, lo, hi, true);
    std::vector<size_t> const &F = model.free();
    const size_t nf = F.size();
    size_t rank = 0;
    if (nf > params_.corrector_dimension_limit) {
      return;
    }
    // The rank of the Jacobian on the free variables, before anything is asked of the policy.
    Vec JF(m_ * nf), left, lengths, right;
    for (size_t i = 0; i < m_; ++i) {
      for (size_t k = 0; k < nf; ++k) {
        JF[i * nf + k] = Jh[i * n_ + F[k]];
      }
    }
    pareto_polish_svd(std::move(JF), m_, nf, left, lengths, right);
    for (auto const &v : lengths) {
      rank += v > largest / 10000;
    }
    if (rank + 1 != m_) {
      return;
    }
    Vec directions(nf * n_, Real(0));
    for (size_t j = 0; j < nf; ++j) {
      directions[j * n_ + F[j]] = 1;
    }
    if (nf > 0 && !products(s, policy, hvp, weights, directions, columns)) {
      return;
    }
    for (size_t j = 0; j < nf; ++j) {
      model.set_column(0, F[j], &columns[j * n_]);
    }
    model.reduced(Jh, m_, one, Real(0), W, G);
    s.curvature = pareto_polish_reduced_curvature(W, G, nf, largest / 10000, rank);
  }

private:
  // The products of the Hessian of the weighted sum of the objectives at s.x. False if they are not finite.
  template <class Policy, class Hvp>
  bool products(pareto_polish_state<Real> &s, Policy &, Hvp &hvp, Vec const &weights, Vec const &directions,
                Vec &out) const {
    s.hvps += directions.size() / n_;
    hvp(objectives_, s.x, weights, directions, out);
    return out.size() == directions.size() && pareto_polish_finite(out);
  }
  template <class Policy>
  bool products(pareto_polish_state<Real> &, Policy &, pareto_polish_no_hvp &, Vec const &, Vec const &, Vec &) const {
    return false;
  }
  // By differences of the Jacobian, s.J being that at s.x: the step is the square root of the precision of the
  // widest of the ranges of the variables which move, to the side on which it stays within the bounds.
  template <class Policy>
  bool products(pareto_polish_state<Real> &s, Policy &policy, pareto_polish_jacobian_difference_hvp &,
                Vec const &weights, Vec const &directions, Vec &out) const {
    using std::abs;
    using std::isfinite;
    out.assign(directions.size(), Real(0));
    Vec x(n_), values, J;
    for (size_t j = 0; (j + 1) * n_ <= directions.size(); ++j) {
      Real const *v = &directions[j * n_];
      Real width = 0, longest = 0;
      for (size_t k = 0; k < n_; ++k) {
        if (v[k] != 0) {
          const Real room = params_.upper_bounds[k] - params_.lower_bounds[k];
          width = (std::max)(width, room);
          longest = (std::max)(longest, abs(v[k]));
        }
      }
      Real h = BOOST_MATH_NAMESPACE::tools::root_epsilon<Real>() * width / longest;
      bool inside = false;
      for (int side = 0; side < 2 && !inside; ++side, h = inside ? h : -h) {
        inside = isfinite(h);
        for (size_t k = 0; k < n_; ++k) {
          x[k] = s.x[k] + h * v[k];
          inside = inside && x[k] >= params_.lower_bounds[k] && x[k] <= params_.upper_bounds[k];
        }
      }
      if (!inside) {
        return false;
      }
      ++s.jacobians;
      policy(objectives_, x, values, J);
      if (J.size() != s.J.size()) {
        return false;
      }
      for (size_t i = 0; i < m_; ++i) {
        for (size_t k = 0; k < n_ && weights[i] != 0; ++k) {
          out[j * n_ + k] += weights[i] * (J[i * n_ + k] - s.J[i * n_ + k]) / h;
        }
      }
    }
    return pareto_polish_finite(out);
  }

  // The Jacobian at s.x, a little inside the bounds if it is not finite on them. False if it is not finite.
  template <class Policy> bool jacobian(pareto_polish_state<Real> &s, Policy &policy) const {
    if (s.have_J || s.no_J) {
      return s.have_J;
    }
    Vec values, x = s.x;
    for (int attempt = 0; attempt < 2 && (attempt == 0 || x != s.x); ++attempt) {
      ++s.jacobians;
      policy(objectives_, x, values, s.J);
      if (s.J.size() == m_ * n_ && pareto_polish_finite(s.J)) {
        s.have_J = true;
        s.inside = attempt > 0;
        return true;
      }
      for (size_t k = 0; k < n_; ++k) {
        const Real h = (params_.upper_bounds[k] - params_.lower_bounds[k]) / 100000000;
        x[k] = (std::min)((std::max)(x[k], params_.lower_bounds[k] + h), params_.upper_bounds[k] - h);
      }
    }
    s.no_J = true;
    return false;
  }

  Func const &objectives_;
  pareto_polish_parameters<ArgumentContainer> const &params_;
  Vec range_;
  size_t n_, m_;
  std::vector<pareto_polish_event<Real>> *trace_;
  mutable std::mutex mutex_;
};

// A cell, with what the refinement needs; i and j are indices of states.
BOOST_MATH_TEST_EXPORT template <typename Real> struct pareto_polish_edge {
  size_t i = 0, j = 0;
  Real error = 0, length = 0, box = std::numeric_limits<Real>::quiet_NaN();
  pareto_polish_cell_class kind = pareto_polish_cell_class::convex;
  bool bound = false, resolved = true;
  bool flat = false;      // both heights are within the allowances: convex and concave alike
  int stage = 0;          // the attempts made to divide it
  // Three objectives: the third corner of a triangle, and the weights of its corners at the point to divide it at.
  size_t k = pareto_polish_created;
  Real b[3] = {Real(1) / 2, Real(1) / 2, Real(0)};
  std::vector<Real> apex; // the target of the second search beyond an edge where the triangles end
  std::vector<Real> normals; // the multipliers of its tangents, those of each corner in turn
};

// The Sandwich approximation of Burkard, Hamacher and Rote on the chord from A to B, two points in objectives
// divided by their ranges, B to the right of and below A: see the documentation for the error. sa and sb are the
// ends of the intervals of multipliers at A and at B; each gives the tangent which faces the other point. ra and
// rb, divided by the lengths of the normals, are the allowances. e.box is the largest distance from the chord of a
// point of the box between A and B; heights of different sign beyond the allowances make the cell mixed.
BOOST_MATH_TEST_EXPORT template <typename Real>
void pareto_polish_sandwich(Real const *A, Real const *B, Real const *sa, Real const *sb, Real ra, Real rb,
                            pareto_polish_edge<Real> &e) {
  using std::abs;
  using std::sqrt;
  const Real dx = B[0] - A[0], dy = B[1] - A[1];
  e.normals.assign(4, Real(0));
  auto facing = [](Real const *s, Real x, Real y, Real *normal, Real &allowance, Real *multipliers) {
    Real best = std::numeric_limits<Real>::infinity(), length = 1;
    for (int k = 0; k < 2; ++k) {
      const Real norm = sqrt(s[k] * s[k] + (1 - s[k]) * (1 - s[k]));
      const Real g = (s[k] * x + (1 - s[k]) * y) / norm;
      if (g < best) {
        best = g;
        length = norm;
        normal[0] = s[k] / norm;
        normal[1] = (1 - s[k]) / norm;
        multipliers[0] = s[k];
        multipliers[1] = 1 - s[k];
      }
    }
    allowance /= length;
    return best;
  };
  Real na[2] = {0, 1}, nb[2] = {1, 0}, ea = ra, eb = rb;
  const Real gab = facing(sa, dx, dy, na, ea, &e.normals[0]), gba = facing(sb, -dx, -dy, nb, eb, &e.normals[2]);
  const Real L = sqrt(dx * dx + dy * dy), P = abs(gab) + ea, Q = abs(gba) + eb;
  const bool convex = gab >= -ea && gba >= -eb, concave = gab <= ea && gba <= eb;
  e.length = L;
  e.box = L > 0 ? abs(dx * dy) / L : Real(0);
  e.error = e.box;
  e.bound = convex || concave;
  e.flat = convex && concave;
  if (!e.bound) {
    e.kind = pareto_polish_cell_class::mixed;
    return;
  }
  using kind = pareto_polish_cell_class;
  e.kind = convex && (!concave || gab + gba >= 0) ? kind::convex : kind::concave;
  const Real det = nb[0] * na[1] - nb[1] * na[0];
  if (abs(det) * L > 0) {
    e.error = (std::min)(e.error, P * Q / (L * abs(det)));
  }
  const Real heights[2][2] = {{abs(gab), P}, {abs(gba), Q}};
  for (auto const &h : heights) {
    if (h[0] < L) {
      e.error = (std::min)(e.error, h[1] * L / sqrt(L * L - h[0] * h[0]));
    }
  }
}

} // namespace detail

BOOST_MATH_EXPORT template <typename ArgumentContainer, typename ObjectiveContainer>
void validate_pareto_polish_parameters(pareto_polish_parameters<ArgumentContainer> const &params,
                                       std::vector<std::pair<ArgumentContainer, ObjectiveContainer>> const &front) {
  std::ostringstream oss;
  detail::validate_bounds(params.lower_bounds, params.upper_bounds);
  if (!(params.tolerance > 0) || !(params.refinement_tolerance >= 0) || !(params.end_perturbation >= 0) ||
      !(params.max_edge_length > 0)) {
    oss << __FILE__ << ":" << __LINE__ << ":" << __func__;
    oss << ": The tolerance and the largest edge length must be positive, and the refinement tolerance and end "
           "perturbation not negative.";
    throw std::domain_error(oss.str());
  }
  if (params.threads == 0) {
    oss << __FILE__ << ":" << __LINE__ << ":" << __func__;
    oss << ": There must be at least one thread.";
    throw std::invalid_argument(oss.str());
  }
  const bool two = front.empty() || std::distance(std::begin(front[0].second), std::end(front[0].second)) == 2;
  if (params.max_edge_length < (two ? params.refinement_tolerance / 4 : 2 * params.refinement_tolerance)) {
    oss << __FILE__ << ":" << __LINE__ << ":" << __func__;
    oss << ": The largest edge length must be at least a quarter of the refinement tolerance, or twice it for three "
           "or more objectives: a shorter cell cannot be divided.";
    throw std::domain_error(oss.str());
  }
  for (auto const &p : front) {
    const size_t m = static_cast<size_t>(std::distance(std::begin(p.second), std::end(p.second)));
    const size_t m0 = static_cast<size_t>(std::distance(std::begin(front[0].second), std::end(front[0].second)));
    if (m != m0 || m == 0 || m > 10) {
      oss << __FILE__ << ":" << __LINE__ << ":" << __func__;
      oss << ": Every point of the front must have the same number of objectives, between 1 and 10, but got " << m
          << " and " << m0 << ".";
      throw std::domain_error(oss.str());
    }
    bool inside = p.first.size() == params.lower_bounds.size();
    for (size_t k = 0; inside && k < p.first.size(); ++k) {
      inside = p.first[k] >= params.lower_bounds[k] && p.first[k] <= params.upper_bounds[k];
    }
    if (!inside) {
      oss << __FILE__ << ":" << __LINE__ << ":" << __func__;
      oss << ": Every point of the front must have the dimension of the bounds and lie within them.";
      throw std::domain_error(oss.str());
    }
  }
}

// Polishes a front, as nsga2 returns it, with the derivatives of the objectives, and refines it. Of the points
// with finite objectives, those of the result are mutually nondominated, and every one of the front is weakly
// dominated by one of them; a point whose objectives are not finite is returned as it is. See the documentation.
BOOST_MATH_EXPORT template <class Func, typename ArgumentContainer, typename ObjectiveContainer,
                            class JacobianPolicy =
                                pareto_polish_reverse_mode_policy<typename ArgumentContainer::value_type>,
                            class HvpPolicy = pareto_polish_no_hvp>
pareto_polish_result<ArgumentContainer, ObjectiveContainer>
pareto_polish(const Func objectives, std::vector<std::pair<ArgumentContainer, ObjectiveContainer>> const &front,
              pareto_polish_parameters<ArgumentContainer> const &params, JacobianPolicy policy = JacobianPolicy(),
              std::vector<pareto_polish_event<typename ArgumentContainer::value_type>> *trace = nullptr,
              HvpPolicy hvp = HvpPolicy()) {
  using Real = typename ArgumentContainer::value_type;
  using Vec = std::vector<Real>;
  using std::abs;
  using std::isfinite;
  using std::sqrt;
  validate_pareto_polish_parameters(params, front);
  pareto_polish_result<ArgumentContainer, ObjectiveContainer> result;
  const Real tol = params.refinement_tolerance;
  if (front.empty()) {
    result.status = tol > 0 ? pareto_polish_status::no_cells : pareto_polish_status::not_refined;
    return result;
  }
  const size_t n = params.lower_bounds.size();
  const size_t m = static_cast<size_t>(std::distance(std::begin(front[0].second), std::end(front[0].second)));
  using kind = pareto_polish_event_kind;
  using cell_class = pareto_polish_cell_class;
  using state = detail::pareto_polish_state<Real>;
  using edge = detail::pareto_polish_edge<Real>;
  std::vector<state> states(front.size());
  for (size_t i = 0; i < front.size(); ++i) {
    states[i].x.assign(std::begin(front[i].first), std::end(front[i].first));
    states[i].f.assign(std::begin(front[i].second), std::end(front[i].second));
    states[i].origin = states[i].id = i;
  }
  // The ranges of the objectives (one where an objective is constant) over the points which would be returned now:
  // those with finite objectives which are not `out`, that is, which no other dominates and which do not repeat
  // the objectives of an earlier one.
  std::vector<char> out, finite;
  auto ranges = [&]() {
    const size_t N = states.size();
    Vec low(m, std::numeric_limits<Real>::infinity()), range(m, Real(1));
    Vec high(m, -std::numeric_limits<Real>::infinity());
    out.assign(N, 0);
    finite.resize(N);
    for (size_t e = 0; e < N; ++e) {
      finite[e] = detail::pareto_polish_finite(states[e].f);
    }
    for (size_t e = 0; e < N; ++e) {
      for (size_t k = 0; k < N && finite[e] && !out[e]; ++k) {
        out[e] = k != e && finite[k] &&
                 (detail::dominates(states[k].f, states[e].f) || (k < e && states[k].f == states[e].f));
      }
      for (size_t i = 0; i < m && finite[e] && !out[e]; ++i) {
        low[i] = (std::min)(low[i], states[e].f[i]);
        high[i] = (std::max)(high[i], states[e].f[i]);
      }
    }
    for (size_t i = 0; i < m; ++i) {
      if (high[i] > low[i] && isfinite(high[i] - low[i])) {
        range[i] = high[i] - low[i];
      }
    }
    return range;
  };
  const Vec range0 = ranges();
  detail::pareto_polisher<Func, ArgumentContainer> polisher(objectives, params, range0, trace);
  for (auto const &s : states) {
    polisher.record(s, kind::start, s.id, s.id);
  }
  size_t next_id = states.size();
  using workers_type = detail::pareto_polish_workers<JacobianPolicy, HvpPolicy>;
  using policies_type = typename workers_type::policies_type;
  workers_type workers(policy, hvp, params.threads);
  // Polishes the states which are waiting.
  auto work = [&]() {
    std::vector<state *> pending, leaders;
    for (auto &s : states) {
      if (s.pending && detail::pareto_polish_finite(s.f)) {
        pending.push_back(&s);
      }
      s.pending = false;
    }
    auto run = [&](std::vector<state *> const &batch) {
      workers.run(batch.size(), [&](size_t j, policies_type &mine) {
        polisher.polish(*batch[j], mine.jacobian, mine.hvp);
      });
    };
    const size_t stride = 8;
    if (!std::is_same<HvpPolicy, pareto_polish_no_hvp>::value && pending.size() > stride / 2) {
      // With the corrector of second order every eighth point, in the order of the first objective, is polished
      // first, and each of the others then starts from the model of the nearest of these which has one.
      std::stable_sort(pending.begin(), pending.end(), [](state const *a, state const *b) { return a->f[0] < b->f[0]; });
      std::vector<state *> followers;
      for (size_t j = 0; j < pending.size(); ++j) {
        pending[j]->leader = j % stride == stride / 2;
        (pending[j]->leader ? leaders : followers).push_back(pending[j]);
      }
      run(leaders);
      for (state *s : followers) {
        Real nearest = std::numeric_limits<Real>::infinity();
        state const *from = nullptr;
        for (state const *leader : leaders) {
          Real d2 = 0;
          for (size_t k = 0; k < n; ++k) {
            const Real width = params.upper_bounds[k] - params.lower_bounds[k];
            d2 += width > 0 ? (s->x[k] - leader->x[k]) * (s->x[k] - leader->x[k]) / (width * width) : Real(0);
          }
          if (!leader->model.empty() && d2 < nearest) {
            nearest = d2;
            from = leader;
          }
        }
        if (from) {
          s->model = from->model;
        }
      }
      pending = followers;
    }
    run(pending);
  };
  work();
  Vec range = ranges();
  auto distance2 = [&](Vec const &a, Vec const &b) {
    Real d = 0;
    for (size_t i = 0; i < m; ++i) {
      d += (a[i] - b[i]) * (a[i] - b[i]) / (range[i] * range[i]);
    }
    return d;
  };
  auto spent = [&](state const &s) {
    result.objective_evaluations += s.evaluations;
    result.jacobian_evaluations += s.jacobians;
    result.hvp_evaluations += s.hvps;
  };
  // The new points from `first` on which are finite, not dominated by an earlier point, and not within `close` of
  // one, stay.
  auto admit = [&](size_t first, Real close) {
    size_t kept = first;
    for (size_t e = first; e < states.size(); ++e) {
      bool ok = detail::pareto_polish_finite(states[e].f);
      for (size_t k = 0; ok && k < kept; ++k) {
        ok = !detail::pareto_polish_finite(states[k].f) ||
             (!detail::dominates(states[k].f, states[e].f) && distance2(states[k].f, states[e].f) >= close * close);
      }
      if (!ok) {
        polisher.record(states[e], kind::removed, states[e].id, states[e].id);
        spent(states[e]);
      } else if (kept++ != e) {
        states[kept - 1] = std::move(states[e]);
      }
    }
    states.resize(kept);
    return kept - first;
  };
  size_t created = 0;
  if (tol > 0 && m > 1 && params.max_new_points >= m) {
    // The ends of the front: from its best point in each objective, that objective alone is decreased (the others
    // are anchored far above), and the steps of the second kind then take back what the others can give.
    const size_t N = states.size();
    auto anchor = [&](state &s, size_t i) {
      s.a = s.f;
      for (size_t j = 0; j < m; ++j) {
        if (j != i) {
          s.a[j] += 1000000 * range0[j];
        }
      }
    };
    std::vector<size_t> best(m, N);
    for (size_t i = 0; i < m; ++i) {
      for (size_t k = 0; k < N; ++k) {
        if (detail::pareto_polish_finite(states[k].f) && (best[i] == N || states[k].f[i] < states[best[i]].f[i])) {
          best[i] = k;
        }
      }
      state s;
      s.id = next_id++;
      s.x = states[best[i] < N ? best[i] : 0].x;
      s.f = states[best[i] < N ? best[i] : 0].f;
      anchor(s, i);
      polisher.record(s, kind::end_extension, states[best[i] < N ? best[i] : 0].id, i);
      states.push_back(std::move(s));
    }
    work();
    for (size_t i = 0; i < m && params.end_perturbation > 0; ++i) {
      state &s = states[N + i];
      if (best[i] < N && distance2(s.f, states[best[i]].f) < tol * tol) {
        // It did not move: once more from a little inside the bounds.
        // Nothing to do where no variable is on a bound.
        Vec f, x = s.x;
        for (size_t k = 0; k < n; ++k) {
          const Real h = params.end_perturbation * (params.upper_bounds[k] - params.lower_bounds[k]);
          x[k] = x[k] <= params.lower_bounds[k] ? x[k] + h : (x[k] >= params.upper_bounds[k] ? x[k] - h : x[k]);
        }
        if (x == s.x) {
          continue;
        }
        ++s.evaluations;
        if (polisher.evaluate(x, f)) {
          s.x = x;
          s.f = f;
          anchor(s, i);
          s.have_J = s.no_J = false;
          s.pending = true;
          polisher.record(s, kind::end_extension, states[best[i]].id, i);
        }
      }
    }
    work();
    created += admit(N, tol);
  }
  // The cells of the points which are returned as things stand, in objectives divided by their ranges, which the
  // points are certified in: see the documentation, and pareto_polish_sandwich for two objectives.
  std::map<std::pair<size_t, size_t>, int> attempts;
  // Three objectives: the triangles, and what the searches into them have found.
  detail::pareto_polish_triangulation<Real> mesh;
  const size_t none = pareto_polish_created;
  bool planar = false; // the cells of this round are those of the triangulation
  const int max_attempts = m == 2 ? 3 : 1;
  const Real longest = m == 2 ? params.max_edge_length : (std::min)(params.max_edge_length, Real(1) / 4);
  Vec certified_range;
  auto build = [&]() {
    const size_t N = states.size();
    range = ranges();
    std::vector<size_t> order, waiting;
    for (size_t e = 0; e < N; ++e) {
      if (finite[e] && (range != certified_range || !states[e].certified)) {
        waiting.push_back(e);
      }
      if (finite[e] && !out[e]) {
        order.push_back(e);
      }
    }
    certified_range = range;
    polisher.set_range(range);
    workers.run(waiting.size(), [&](size_t j, policies_type &mine) {
      polisher.certify(states[waiting[j]], mine.jacobian);
    });
    std::vector<edge> edges;
    if (m == 2) {
      std::sort(order.begin(), order.end(), [&](size_t a, size_t b) { return states[a].f[0] < states[b].f[0]; });
      for (size_t q = 0; q + 1 < order.size(); ++q) {
        state const &a = states[order[q]], &b = states[order[q + 1]];
        const Real A[2] = {a.f[0] / range[0], a.f[1] / range[1]}, B[2] = {b.f[0] / range[0], b.f[1] / range[1]};
        const bool tangents = a.lam.size() == 2 && b.lam.size() == 2;
        const Real sa[2] = {a.lam_lo, a.lam_hi}, sb[2] = {b.lam_lo, b.lam_hi};
        const Real apart = detail::pareto_polish_distance(a.x, b.x);
        edge e;
        e.i = order[q];
        e.j = order[q + 1];
        detail::pareto_polish_sandwich(A, B, sa, sb, tangents ? a.slack * apart : Real(0),
                                       tangents ? b.slack * apart : Real(0), e);
        if (!tangents) {
          e.kind = cell_class::mixed;
          e.normals.clear();
        }
        e.bound = e.bound && tangents && a.stationary && b.stationary;
        edges.push_back(e);
      }
      // Where a convex cell meets a concave one the curvature changes sign in one of the two, or between them.
      // Neither is then a bound, nor are the flat cells between them, which otherwise take the class before them.
      for (size_t c = 0, last = edges.size(); c < edges.size(); ++c) {
        if (edges[c].kind == cell_class::mixed) {
          last = edges.size();
        } else if (edges[c].flat && last < edges.size()) {
          edges[c].kind = edges[last].kind;
        } else if (!edges[c].flat) {
          for (size_t k = last; k <= c && last < edges.size() && edges[last].kind != edges[c].kind; ++k) {
            edges[k].resolved = false;
          }
          last = c;
        }
      }
      for (auto &e : edges) {
        e.bound = e.bound && e.resolved;
        if (!e.bound) {
          e.error = e.box;
        }
        e.resolved = e.error <= tol && e.length <= longest;
      }
    }
    if (m == 3) {
      mesh.cells(states, order, range, tol, longest, params.max_edge_length, edges);
    }
    // Three objectives with no triangle (two points, or all of them on one line of the chart): edges, as for more.
    planar = m == 3 && !edges.empty();
    for (size_t p = 0; p < order.size() && m > 2 && !planar; ++p) {
      for (size_t q = p + 1; q < order.size(); ++q) {
        state const &a = states[order[p]], &b = states[order[q]];
        const Real dij = distance2(a.f, b.f);
        bool gabriel = true;
        for (size_t k = 0; gabriel && k < order.size(); ++k) {
          gabriel = k == p || k == q || !(distance2(a.f, states[order[k]].f) + distance2(b.f, states[order[k]].f) < dij);
        }
        if (gabriel) {
          edge e;
          e.i = order[p];
          e.j = order[q];
          e.length = e.error = sqrt(dij);
          e.kind = cell_class::mixed;
          if (a.lam.size() == m && b.lam.size() == m) {
            Real sab = 0, sba = 0;
            for (size_t i = 0; i < m; ++i) {
              sab += a.lam[i] * (b.f[i] - a.f[i]) / range[i];
              sba += b.lam[i] * (a.f[i] - b.f[i]) / range[i];
            }
            e.normals = a.lam;
            e.normals.insert(e.normals.end(), b.lam.begin(), b.lam.end());
            sab /= detail::pareto_polish_norm(a.lam);
            sba /= detail::pareto_polish_norm(b.lam);
            if ((sab < 0) == (sba < 0)) {
              e.kind = sab < 0 ? cell_class::concave : cell_class::convex;
            }
            e.error = e.kind == cell_class::mixed ? (std::max)(abs(sab), abs(sba)) : (std::min)(abs(sab), abs(sba));
          }
          e.resolved = e.error <= tol && e.length <= longest;
          edges.push_back(e);
        }
      }
    }
    Real worst = 0;
    state marker;
    marker.id = edges.size();
    polisher.record(marker, kind::round, created, created);
    for (auto &e : edges) {
      auto const before = attempts.find(std::make_pair(states[e.i].id, states[e.j].id));
      e.stage = planar ? e.stage : (before == attempts.end() ? 0 : before->second);
      if (m == 2 && states[e.i].gap_right && states[e.j].gap_left) {
        e.kind = cell_class::gap;
        e.error = e.box;
      } else if (!planar && !e.resolved && e.stage >= max_attempts) {
        e.kind = cell_class::unresolved;
      }
      if (e.kind == cell_class::gap || e.kind == cell_class::unresolved) {
        e.bound = false;
        e.resolved = true;
      }
      worst = e.kind == cell_class::gap || e.kind == cell_class::boundary ? worst : (std::max)(worst, e.error);
      polisher.record(marker, kind::cell, states[e.i].id, states[e.j].id, e.kind, e.error,
                      e.k == none ? none : states[e.k].id, e.normals);
    }
    result.achieved_tolerance = worst;
    return edges;
  };
  // The refinement: each round gives the cells which are not resolved a new point, or moves an end into them.
  std::vector<edge> edges;
  bool open = false;
  for (;;) {
    edges = build();
    std::vector<edge> todo;
    for (auto const &e : edges) {
      if (!e.resolved) {
        todo.push_back(e);
      }
    }
    open = !todo.empty();
    if (!open || !(tol > 0) || created >= params.max_new_points) {
      break;
    }
    // The worst first. An error which is not a number counts as the largest, so that this is an ordering.
    std::sort(todo.begin(), todo.end(), [](edge const &a, edge const &b) {
      const Real most = std::numeric_limits<Real>::infinity(), ea = a.error == a.error ? a.error : most;
      const Real eb = b.error == b.error ? b.error : most;
      return ea > eb || (ea == eb && (std::make_pair(a.i, a.j) < std::make_pair(b.i, b.j) ||
                                      (a.i == b.i && a.j == b.j && a.k < b.k)));
    });
    todo.resize((std::min)(todo.size(), params.max_new_points - created));
    const size_t N = states.size();
    const Real close = m <= 3 ? tol / 8 : tol;
    // A cell, and 0: a point between its ends, 1: its left end moved, 2: then its right end.
    std::vector<std::pair<size_t, int>> tasks;
    for (auto const &e : todo) {
      int stage = e.stage;
      if (planar) {
        // A point of its own image first, then one sought along r.
        if (e.k != none) {
          mesh.tried[mesh.corners(states, e)] = stage + 1;
        } else {
          mesh.sought[(std::minmax)(states[e.i].id, states[e.j].id)] = stage + 1;
        }
        tasks.emplace_back(&e - todo.data(), stage);
        continue;
      }
      stage += m == 2 && stage == 1 && states[e.i].gap_right;
      stage += m == 2 && stage == 2 && states[e.j].gap_left;
      attempts[std::make_pair(states[e.i].id, states[e.j].id)] = (std::min)(stage + 1, max_attempts);
      if (stage < max_attempts) {
        tasks.emplace_back(&e - todo.data(), stage);
      }
    }
    states.resize(N + tasks.size());
    // A new point whose objectives are finite is polished and certified; any other is not kept.
    auto finish = [&](state &s, bool ok, policies_type &mine) {
      if (ok) {
        polisher.polish(s, mine.jacobian, mine.hvp);
        polisher.certify(s, mine.jacobian);
      } else {
        s.f.assign(m, std::numeric_limits<Real>::quiet_NaN());
      }
    };
    workers.run(tasks.size(), [&](size_t task, policies_type &mine) {
      const edge &e = todo[tasks[task].first];
      const int end = tasks[task].second;
      state const &a = states[e.i], &b = states[e.j];
      state &s = states[N + task];
      s.id = next_id + task;
      s.x.resize(n);
      bool ok = false;
      if (planar) {
        // From the point of the domain with the weights of the cell: as a point of the front is polished, or
        // (the second attempt for a triangle) towards the point of the triangle with those weights, along r.
        state const &c = states[e.k == none ? e.j : e.k];
        Vec target(3);
        for (size_t k = 0; k < n; ++k) {
          s.x[k] = e.b[0] * a.x[k] + e.b[1] * b.x[k] + e.b[2] * c.x[k];
        }
        for (size_t i = 0; i < 3; ++i) {
          target[i] = e.b[0] * a.f[i] + e.b[1] * b.f[i] + e.b[2] * c.f[i];
        }
        ++s.evaluations;
        ok = polisher.evaluate(s.x, s.f);
        if (end == 1) {
          for (size_t i = 0; i < 3 && e.k == none; ++i) {
            target[i] = e.apex[i] * range[i];
          }
          s.a = target;
          s.dir.assign(mesh.chart.begin(), mesh.chart.begin() + 3);
          s.guard = false;
        }
        polisher.record(s, kind::refinement, a.id, b.id, cell_class::convex, Real(0), e.k == none ? none : c.id);
        s.model = a.model;
        return finish(s, ok, mine);
      }
      if (end != 0) {
        s.x = end == 1 ? a.x : b.x;
        s.f = end == 1 ? a.f : b.f;
        s.a = {b.f[0], a.f[1]};
        s.dir = {end == 1 ? Real(1) / 1000 : Real(1), end == 1 ? Real(1) : Real(1) / 1000};
        s.guard = false;
        // In short steps, so that it does not leave the piece of the front which it is on.
        for (size_t k = 0; k < n; ++k) {
          s.radius = (std::min)(s.radius, -abs(b.x[k] - a.x[k]) / 20);
        }
        s.radius = -s.radius;
        ok = true;
      } else {
        // The midpoint of the ends in the domain.
        for (size_t k = 0; k < n; ++k) {
          s.x[k] = a.x[k] + (b.x[k] - a.x[k]) / 2;
        }
        ++s.evaluations;
        ok = polisher.evaluate(s.x, s.f);
      }
      polisher.record(s, kind::refinement, a.id, b.id);
      if (end == 0 && !a.model.empty() && !b.model.empty()) {
        s.model.blend(a.model, Real(1) / 2, b.model, Real(1) / 2);
      } else {
        s.model = (end == 2 || (end == 0 && a.model.empty())) ? b.model : a.model;
      }
      finish(s, ok, mine);
    });
    for (size_t task = 0; task < tasks.size(); ++task) {
      // Where a search from an end converged, that end, or the point which it came to, is the end of a gap.
      const int end = tasks[task].second;
      edge const &e = todo[tasks[task].first];
      state &s = states[N + task], &from = states[end == 1 ? e.i : e.j];
      const size_t ia = states[e.i].id, ib = states[e.j].id;
      const bool found = detail::pareto_polish_finite(s.f);
      Vec y[4];
      for (size_t c = 0; c < 4 && planar; ++c) {
        state const &corner =
            c == 0 ? states[e.i] : (c == 1 ? states[e.j] : (c == 2 ? states[e.k == none ? e.j : e.k] : s));
        y[c].assign(3, Real(0));
        for (size_t i = 0; i < 3 && (c < 3 || found); ++i) {
          y[c][i] = corner.f[i] / range[i];
        }
      }
      if (planar && e.k == none) {
        // An edge where the triangles end: the distance of the new point from it estimates how far the front
        // is from it between its ends, and a quarter of that for the two edges which the point makes of it.
        Real d = std::numeric_limits<Real>::infinity(), along = 0, length2 = 0;
        for (size_t i = 0; i < 3 && found; ++i) {
          along += (y[3][i] - y[0][i]) * (y[1][i] - y[0][i]);
          length2 += (y[1][i] - y[0][i]) * (y[1][i] - y[0][i]);
        }
        if (found) {
          const Real t = (std::min)((std::max)(along / length2, Real(0)), Real(1));
          Vec foot(3);
          for (size_t i = 0; i < 3; ++i) {
            foot[i] = y[0][i] + t * (y[1][i] - y[0][i]);
          }
          d = detail::pareto_polish_distance(y[3], foot);
          mesh.sagitta[(std::minmax)(ia, s.id)] = mesh.sagitta[(std::minmax)(ib, s.id)] = d / 4;
        }
        mesh.sagitta[(std::minmax)(ia, ib)] = d;
        // The search beyond the edge stopped at a point of the front on the side of its triangle, farther from
        // the edge than the tolerance: the front ends there, and the triangle of that point and the edge (the
        // point itself, or the one which it repeats) is presumed outside the front.
        Real side = 0;
        for (size_t w = 0; w < 2 && found; ++w) {
          Real outward = 0, inward = 0;
          for (size_t j = 0; j < 3; ++j) {
            outward += (e.apex[j] - (y[0][j] + y[1][j]) / 2) * mesh.chart[3 + 3 * w + j];
            inward += (y[3][j] - (y[0][j] + y[1][j]) / 2) * mesh.chart[3 + 3 * w + j];
          }
          side += outward * inward;
        }
        if (end == 1 && found && s.converged && d > tol && side < 0) {
          state const *third = &s;
          for (size_t k = 0; k < N; ++k) {
            if (k != e.i && k != e.j && detail::pareto_polish_finite(states[k].f) &&
                distance2(states[k].f, s.f) < close * close) {
              third = &states[k];
            }
          }
          std::array<size_t, 3> key{{ia, ib, third->id}};
          std::sort(key.begin(), key.end());
          mesh.tried[key] = -1;
          Vec outside = states[e.i].f;
          outside.insert(outside.end(), states[e.j].f.begin(), states[e.j].f.end());
          outside.insert(outside.end(), third->f.begin(), third->f.end());
          mesh.beyond.push_back(outside);
        }
      } else if (planar && end == 1 && found && s.converged) {
        // The search along r stopped where its own scalarization is stationary. If that is not on the line through
        // its target, the line does not meet the front: the triangle is presumed outside the front there.
        Real off[2] = {0, 0};
        for (size_t i = 0; i < 3; ++i) {
          const Real d = y[3][i] - (e.b[0] * y[0][i] + e.b[1] * y[1][i] + e.b[2] * y[2][i]);
          off[0] += d * mesh.chart[3 + i];
          off[1] += d * mesh.chart[6 + i];
        }
        // Nor does it if the point is on the line and a point of the front dominates it: it is weakly optimal only.
        bool dominated = false;
        for (size_t k = 0; k < N && !dominated; ++k) {
          dominated = detail::pareto_polish_finite(states[k].f) && detail::dominates(states[k].f, s.f);
        }
        if (dominated || sqrt(off[0] * off[0] + off[1] * off[1]) > e.length / 100) {
          mesh.tried[mesh.corners(states, e)] = -1;
          Vec outside;
          for (size_t c = 0; c < 3; ++c) {
            state const &corner = c == 0 ? states[e.i] : (c == 1 ? states[e.j] : states[e.k]);
            outside.insert(outside.end(), corner.f.begin(), corner.f.end());
          }
          mesh.beyond.push_back(outside);
        }
      }
      if (m == 2 && end != 0 && found) {
        (end == 1 ? s.gap_right : s.gap_left) = s.converged;
        if (s.converged && distance2(s.f, from.f) < close * close) {
          (end == 1 ? from.gap_right : from.gap_left) = true;
        }
        attempts[end == 1 ? std::make_pair(s.id, ib) : std::make_pair(ia, s.id)] = end + 1;
      }
    }
    next_id += tasks.size();
    created += admit(N, close);
  }
  if (params.curvature_diagnostic && !std::is_same<HvpPolicy, pareto_polish_no_hvp>::value) {
    workers.run(states.size(), [&](size_t j, policies_type &mine) {
      if (!out[j]) {
        polisher.curvature(states[j], mine.jacobian, mine.hvp);
      }
    });
  }
  // What dominates a point which is not returned is itself returned or dominated, so every point of the front
  // with finite objectives is weakly dominated by one that is returned.
  std::vector<size_t> index(states.size(), 0);
  result.ranges = range;
  if (mesh.chart.size() >= 3) {
    result.chart_direction.assign(mesh.chart.begin(), mesh.chart.begin() + 3);
  }
  for (size_t e = 0; e < states.size(); ++e) {
    auto const &s = states[e];
    spent(s);
    typename pareto_polish_result<ArgumentContainer, ObjectiveContainer>::point_type r;
    r.x = front[0].first;
    r.f = front[0].second;
    std::copy(s.x.begin(), s.x.end(), std::begin(r.x));
    std::copy(s.f.begin(), s.f.end(), std::begin(r.f));
    r.weak_kkt_residual = s.weak;
    r.multipliers = s.lam;
    r.linearized_nonworsening_gain = s.gain;
    r.termination = s.termination;
    r.origin = s.origin;
    r.objective_evaluations = s.evaluations;
    r.jacobian_evaluations = s.jacobians;
    r.hvp_evaluations = s.hvps;
    r.corrector_steps = s.corrector_steps;
    r.corrector_fallbacks = s.fallbacks;
    r.least_reduced_curvature = s.curvature;
    r.jacobian_from_inside = s.inside;
    const bool input = s.origin != pareto_polish_created;
    const bool polished = input && s.termination == pareto_polish_termination::stationarity_tolerance;
    result.polished_inputs += polished;
    if (!out[e]) {
      index[e] = result.points.size();
      result.points.push_back(std::move(r));
      result.new_points += !input;
      continue;
    }
    polisher.record(s, kind::removed, s.id, s.id);
    bool dominated = false;
    for (size_t k = 0; k < states.size() && !dominated && input; ++k) {
      dominated = finite[k] && detail::dominates(states[k].f, s.f);
    }
    result.dominated_inputs += dominated && polished;
  }
  for (auto const &e : edges) {
    result.cells.push_back({{index[e.i], index[e.j]}, e.error, e.box, e.kind, e.bound, e.normals});
    if (e.k != none) {
      result.cells.back().vertices.push_back(index[e.k]);
    }
    result.gap_cells += e.kind == cell_class::gap;
    result.unresolved_cells += e.kind == cell_class::unresolved;
    result.boundary_cells += e.kind == cell_class::boundary;
  }
  using status = pareto_polish_status;
  result.status = !(tol > 0)      ? status::not_refined
                  : edges.empty() ? status::no_cells
                  : open          ? status::budget_exhausted
                                  : (result.unresolved_cells > 0 ? status::cells_unresolved
                                                                 : status::discovered_cells_resolved);
  return result;
}

} BOOST_MATH_NAMESPACE_END // namespace boost::math::optimization
#endif
