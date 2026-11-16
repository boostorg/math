/*
 * Copyright Nick Thompson, 2026
 * Use, modification and distribution are subject to the
 * Boost Software License, Version 1.0. (See accompanying file
 * LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef BOOST_MATH_OPTIMIZATION_DETAIL_PARETO_POLISH_CORRECTOR_HPP
#define BOOST_MATH_OPTIMIZATION_DETAIL_PARETO_POLISH_CORRECTOR_HPP
#include <boost/math/tools/config.hpp>
#include <boost/math/tools/detail/linear_algebra/symmetric_eigen.hpp>
#include <boost/math/optimization/detail/pareto_polish_cells.hpp>
#include <boost/math/optimization/detail/pareto_polish_subproblem.hpp>
#include <boost/math/optimization/detail/rdiff_optimization_policies.hpp>
#ifndef BOOST_MATH_BUILD_MODULE
#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>
#endif

BOOST_MATH_NAMESPACE_BEGIN namespace optimization {

// No Hessian-vector products: the polish is of first order.
BOOST_MATH_EXPORT struct pareto_polish_no_hvp {};

// Hessian-vector products by forward differences of the Jacobian policy, within the bounds: one more Jacobian for
// each, counted as such, and nothing more is asked of the objectives. They are good to about the square root of
// the precision.
BOOST_MATH_EXPORT struct pareto_polish_jacobian_difference_hvp {};

namespace detail {

// The arguments of a reverse mode policy, recorded once between two checkpoints on the tapes of the thread which
// uses them: what was recorded before is kept, and release() leaves the tapes as they were. A copy starts with
// none, as those of the original are on the tapes of the thread which recorded them.
template <typename Real, std::size_t Order> class pareto_polish_arguments {
public:
  pareto_polish_arguments() = default;
  pareto_polish_arguments(pareto_polish_arguments const &) noexcept {}
  pareto_polish_arguments &operator=(pareto_polish_arguments const &) noexcept { return *this; }
  std::vector<rdiff::rvar<Real, Order>> &at(std::vector<Real> const &x) {
    if (tape_ != &rdiff::get_active_tape<Real, Order>() || args_.size() != x.size()) {
      release(); // of another dimension: that recording goes first
      checkpoint();
      args_.assign(x.begin(), x.end());
      checkpoint();
      tape_ = &rdiff::get_active_tape<Real, Order>();
    }
    for (size_t k = 0; k < x.size(); ++k) {
      set(args_[k], x[k]);
    }
    return args_;
  }
  // For the thread which recorded them, once nothing recorded after them is in use.
  void release() {
    for (int twice = 0; twice < 2 && tape_ == &rdiff::get_active_tape<Real, Order>(); ++twice) {
      for_each_tape<Real>([](auto &tape) { tape.pop_checkpoint(); }, std::integral_constant<std::size_t, Order>{});
    }
    args_.clear();
    tape_ = nullptr;
  }

private:
  static void checkpoint() {
    for_each_tape<Real>([](auto &tape) { tape.add_checkpoint(); }, std::integral_constant<std::size_t, Order>{});
  }
  static void set(rdiff::rvar<Real, 1> &a, Real const &v) { a.set_item(v); }
  static void set(rdiff::rvar<Real, 2> &a, Real const &v) { a.get_value().get_value() = v; }
  std::vector<rdiff::rvar<Real, Order>> args_;
  void const *tape_ = nullptr;
};

} // namespace detail

// Hessian-vector products by reverse over reverse autodiff of the objectives, which are then called with a
// std::vector<rvar<Real, 2>>. Such a policy is called as policy(objectives, x, weights, directions, products) to
// fill the products H v, one after the other, of the Hessian H of sum_i weights[i] f_i at x with the vectors v
// which follow one another in `directions`. For the tapes see pareto_polish_reverse_mode_policy.
BOOST_MATH_EXPORT template <typename Real> struct pareto_polish_reverse_mode_hvp_policy {
  template <class Objective>
  void operator()(Objective &&objectives, std::vector<Real> const &x, std::vector<Real> const &weights,
                  std::vector<Real> const &directions, std::vector<Real> &products) {
    const size_t n = x.size();
    products.clear();
    if (n == 0 || directions.size() % n != 0) {
      return; // no products: pareto_polish takes that as a failure of the policy
    }
    auto &args = args_.at(x);
    auto scalar = [&](std::vector<rdiff::rvar<Real, 2>> const &arguments) {
      auto y = objectives(arguments);
      rdiff::rvar<Real, 2> sum = Real(0);
      size_t i = 0;
      for (auto const &yi : y) {
        sum += weights[i++] * yi;
      }
      return sum;
    };
    // Only what is recorded here, after the checkpoints of the arguments, is cleared and swept: the adjoints of
    // what the caller recorded before are nodes of the first order tape which must stay as they are.
    using first = rdiff::rvar<Real, 1>;
    auto &tape = rdiff::get_active_tape<Real, 2>();
    for_each_tape<Real>([](auto &t) { t.rewind_to_last_checkpoint(); }, std::integral_constant<std::size_t, 2>{});
    rdiff::rvar<Real, 2> value = scalar(args);
    for (auto it = tape.last_checkpoint(); it != tape.end(); ++it) {
      *it->get_adjoint_ptr() = first(Real(0));
    }
    for (auto &xi : args) {
      xi.adjoint() = first(Real(0));
    }
    value.adjoint() = first(Real(1));
    for (auto it = tape.end(), begin = tape.last_checkpoint(); it != begin;) {
      (--it)->backward();
    }
    std::vector<Real> v(n), Hv;
    for (size_t j = 0; (j + 1) * n <= directions.size(); ++j) {
      v.assign(directions.begin() + j * n, directions.begin() + (j + 1) * n);
      reverse_mode_hessian_vector_product_policy<Real>{}(scalar, args, v, Hv);
      products.insert(products.end(), Hv.begin(), Hv.end());
    }
  }
  void release() { args_.release(); }

private:
  detail::pareto_polish_arguments<Real, 2> args_;
};

namespace detail {

// The least eigenvalue of the symmetric n by n matrix W on the subspace orthogonal to the rows of the k by n matrix
// G, whose singular values above `floor` are counted in `rank`. Infinity if
// that subspace is the origin, and NaN if W or G is not finite. Both are scaled first, so that no square overflows.
BOOST_MATH_TEST_EXPORT template <typename Real>
Real pareto_polish_reduced_curvature(std::vector<Real> W, std::vector<Real> G, size_t n, Real floor, size_t &rank) {
  using std::abs;
  using std::isfinite;
  using std::sqrt;
  const size_t k = n == 0 ? 0 : G.size() / n;
  std::vector<Real> Q, Wq(n), eigenvalues;
  Real scale = 0, longest_entry = 0;
  bool finite = true;
  for (auto const &v : W) {
    scale = (std::max)(scale, abs(v));
    finite = finite && isfinite(v);
  }
  for (auto const &v : G) {
    longest_entry = (std::max)(longest_entry, abs(v));
    finite = finite && isfinite(v);
  }
  rank = 0;
  if (!finite) {
    return std::numeric_limits<Real>::quiet_NaN();
  }
  scale = scale > 0 ? scale : Real(1);
  longest_entry = longest_entry > 0 ? longest_entry : Real(1);
  for (auto &v : W) {
    v /= scale;
  }
  for (auto &v : G) {
    v /= longest_entry;
  }
  floor /= longest_entry;
  // The rows of G span what its right singular vectors of singular values above the floor do.
  std::vector<Real> left, lengths;
  pareto_polish_svd(std::move(G), k, n, left, lengths, Q);
  while (rank < k && rank < n && lengths[rank] > floor) {
    ++rank;
  }
  Q.resize(rank * n);
  if (rank == n) {
    return std::numeric_limits<Real>::infinity();
  }
  // (I - qq')W(I - qq') for each q, and then the rows of Q are moved above every eigenvalue of that.
  Real frobenius = 0;
  for (auto const &v : W) {
    frobenius += v * v;
  }
  const Real sigma = 1 + 2 * sqrt(frobenius);
  for (size_t r = 0; r < rank; ++r) {
    Real const *q = &Q[r * n];
    Real qWq = 0;
    for (size_t i = 0; i < n; ++i) {
      Wq[i] = 0;
      for (size_t j = 0; j < n; ++j) {
        Wq[i] += W[i * n + j] * q[j];
      }
      qWq += q[i] * Wq[i];
    }
    for (size_t i = 0; i < n; ++i) {
      for (size_t j = 0; j < n; ++j) {
        W[i * n + j] += -q[i] * Wq[j] - Wq[i] * q[j] + qWq * q[i] * q[j];
      }
    }
  }
  for (size_t r = 0; r < rank; ++r) {
    for (size_t i = 0; i < n; ++i) {
      for (size_t j = 0; j < n; ++j) {
        W[i * n + j] += sigma * Q[r * n + i] * Q[r * n + j];
      }
    }
  }
  if (!tools::detail::symmetric_eigendecomposition(std::move(W), n, Q, eigenvalues)) {
    return std::numeric_limits<Real>::quiet_NaN();
  }
  return scale * *std::min_element(eigenvalues.begin(), eigenvalues.end());
}

// The second order model of the corrector of pareto_polish: approximations of the Hessians H_i of the objectives,
// column by column. The columns come from Hessian-vector products; from then on each is kept by the symmetric rank
// one update from the changes of the gradient of its objective.
template <typename Real> class pareto_polish_corrector {
public:
  using Vec = std::vector<Real>;
  void resize(size_t n, size_t m) {
    if (n != n_ || m != H_.size()) {
      n_ = n;
      H_.assign(m, Vec()); // each n by n, once a column of it is known
      known_.assign(m, std::vector<char>(n, 0));
    }
  }
  // The model between two points: a combination of theirs, in the columns which both know.
  void blend(pareto_polish_corrector const &a, Real wa, pareto_polish_corrector const &b, Real wb) {
    if (a.n_ != b.n_ || a.H_.size() != b.H_.size()) {
      return;
    }
    resize(a.n_, a.H_.size());
    for (size_t i = 0; i < H_.size(); ++i) {
      const bool both = !a.H_[i].empty() && !b.H_[i].empty();
      for (size_t k = 0; k < n_; ++k) {
        known_[i][k] = both && a.known_[i][k] && b.known_[i][k];
      }
      H_[i].assign(both ? n_ * n_ : 0, Real(0));
      for (size_t k = 0; k < H_[i].size(); ++k) {
        H_[i][k] = wa * a.H_[i][k] + wb * b.H_[i][k];
      }
    }
  }
  bool empty() const {
    for (auto const &row : known_) {
      if (std::find(row.begin(), row.end(), char(1)) != row.end()) {
        return false;
      }
    }
    return true;
  }
  void forget() {
    for (auto &row : known_) {
      std::fill(row.begin(), row.end(), 0);
    }
  }

  // The free variables: those with room which are not on a bound that w = J'lam pushes them into (every variable
  // on a bound, if `strict`). True if they are the ones of the last call.
  bool free_variables(Vec const &w, Vec const &lo, Vec const &hi, bool strict = false) {
    std::vector<size_t> F;
    for (size_t k = 0; k < n_; ++k) {
      const bool lower = lo[k] == 0 && (strict || w[k] > 0), upper = hi[k] == 0 && (strict || w[k] < 0);
      if (lo[k] < hi[k] && !lower && !upper) {
        F.push_back(k);
      }
    }
    const bool same = F == free_;
    free_ = std::move(F);
    return same;
  }
  std::vector<size_t> const &free() const { return free_; }

  // The free variables whose column of H_i is not known.
  std::vector<size_t> missing(size_t i) const {
    std::vector<size_t> out;
    for (size_t k : free_) {
      if (!known_[i][k]) {
        out.push_back(k);
      }
    }
    return out;
  }
  void set_column(size_t i, size_t k, Real const *column) {
    H_[i].resize(n_ * n_, Real(0));
    std::copy(column, column + n_, H_[i].begin() + k * n_);
    known_[i][k] = 1;
  }

  // H_i s = y for the step s, if that is a step in the variables whose columns are known.
  void secant(size_t i, Vec const &s, Vec const &y) {
    using std::abs;
    using std::sqrt;
    Vec &H = H_[i];
    std::vector<char> const &known = known_[i];
    if (H.empty()) {
      return;
    }
    Vec r(n_, Real(0));
    Real rs = 0, rr = 0, ss = 0;
    bool any = false;
    for (size_t k = 0; k < n_; ++k) {
      if (!known[k] && s[k] != 0) {
        return;
      }
      any = any || known[k];
    }
    for (size_t k = 0; k < n_ && any; ++k) {
      if (known[k]) {
        r[k] = y[k];
        for (size_t j = 0; j < n_; ++j) {
          r[k] -= known[j] ? (H[k * n_ + j] + H[j * n_ + k]) / 2 * s[j] : Real(0);
        }
        rs += r[k] * s[k];
        rr += r[k] * r[k];
        ss += s[k] * s[k];
      }
    }
    if (!(abs(rs) > sqrt(rr * ss) / 100000000)) {
      return;
    }
    for (size_t k = 0; k < n_; ++k) {
      for (size_t j = 0; j <= k; ++j) {
        if (known[k] && known[j]) {
          H[k * n_ + j] = H[j * n_ + k] = (H[k * n_ + j] + H[j * n_ + k]) / 2 + r[k] * r[j] / rs;
        }
      }
    }
  }

  // W = sum_i weights[i] H_i on the free variables, with `shift` added to its diagonal, and the m rows of J on them.
  void reduced(Vec const &J, size_t m, Vec const &weights, Real shift, Vec &W, Vec &G) const {
    const size_t nf = free_.size();
    W.assign(nf * nf, Real(0));
    G.assign(m * nf, Real(0));
    for (size_t a = 0; a < nf; ++a) {
      for (size_t b = 0; b < nf; ++b) {
        for (size_t i = 0; i < weights.size(); ++i) {
          if (weights[i] != 0 && !H_[i].empty()) {
            W[a * nf + b] += weights[i] * (H_[i][free_[a] * n_ + free_[b]] + H_[i][free_[b] * n_ + free_[a]]) / 2;
          }
        }
      }
      W[a * nf + a] += shift;
      for (size_t i = 0; i < m; ++i) {
        G[i * nf + a] = J[i * n_ + free_[a]];
      }
    }
  }

  // The step of Newton's method for
  //   minimize t + d'(W + shift I)d / 2  subject to  c_i + g_i.d <= t,  d = 0 off the free variables,
  // g_i the rows of J and W as in `reduced`: the objectives with positive multipliers hold with equality, starting
  // from those of lam0,
  // one being added or left out at a time until the multipliers are not negative and the others are below t. The
  // step is shortened to stay within lo <= d <= hi. False if there is no such step.
  bool step(Vec const &J, size_t m, Vec const &c, Vec const &lo, Vec const &hi, Vec const &lam0, Vec const &weights,
            Real shift, Vec &d, Vec &lam) const {
    using std::abs;
    const size_t nf = free_.size();
    Vec W, G, a, sol;
    reduced(J, m, weights, shift, W, G);
    std::vector<char> active(m, 0);
    const Real largest = *std::max_element(lam0.begin(), lam0.end());
    for (size_t i = 0; i < m; ++i) {
      active[i] = lam0[i] > largest / 1000000;
    }
    bool found = false;
    for (size_t iteration = 0; iteration < 2 * m + 2 && !found; ++iteration) {
      std::vector<size_t> A;
      for (size_t i = 0; i < m; ++i) {
        if (active[i]) {
          A.push_back(i);
        }
      }
      // Unknowns d, lam_A, t: W d + G_A' lam = 0, G_A d - t = -c_A, sum lam = 1.
      const size_t k = A.size(), K = nf + k + 1, width = K + 1;
      if (k == 0) {
        return false;
      }
      a.assign(K * width, Real(0));
      for (size_t r = 0; r < nf; ++r) {
        for (size_t b = 0; b < nf; ++b) {
          a[r * width + b] = W[r * nf + b];
        }
        for (size_t i = 0; i < k; ++i) {
          a[r * width + nf + i] = a[(nf + i) * width + r] = G[A[i] * nf + r];
        }
      }
      for (size_t i = 0; i < k; ++i) {
        a[(nf + i) * width + nf + k] = a[(nf + k) * width + nf + i] = -1;
        a[(nf + i) * width + K] = -c[A[i]];
      }
      a[(nf + k) * width + K] = -1;
      if (!pareto_polish_solve(a, K, sol)) {
        return false;
      }
      const Real t = sol[nf + k];
      size_t worst = m;
      Real most = 0, size = abs(t);
      for (size_t i = 0; i < k; ++i) {
        if (sol[nf + i] < most && k > 1) {
          most = sol[nf + i];
          worst = A[i];
        }
      }
      if (worst < m) {
        active[worst] = 0;
        continue;
      }
      for (size_t i = 0; i < m; ++i) {
        Real v = c[i];
        for (size_t r = 0; r < nf; ++r) {
          v += G[i * nf + r] * sol[r];
        }
        size = (std::max)(size, abs(c[i]));
        if (!active[i] && v - t > most) {
          most = v - t;
          worst = i;
        }
      }
      if (worst < m && most > 64 * std::numeric_limits<Real>::epsilon() * size) {
        active[worst] = 1;
        continue;
      }
      lam.assign(m, Real(0));
      for (size_t i = 0; i < k; ++i) {
        lam[A[i]] = (std::max)(sol[nf + i], Real(0));
      }
      found = true;
    }
    if (!found) {
      return false;
    }
    Real alpha = 1;
    d.assign(n_, Real(0));
    for (size_t r = 0; r < nf; ++r) {
      const size_t v = free_[r];
      d[v] = sol[r];
      if (d[v] < lo[v]) {
        alpha = (std::min)(alpha, lo[v] / d[v]);
      } else if (d[v] > hi[v]) {
        alpha = (std::min)(alpha, hi[v] / d[v]);
      }
    }
    if (!(alpha >= Real(1) / 4)) {
      return false;
    }
    Real norm = 0;
    for (size_t r = 0; r < nf; ++r) {
      const size_t v = free_[r];
      d[v] = (std::min)((std::max)(alpha * d[v], lo[v]), hi[v]);
      norm += d[v] * d[v];
    }
    using std::isfinite;
    return isfinite(norm) && norm > 0;
  }

private:
  size_t n_ = 0;
  std::vector<Vec> H_; // column k of H_i in H_[i][k n, (k + 1) n)
  std::vector<std::vector<char>> known_;
  std::vector<size_t> free_;
};

} // namespace detail
} BOOST_MATH_NAMESPACE_END // namespace boost::math::optimization
#endif
