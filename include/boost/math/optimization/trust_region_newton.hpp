//  (C) Copyright Nick Thompson 2026.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#ifndef BOOST_MATH_OPTIMIZATION_TRUST_REGION_NEWTON_HPP
#define BOOST_MATH_OPTIMIZATION_TRUST_REGION_NEWTON_HPP
#include <boost/math/optimization/detail/differentiable_opt_utilties.hpp>
#include <boost/math/optimization/detail/gradient_opt_base.hpp>
#include <boost/math/optimization/detail/rdiff_optimization_policies.hpp>
#ifndef BOOST_MATH_BUILD_MODULE
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>
#endif

BOOST_MATH_NAMESPACE_BEGIN
namespace optimization {

/** @brief> sets the value of an argument, without recording it on a tape */
BOOST_MATH_EXPORT template<typename RealType>
struct trust_region_newton_update_policy
{
  template<typename T, std::size_t Order>
  void operator()(rdiff::rvar<T, Order>& x, RealType v)
  {
    x.template get_value_at<1>().get_value() = v;
  }
  void operator()(RealType& x, RealType v) { x = v; }
};

/**
 * @brief Trust-region Newton optimizer
 *
 * Minimizes the quadratic model g.p + p.Hp/2 over |p| <= radius by the
 * Steihaug-Toint truncated conjugate gradient method, which needs only
 * Hessian-vector products and stops at the boundary on negative curvature,
 * then accepts or rejects the step and updates the radius by the ratio of the
 * actual to the predicted reduction.
 *
 * Jorge Nocedal and Stephen J. Wright, Numerical Optimization, 2nd Edition,
 * Springer, 2006. Algorithms 4.1 and 7.2.
 *
 * @tparam> HessianVectorPolicy: computes Hv at the x of the last gradient
 * evaluation
 */
template<typename ArgumentContainer,
         typename RealType,
         class Objective,
         class InitializationPolicy,
         class ObjectiveEvalPolicy,
         class GradEvalPolicy,
         class HessianVectorPolicy>
class trust_region_newton
  : public abstract_optimizer<ArgumentContainer,
                              RealType,
                              Objective,
                              InitializationPolicy,
                              ObjectiveEvalPolicy,
                              GradEvalPolicy,
                              trust_region_newton_update_policy<RealType>,
                              trust_region_newton<ArgumentContainer,
                                                  RealType,
                                                  Objective,
                                                  InitializationPolicy,
                                                  ObjectiveEvalPolicy,
                                                  GradEvalPolicy,
                                                  HessianVectorPolicy>>
{
  using base_opt = abstract_optimizer<ArgumentContainer,
                                      RealType,
                                      Objective,
                                      InitializationPolicy,
                                      ObjectiveEvalPolicy,
                                      GradEvalPolicy,
                                      trust_region_newton_update_policy<RealType>,
                                      trust_region_newton<ArgumentContainer,
                                                          RealType,
                                                          Objective,
                                                          InitializationPolicy,
                                                          ObjectiveEvalPolicy,
                                                          GradEvalPolicy,
                                                          HessianVectorPolicy>>;

  HessianVectorPolicy hvp_;
  RealType radius_;
  ArgumentContainer x_nodes_;     // the arguments as initialized, to restore projected ones
  std::vector<RealType> x_eval_;  // where g, and the tape, were evaluated
  std::vector<RealType> xv_, x_old_, g_old_, p_, r_, d_, Hd_;

  // Steihaug-Toint CG: p_ approximately minimizes g.p + p.Hp/2 subject to
  // |p| <= radius_, and r_ = g + Hp.
  void steihaug(const std::vector<RealType>& g)
  {
    using std::sqrt;
    const size_t n = g.size();
    p_.assign(n, RealType(0));
    r_ = g;
    d_.resize(n);
    std::transform(g.begin(), g.end(), d_.begin(), [](RealType gi) { return -gi; });
    RealType rr = dot(r_, r_);
    const RealType tol = (std::min)(RealType(0.1), sqrt(sqrt(rr))) * sqrt(rr);
    for (size_t j = 0; j < n && sqrt(rr) > tol; ++j) {
      hvp_(this->objective_, this->x_, d_, Hd_);
      const RealType dHd = dot(d_, Hd_);
      const RealType alpha = rr / dHd;
      const RealType pd = dot(p_, d_), dd = dot(d_, d_), pp = dot(p_, p_);
      const RealType R2 = radius_ * radius_;
      if (!(dHd > 0) || pp + alpha * (2 * pd + alpha * dd) >= R2) {
        // the positive root of |p + tau d| = radius_
        const RealType q = sqrt(pd * pd + dd * (R2 - pp));
        const RealType tau = pd > 0 ? (R2 - pp) / (pd + q) : (q - pd) / dd;
        axpy(tau, d_, p_);
        axpy(tau, Hd_, r_);
        return;
      }
      axpy(alpha, d_, p_);
      axpy(alpha, Hd_, r_);
      const RealType rr_next = dot(r_, r_);
      const RealType beta = rr_next / rr;
      rr = rr_next;
      for (size_t i = 0; i < n; ++i) {
        d_[i] = beta * d_[i] - r_[i];
      }
    }
  }

public:
  trust_region_newton(Objective&& objective,
                      ArgumentContainer& x,
                      RealType initial_radius,
                      InitializationPolicy&& ip,
                      ObjectiveEvalPolicy&& oep,
                      GradEvalPolicy&& gep,
                      HessianVectorPolicy&& hvp)
    : base_opt(std::forward<Objective>(objective),
               x,
               std::forward<InitializationPolicy>(ip),
               std::forward<ObjectiveEvalPolicy>(oep),
               std::forward<GradEvalPolicy>(gep),
               trust_region_newton_update_policy<RealType>{})
    , hvp_(std::forward<HessianVectorPolicy>(hvp))
    , radius_(initial_radius)
    , x_nodes_(x)
  {
    if (!(initial_radius > 0 && initial_radius <= (std::numeric_limits<RealType>::max)())) {
      throw std::domain_error("trust_region_newton: the initial radius must be positive and finite");
    }
  }

  RealType radius() const { return radius_; }

  /** One trial step. On rejection x, the objective value and the gradient
   * are left unchanged and the radius shrinks. */
  void step()
  {
    auto& x = this->arguments();
    auto& g = this->gradients();
    auto& f = this->objective_value();
    const size_t n = x.size();
    xv_.resize(n);
    for (size_t i = 0; i < n; ++i) {
      xv_[i] = static_cast<RealType>(x[i]);
    }
    if (xv_ != x_eval_) { // first step, rejected step, or projected by minimize
      for (size_t i = 0; i < n; ++i) {
        x[i] = x_nodes_[i];
        this->update_(x[i], xv_[i]);
      }
      this->grad_eval_(this->objective_, x, this->obj_eval_, f, g);
      x_eval_ = xv_;
    }
    steihaug(g);
    const RealType predicted = -(dot(g, p_) + dot(r_, p_)) / 2;
    const RealType step_norm = norm_2(p_);
    if (!(step_norm > 0 && predicted > 0)) {
      return; // a stationary point, or a non-finite gradient or Hessian
    }

    const RealType f_old = f;
    x_old_ = xv_;
    g_old_ = g;
    for (size_t i = 0; i < n; ++i) {
      x_eval_[i] = x_old_[i] + p_[i];
      this->update_(x[i], x_eval_[i]);
    }
    this->grad_eval_(this->objective_, x, this->obj_eval_, f, g);

    // A change within the rounding error of f does not measure the model: trust it.
    using std::abs;
    const RealType actual = f_old - f;
    const RealType noise = 2 * std::numeric_limits<RealType>::epsilon() * abs(f_old);
    const RealType rho = abs(actual) <= noise && predicted <= noise ? RealType(1) : actual / predicted;
    if (!(rho >= RealType(0.25))) {
      radius_ = step_norm / 4;
    } else if (rho > RealType(0.75) && step_norm >= radius_ * (1 - RealType(1) / 1024)) {
      radius_ *= 2;
    }
    if (!(rho > RealType(1e-4))) {
      for (size_t i = 0; i < n; ++i) {
        this->update_(x[i], x_old_[i]);
      }
      f = f_old;
      g = g_old_;
    }
  }
};

/** x holds rdiff::rvar<RealType, 2> for the default policies. */
BOOST_MATH_EXPORT template<
  class Objective,
  typename ArgumentContainer,
  typename RealType = typename argument_container_t<ArgumentContainer>::type,
  class InitializationPolicy = tape_initializer_rvar<RealType, 2>,
  class FunctionEvalPolicy = reverse_mode_function_eval_policy<RealType, 2>,
  class GradientEvalPolicy = reverse_mode_gradient_evaluation_policy<RealType, 2>,
  class HessianVectorPolicy = reverse_mode_hessian_vector_product_policy<RealType>>
auto
make_trust_region_newton(Objective&& obj,
                         ArgumentContainer& x,
                         typename argument_container_t<ArgumentContainer>::type initial_radius = 1,
                         InitializationPolicy&& ip = InitializationPolicy{},
                         FunctionEvalPolicy&& fep = FunctionEvalPolicy{},
                         GradientEvalPolicy&& gep = GradientEvalPolicy{},
                         HessianVectorPolicy&& hvp = HessianVectorPolicy{})
{
  return trust_region_newton<ArgumentContainer,
                             typename argument_container_t<ArgumentContainer>::type,
                             std::decay_t<Objective>,
                             InitializationPolicy,
                             FunctionEvalPolicy,
                             GradientEvalPolicy,
                             HessianVectorPolicy>(
    std::forward<Objective>(obj),
    x,
    initial_radius,
    std::forward<InitializationPolicy>(ip),
    std::forward<FunctionEvalPolicy>(fep),
    std::forward<GradientEvalPolicy>(gep),
    std::forward<HessianVectorPolicy>(hvp));
}

} // namespace optimization
BOOST_MATH_NAMESPACE_END
#endif
