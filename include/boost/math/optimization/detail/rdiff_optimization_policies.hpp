//           Copyright Maksym Zhelyenzyakov 2025-2026.
// Distributed under the Boost Software License, Version 1.0.
//      (See accompanying file LICENSE_1_0.txt or copy at
//           https://www.boost.org/LICENSE_1_0.txt)
#ifndef BOOST_MATH_OPTIMIZATION_DETAIL_RDIFF_OPTIMIZATION_POLICIES_HPP
#define BOOST_MATH_OPTIMIZATION_DETAIL_RDIFF_OPTIMIZATION_POLICIES_HPP

#include <boost/math/differentiation/autodiff_reverse.hpp>
#ifndef BOOST_MATH_BUILD_MODULE
#include <random>
#include <type_traits>
#endif

BOOST_MATH_NAMESPACE_BEGIN
namespace optimization {

namespace rdiff = BOOST_MATH_NAMESPACE::differentiation::reverse_mode;

/** @brief> applies f to the active tapes of orders Order, Order - 1, ..., 1
 */
template<typename RealType, class F>
void
for_each_tape(F&&, std::integral_constant<std::size_t, 0>)
{
}
template<typename RealType, class F, std::size_t Order>
void
for_each_tape(F&& f, std::integral_constant<std::size_t, Order>)
{
  f(rdiff::get_active_tape<RealType, Order>());
  for_each_tape<RealType>(f, std::integral_constant<std::size_t, Order - 1>{});
}
/******************************************************************/
/**
 * @brief> function evaluation policy for reverse mode autodiff
 * @arg objective> objective function to evaluate
 * @arg x> argument list
 */
template<typename RealType, std::size_t Order = 1>
struct reverse_mode_function_eval_policy
{
  template<typename Objective, class ArgumentContainer>
  rdiff::rvar<RealType, Order> operator()(Objective&& objective,
                                          ArgumentContainer& x)
  {
    rdiff::get_active_tape<RealType, Order>().zero_grad();
    for_each_tape<RealType>([](auto& tape) { tape.rewind_to_last_checkpoint(); },
                            std::integral_constant<std::size_t, Order>{});

    return objective(x);
  }
};
/******************************************************************/

/**
 * @brief> gradient evaluation policy
 * @arg obj_f> objective
 * @arg x> argument list
 * @arg f_eval_olicy> funciton evaluation policy. These need to be
 *                    done in tandem
 * @arg obj_v> reference to variable inside gradient class
 */
template<typename RealType, std::size_t Order = 1>
struct reverse_mode_gradient_evaluation_policy
{
  template<class Objective,
           class ArgumentContainer,
           class FunctionEvaluationPolicy>
  void operator()(Objective&& obj_f,
                  ArgumentContainer& x,
                  FunctionEvaluationPolicy&& f_eval_pol,
                  RealType& obj_v,
                  std::vector<RealType>& g)
  {
    // compute objective via eval policy that takes care of tape
    rdiff::rvar<RealType, Order> v = f_eval_pol(obj_f, x);
    // For Order > 1 the adjoints of x are rvars on the rewound lower order
    // tape: give them fresh nodes before accumulating into them.
    for (auto& xi : x) {
      xi.adjoint() = typename std::decay<decltype(xi.adjoint())>::type(RealType(0));
    }
    v.backward();
    obj_v = v.item();
    g.resize(x.size());

    // copy gradients into gradient vector
    for (size_t i = 0; i < x.size(); ++i) {
      g[i] = static_cast<RealType>(x[i].adjoint());
    }
  }
};

/**
 * @brief> Hessian-vector product policy for reverse mode autodiff
 *
 * Computes Hv = H(x) v at the x of the last gradient evaluation by
 * reverse_mode_gradient_evaluation_policy<RealType, 2>: one reverse sweep of
 * the first order tape, which records the gradient, seeded with v.
 */
template<typename RealType>
struct reverse_mode_hessian_vector_product_policy
{
  template<class Objective, class ArgumentContainer>
  void operator()(Objective&& /*obj_f*/,
                  ArgumentContainer& x,
                  const std::vector<RealType>& v,
                  std::vector<RealType>& Hv)
  {
    // Only the nodes recorded since the evaluation, and the arguments,
    // take part in the sweep.
    auto& tape = rdiff::get_active_tape<RealType, 1>();
    for (auto it = tape.last_checkpoint(); it != tape.end(); ++it) {
      *it->get_adjoint_ptr() = RealType(0);
    }
    for (size_t i = 0; i < x.size(); ++i) {
      x[i].get_value().adjoint() = RealType(0);
      x[i].adjoint().adjoint() += v[i];
    }
    for (auto it = tape.end(), first = tape.last_checkpoint(); it != first;) {
      (--it)->backward();
    }
    Hv.resize(x.size());
    for (size_t i = 0; i < x.size(); ++i) {
      Hv[i] = x[i].get_value().adjoint();
    }
  }
};
/******************************************************************
 * init policies
 */
template<typename RealType, std::size_t Order = 1>
struct tape_initializer_rvar
{
    template<class ArgumentContainer>
    void operator()(ArgumentContainer& /*x*/) const noexcept
    {
        static_assert(std::is_same<typename ArgumentContainer::value_type,
                                   rdiff::rvar<RealType, Order>>::value,
                      "ArgumentContainer::value_type must be rdiff::rvar<RealType,Order>");
        for_each_tape<RealType>([](auto& tape) { tape.add_checkpoint(); },
                                std::integral_constant<std::size_t, Order>{});
    }
};

template<typename RealType>
struct random_uniform_initializer_rvar
{
  RealType low_, high_;
  size_t seed_;
  random_uniform_initializer_rvar(RealType low = 0,
                                  RealType high = 1,
                                  size_t seed = std::random_device{}())
    : low_(low)
    , high_(high)
    , seed_(seed) {};

  template<class ArgumentContainer>
  void operator()(ArgumentContainer& x) const
  {
    static std::mt19937 gen{ std::random_device{}() };
    static std::uniform_real_distribution<double> dist(0.0, 1.0);
    for (auto& xi : x) {
      xi = rdiff::rvar<RealType, 1>(static_cast<RealType>(dist(gen)));
    }
    auto& tape = rdiff::get_active_tape<RealType, 1>();
    tape.add_checkpoint();
  }
};

template<typename RealType>
struct costant_initializer_rvar
{
  RealType constant;
  explicit costant_initializer_rvar(RealType v = 0)
    : constant(v) {};
  template<class ArgumentContainer>
  void operator()(ArgumentContainer& x) const
  {
    for (auto& xi : x) {
      xi = rdiff::rvar<RealType, 1>(constant);
    }
    auto& tape = rdiff::get_active_tape<RealType, 1>();
    tape.add_checkpoint();
  }
};
} // namespace optimization
BOOST_MATH_NAMESPACE_END

#endif
