#ifndef GRADIENT_DESCENT_HPP
#define GRADIENT_DESCENT_HPP
#include <boost/math/optimization/detail/differentiable_opt_utilties.hpp>
#include <boost/math/optimization/detail/gradient_opt_base.hpp>
#include <boost/math/optimization/detail/rdiff_optimization_policies.hpp>
#include <functional>
#include <vector>
namespace rdiff = boost::math::differentiation::reverse_mode;

namespace boost {
namespace math {
namespace optimization {

template<typename RealType>
struct gradient_descent_update_policy
{
    RealType lr_;
    gradient_descent_update_policy(RealType lr)
        : lr_(lr){};

    template<typename ArgumentType,
             typename = typename std::enable_if<boost::math::differentiation::reverse_mode::detail::
                                                    is_expression<ArgumentType>::value>::type>
    void operator()(ArgumentType &x, RealType &g)
    {
        // this update effectively "mutes" the tape
        // TODO: add a tape scope guard method so that
        // you can do math on autodiff types without
        // accumulating gradients
        x.get_value() -= lr_ * g;
    }
    template<
        typename ArgumentType,
        typename std::enable_if<
            !boost::math::differentiation::reverse_mode::detail::is_expression<ArgumentType>::value,
            int>::type
        = 0>
    void operator()(ArgumentType &x, RealType &g) const
    {
        x -= lr_ * g;
    }
};
template<typename ArgumentContainer,
         typename RealType,
         class Objective,
         class InitializationPolicy,
         class ObjectiveEvalPolicy,
         class GradEvalPolicy>
class gradient_descent : public abstract_optimizer<ArgumentContainer,
                                                   RealType,
                                                   Objective,
                                                   InitializationPolicy,
                                                   ObjectiveEvalPolicy,
                                                   GradEvalPolicy,
                                                   gradient_descent_update_policy<RealType>,
                                                   gradient_descent<ArgumentContainer,
                                                                    RealType,
                                                                    Objective,
                                                                    InitializationPolicy,
                                                                    ObjectiveEvalPolicy,
                                                                    GradEvalPolicy>>
{
    using base_opt = abstract_optimizer<ArgumentContainer,
                                        RealType,
                                        Objective,
                                        InitializationPolicy,
                                        ObjectiveEvalPolicy,
                                        GradEvalPolicy,
                                        gradient_descent_update_policy<RealType>,
                                        gradient_descent<ArgumentContainer,
                                                         RealType,
                                                         Objective,
                                                         InitializationPolicy,
                                                         ObjectiveEvalPolicy,
                                                         GradEvalPolicy>>;

public:
    using base_opt::base_opt;
};
template<class Objective, typename ArgumentContainer, typename RealType>
auto make_gradient_descent(Objective &&obj, ArgumentContainer &x, RealType lr = RealType{0.01})
{
    return gradient_descent<ArgumentContainer,
                            RealType,
                            Objective,
                            tape_initializer_rvar<RealType>,
                            reverse_mode_function_eval_policy<RealType>,
                            reverse_mode_gradient_evaluation_policy<RealType>>(
        std::forward<Objective>(obj),
        x,
        tape_initializer_rvar<RealType>{},
        reverse_mode_function_eval_policy<RealType>{},
        reverse_mode_gradient_evaluation_policy<RealType>{},
        gradient_descent_update_policy<RealType>(lr));
}
template<class Objective, typename ArgumentContainer, typename RealType, class InitializationPolicy>
auto make_gradient_descent(Objective            &&obj,
                           ArgumentContainer     &x,
                           RealType               lr,
                           InitializationPolicy &&ip)
{
    return gradient_descent<ArgumentContainer,
                            RealType,
                            Objective,
                            InitializationPolicy,
                            reverse_mode_function_eval_policy<RealType>,
                            reverse_mode_gradient_evaluation_policy<RealType>>(
        std::forward<Objective>(obj),
        x,
        std::forward<InitializationPolicy>(ip),
        reverse_mode_function_eval_policy<RealType>{},
        reverse_mode_gradient_evaluation_policy<RealType>{},
        gradient_descent_update_policy<RealType>(lr));
}
template<typename ArgumentContainer,
         typename RealType,
         class Objective,
         class InitializationPolicy,
         class ObjectiveEvalPolicy,
         class GradEvalPolicy>
auto make_gradient_descent(Objective            &&obj,
                           ArgumentContainer     &x,
                           RealType              &lr,
                           InitializationPolicy &&ip,
                           ObjectiveEvalPolicy  &&oep,
                           GradEvalPolicy       &&gep)
{
    return gradient_descent<ArgumentContainer,
                            RealType,
                            Objective,
                            InitializationPolicy,
                            ObjectiveEvalPolicy,
                            GradEvalPolicy>(std::forward<Objective>(obj),
                                            x,
                                            std::forward<InitializationPolicy>(ip),
                                            std::forward<ObjectiveEvalPolicy>(oep),
                                            std::forward<GradEvalPolicy>(gep),
                                            gradient_descent_update_policy<RealType>{lr});
}
} // namespace optimization
} // namespace math
} // namespace boost
#endif
