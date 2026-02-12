//           Copyright Maksym Zhelyenzyakov 2025-2026.
// Distributed under the Boost Software License, Version 1.0.
//      (See accompanying file LICENSE_1_0.txt or copy at
//           https://www.boost.org/LICENSE_1_0.txt)
#include "test_autodiff_reverse.hpp"
#include "test_functions_for_optimization.hpp"
#include <boost/math/differentiation/autodiff_reverse.hpp>
#include <boost/math/optimization/lbfgs.hpp>
#include <boost/math/optimization/minimizer.hpp>
#include <boost/math/optimization/gradient_descent.hpp>

namespace rdiff = boost::math::differentiation::reverse_mode;
namespace bopt = boost::math::optimization;

BOOST_AUTO_TEST_SUITE(basic_lbfgs)

BOOST_AUTO_TEST_CASE(default_lbfgs_test) //, T, all_float_types)
{
    using T                = double;
    constexpr size_t NITER = 10;
    constexpr size_t M     = 10;
    const T          eps   = T{1e-5};

    RandomSample<T>                  rng{T(-10), T(10)};
    std::array<rdiff::rvar<T, 1>, 2> x;
    x[0] = rng.next();
    x[1] = rng.next();

    auto opt = bopt::make_lbfgs(&rosenbrock_saddle<rdiff::rvar<T, 1>>, x, M);
    auto constraint = bopt::unconstrained_policy<std::array<rdiff::rvar<T, 1>,2>>{};
    auto convergence_policy = bopt::gradient_norm_convergence_policy<T>(T{ 1e-20 });

    auto result = bopt::minimize(opt, constraint, convergence_policy);
    for (auto& xi : x) {
        BOOST_REQUIRE_CLOSE(xi, T{1.0}, eps);
    }
}

// Custom initialization policy that zeros out the parameters
template<typename RealType>
struct zero_init_policy
{
    void operator()(std::vector<RealType>& x) const noexcept
    {
        std::fill(x.begin(), x.end(), RealType{0});
    }
};

template<typename RealType>
struct analytic_objective_eval_pol
{
    template<typename Objective, typename ArgumentContainer>
    RealType operator()(Objective&& objective, ArgumentContainer& x)
    {
        return objective(x);
    }
};
template<typename RealType>
struct analytic_gradient_eval_pol
{
    template<class Objective, class ArgumentContainer, class FunctionEvaluationPolicy>
    void operator()(Objective&&                obj_f,
                    ArgumentContainer&         x,
                    FunctionEvaluationPolicy&& f_eval_pol,
                    RealType&                  obj_v,
                    std::vector<RealType>&     grad_container)
    {
        RealType v = f_eval_pol(obj_f, x);
        obj_v      = v;
        grad_container.resize(x.size());
        for (size_t i = 0; i < x.size(); ++i) {
            grad_container[i] = 2 * x[i];
        }
    }
};

BOOST_AUTO_TEST_CASE_TEMPLATE(custom_init_lbfgs_test, T, all_float_types)
{
    constexpr size_t M   = 8;
    const T          eps = T{1e-6};

    RandomSample<T>                  rng{T(-5), T(5)};
    std::array<rdiff::rvar<T, 1>, 2> x;
    x[0] = rng.next();
    x[1] = rng.next();

    auto opt    = bopt::make_lbfgs(&rosenbrock_saddle<rdiff::rvar<T, 1>>,
                                x,
                                M,
                                bopt::costant_initializer_rvar<T>(0.0));
    auto constraint = bopt::unconstrained_policy<std::array<rdiff::rvar<T, 1>,2>>{};
    auto convergence_policy = bopt::gradient_norm_convergence_policy<T>(T{ 1e-8 });
    auto result = minimize(opt, constraint, convergence_policy);

    for (auto& xi : x) {
        BOOST_REQUIRE_CLOSE(xi, T{1.0}, eps);
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(analytic_lbfgs_test, T, all_float_types)
{
    constexpr size_t M   = 10;
    const T          eps = T{1e-3};

    RandomSample<T> rng{T(-5), T(5)};
    std::vector<T>  x(3);
    for (auto& xi : x)
        xi = rng.next();

    auto opt = bopt::make_lbfgs(&quadratic<T>,                    // Objective
                                x,                                // Arguments
                                M,                                // History size
                                zero_init_policy<T>{},            // Initialization
                                analytic_objective_eval_pol<T>{}, // Function eval
                                analytic_gradient_eval_pol<T>{},  // Gradient eval
                                bopt::armijo_line_search_policy<T>{});

    auto result = minimize(opt);

    for (auto& xi : x) {
        BOOST_REQUIRE_SMALL(xi, eps);
    }
}
BOOST_AUTO_TEST_SUITE_END()
