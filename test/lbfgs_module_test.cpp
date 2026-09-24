//           Copyright Maksym Zhelyenzyakov 2025-2026.
//           Copyright Matt Borland 2026.
// Distributed under the Boost Software License, Version 1.0.
//      (See accompanying file LICENSE_1_0.txt or copy at
//           https://www.boost.org/LICENSE_1_0.txt)
//
// Module build check for boost::math::optimization L-BFGS. It exercises the
// exported make_lbfgs / minimize entry points when built against the boost.math
// module. The reverse-mode-autodiff paths of test_lbfgs.cpp are not usable here
// (autodiff_reverse is not part of the module export surface), so this test
// drives the fully analytic path: minimizing the isotropic quadratic
// f(x) = sum x_i^2 whose minimizer is the origin with value 0.

#ifndef BOOST_MATH_BUILD_MODULE
#include <boost/math/optimization/lbfgs.hpp>
#include <boost/math/optimization/minimizer.hpp>
#else
import boost.math;
#endif

#ifndef BOOST_MATH_BUILD_MODULE
#include <vector>
#include <cstddef>
#endif
#include "math_unit_test.hpp"

namespace bopt = boost::math::optimization;

// Objective: isotropic quadratic, minimum 0 at the origin.
template<typename Real>
Real quadratic(std::vector<Real>& x)
{
    Real res{ 0 };
    for (auto& item : x)
    {
        res += item * item;
    }
    return res;
}

// No-op initializer: keep the caller-supplied starting point.
template<typename Real>
struct identity_init_policy
{
    void operator()(std::vector<Real>&) const noexcept {}
};

// Objective evaluation policy: just call the objective.
template<typename Real>
struct analytic_objective_eval_policy
{
    template<typename Objective, typename ArgumentContainer>
    Real operator()(Objective&& objective, ArgumentContainer& x)
    {
        return objective(x);
    }
};

// Gradient evaluation policy: analytic gradient of sum x_i^2 is 2*x.
template<typename Real>
struct analytic_gradient_eval_policy
{
    template<class Objective, class ArgumentContainer, class FunctionEvaluationPolicy>
    void operator()(Objective&&                obj_f,
                    ArgumentContainer&         x,
                    FunctionEvaluationPolicy&& f_eval_pol,
                    Real&                      obj_v,
                    std::vector<Real>&         grad_container)
    {
        obj_v = f_eval_pol(obj_f, x);
        grad_container.resize(x.size());
        for (std::size_t i{ 0 }; i < x.size(); ++i)
        {
            grad_container[i] = 2 * x[i];
        }
    }
};

// Self-contained backtracking Armijo line search. The library's own
// armijo_line_search_policy and its dot/axpy helpers are not exported by the
// module, so the search is reimplemented here with plain loops.
template<typename Real>
struct armijo_line_search_policy
{
    template<class Objective, class ObjectiveEvalPolicy, class GradientEvalPolicy, class ArgumentContainer>
    Real operator()(Objective&                   objective,
                    ObjectiveEvalPolicy&         obj_eval,
                    GradientEvalPolicy& /*grad_eval*/,
                    ArgumentContainer&           x,
                    const std::vector<Real>&     g,
                    const std::vector<Real>&     p,
                    Real                         f_x) const
    {
        const Real c{ static_cast<Real>(1e-4) };
        const Real rho{ static_cast<Real>(0.5) };

        Real gTp{ 0 };
        for (std::size_t i{ 0 }; i < g.size(); ++i)
        {
            gTp += g[i] * p[i];
        }

        Real alpha{ 1 };
        for (int iter{ 0 }; iter < 20; ++iter)
        {
            ArgumentContainer x_trial{ x };
            for (std::size_t i{ 0 }; i < x_trial.size(); ++i)
            {
                x_trial[i] += alpha * p[i];
            }
            const Real f_trial{ obj_eval(objective, x_trial) };
            if (f_trial <= f_x + c * alpha * gTp)
            {
                return alpha;
            }
            alpha *= rho;
        }
        return alpha;
    }
};

template<typename Real>
void test_quadratic_minimization()
{
    constexpr std::size_t M{ 10 };
    const Real            tol{ static_cast<Real>(1e-3) };

    // Non-trivial, distinct starting coordinates.
    std::vector<Real> x{ static_cast<Real>(2), static_cast<Real>(-1.5), static_cast<Real>(3) };

    auto opt = bopt::make_lbfgs(&quadratic<Real>,                     // Objective
                                x,                                    // Arguments
                                M,                                    // History size
                                identity_init_policy<Real>{},         // Initialization
                                analytic_objective_eval_policy<Real>{}, // Function eval
                                analytic_gradient_eval_policy<Real>{},  // Gradient eval
                                armijo_line_search_policy<Real>{});   // Line search

    auto result = bopt::minimize(opt);

    // The minimizer of sum x_i^2 is the origin; each coordinate converges to 0.
    for (auto& xi : x)
    {
        CHECK_ABSOLUTE_ERROR(static_cast<Real>(0), xi, tol);
    }

    // Objective value at the minimum is 0, and the run reports convergence.
    CHECK_ABSOLUTE_ERROR(static_cast<Real>(0), result.objective_value, tol);
    CHECK_EQUAL(true, result.converged);
}

int main()
{
    test_quadratic_minimization<double>();
    test_quadratic_minimization<float>();

    return boost::math::test::report_errors();
}
