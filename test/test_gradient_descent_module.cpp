//           Copyright Maksym Zhelyenzyakov 2025-2026.
//           Copyright Matt Borland 2026.
// Distributed under the Boost Software License, Version 1.0.
//      (See accompanying file LICENSE_1_0.txt or copy at
//           https://www.boost.org/LICENSE_1_0.txt)
//
// Module build check for boost::math::optimization::gradient_descent. When
// BOOST_MATH_BUILD_MODULE is defined the component is consumed through
// `import boost.math;` and the standard library through `import std;` (performed
// by math_unit_test.hpp). Nothing standard or Boost is included textually in
// that mode.
//
// The optimizer minimizes f(x) = sum_i x_i^2, whose gradient is 2 x_i. With a
// fixed learning rate lr a single descent step maps each coordinate x_i to
// x_i * (1 - 2*lr), and repeated steps drive the minimizer to the origin. Both
// the one-step and the limiting values are mathematically certain; the setup
// mirrors the analytic_derivative_policies case in
// test_gradient_descent_optimizer.cpp.

#ifndef BOOST_MATH_BUILD_MODULE
#include <boost/math/optimization/gradient_descent.hpp>
#else
import boost.math;
#endif

#ifndef BOOST_MATH_BUILD_MODULE
#include <vector>
#include <cstddef>
#endif
#include "math_unit_test.hpp"

namespace bopt = boost::math::optimization;

// f(x) = sum_i x_i^2. Its gradient, 2 x_i, is supplied analytically by the
// gradient policy below, so one descent step scales each x_i by (1 - 2*lr).
template <typename Real>
Real sum_of_squares(const std::vector<Real>& x)
{
    Real s {0};
    for (const auto& xi : x)
    {
        s += xi * xi;
    }
    return s;
}

// Leave the caller-supplied starting point untouched.
template <typename Real>
struct no_init_policy
{
    void operator()(std::vector<Real>&) const noexcept {}
};

// Evaluate the objective directly.
template <typename Real>
struct analytic_objective_eval_pol
{
    template <typename Objective, typename ArgumentContainer>
    Real operator()(Objective&& objective, ArgumentContainer& x)
    {
        return objective(x);
    }
};

// Supply the analytic gradient of sum_of_squares: grad_i = 2 x_i.
template <typename Real>
struct analytic_gradient_eval_pol
{
    template <class Objective, class ArgumentContainer, class FunctionEvaluationPolicy>
    void operator()(Objective&& obj_f,
                    ArgumentContainer& x,
                    FunctionEvaluationPolicy&& f_eval_pol,
                    Real& obj_v,
                    std::vector<Real>& grad_container)
    {
        obj_v = f_eval_pol(obj_f, x);
        grad_container.resize(x.size());
        for (std::size_t i = 0; i < x.size(); ++i)
        {
            grad_container[i] = 2 * x[i];
        }
    }
};

template <typename Real>
void test_gradient_descent(Real, const char* type_name)
{
    std::cout << "Testing gradient_descent for type " << type_name << std::endl;

    // One step with lr = 0.01 scales each coordinate by (1 - 2*lr) = 0.98:
    // {3, -4, 5} -> {2.94, -3.92, 4.9}.
    {
        Real lr = Real {1} / Real {100};
        std::vector<Real> x = {Real {3}, Real {-4}, Real {5}};

        auto gdopt = bopt::make_gradient_descent(&sum_of_squares<Real>,
                                                 x,
                                                 lr,
                                                 no_init_policy<Real> {},
                                                 analytic_objective_eval_pol<Real> {},
                                                 analytic_gradient_eval_pol<Real> {});
        gdopt.step();

        CHECK_ULP_CLOSE(2.94, x[0], 4);
        CHECK_ULP_CLOSE(-3.92, x[1], 4);
        CHECK_ULP_CLOSE(4.9, x[2], 4);
    }

    // Repeated steps drive the minimizer of f(x) = sum_i x_i^2 to the origin.
    {
        Real lr = Real {1} / Real {100};
        std::vector<Real> x = {Real {3}, Real {-4}, Real {5}};

        auto gdopt = bopt::make_gradient_descent(&sum_of_squares<Real>,
                                                 x,
                                                 lr,
                                                 no_init_policy<Real> {},
                                                 analytic_objective_eval_pol<Real> {},
                                                 analytic_gradient_eval_pol<Real> {});
        for (std::size_t i = 0; i < 2000; ++i)
        {
            gdopt.step();
        }

        CHECK_ABSOLUTE_ERROR(Real(0), x[0], Real(1e-3));
        CHECK_ABSOLUTE_ERROR(Real(0), x[1], Real(1e-3));
        CHECK_ABSOLUTE_ERROR(Real(0), x[2], Real(1e-3));
    }
}

int main()
{
    test_gradient_descent(0.0, "double");
    test_gradient_descent(0.0F, "float");

    return boost::math::test::report_errors();
}
