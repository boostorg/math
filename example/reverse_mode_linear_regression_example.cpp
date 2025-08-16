//           Copyright Maksym Zhelyenzyakov 2025-2026.
// Distributed under the Boost Software License, Version 1.0.
//      (See accompanying file LICENSE_1_0.txt or copy at
//           https://www.boost.org/LICENSE_1_0.txt)
#include <array>
#include <boost/math/differentiation/autodiff_reverse.hpp>
#include <iostream>
#include <random>
#include <vector>

using namespace boost::math::differentiation::reverse_mode;
double random_double(double min_val, double max_val)
{
    static std::random_device              rd;
    static std::mt19937                    gen(rd());
    std::uniform_real_distribution<double> dist(min_val, max_val);
    return dist(gen);
}

template<size_t N>
rvar<double, 1> loss(std::array<double, N>& y_target, std::array<rvar<double, 1>, N>& y_fit)
{
    rvar<double, 1> loss_v = make_rvar<double, 1>(0.0);
    for (size_t i = 0; i < N; ++i) {
        loss_v += pow(abs(y_target[i] - y_fit[i]), 2) / N;
    }
    return loss_v;
}
double noisy_linear_function(double intercept, double slope, double x)
{
    return intercept + slope * x + random_double(-0.1, 0.1);
}

template<size_t N>
std::array<rvar<double, 1>, N> model(rvar<double, 1>&       a,
                                     rvar<double, 1>&       b,
                                     std::array<double, N>& x)
{
    std::array<rvar<double, 1>, N> ret;
    for (size_t i = 0; i < N; ++i) {
        ret[i] = a * x[i] + b;
    }
    return ret;
}
int main()
{
    double                               slope            = random_double(-5, 5);
    double                               intercept        = random_double(-5, 5);

    const size_t                         num_data_samples = 100;
    /**/
    std::array<double, num_data_samples> noisy_data_x;
    std::array<double, num_data_samples> noisy_data_y;
    for (size_t i = 0; i < num_data_samples; i++) {
        double x        = random_double(-1, 1);
        double y        = noisy_linear_function(intercept, slope, x);
        noisy_data_x[i] = x;
        noisy_data_y[i] = y;
    }

    double                    slope_guess     = random_double(-5, 5);
    double                    intercept_guess = random_double(-5, 5);

    rvar<double, 1>           a               = make_rvar<double, 1>(slope_guess);
    rvar<double, 1>           b               = make_rvar<double, 1>(intercept_guess);

    gradient_tape<double, 1>& tape            = get_active_tape<double, 1>();
    tape.add_checkpoint();

    auto            y_fit         = model(a, b, noisy_data_x);
    rvar<double, 1> loss_v        = loss(noisy_data_y, y_fit);

    double          learning_rate = 1e-3;
    while (loss_v > 0.005) {
        tape.rewind_to_last_checkpoint();
        y_fit    = model(a, b, noisy_data_x);
        loss_v   = loss(noisy_data_y, y_fit);
        auto gv  = grad(loss_v, &a, &b);
        a       -= gv[0] * learning_rate;
        b       -= gv[1] * learning_rate;
    }

    double slope_error              = std::abs(slope - a.item());
    double intercept_error          = std::abs(intercept - b.item());
    double relative_slope_error     = slope_error / std::abs(slope);
    double relative_intercept_error = intercept_error / std::abs(intercept);

    std::cout << "Autodiff Linear Regression Summary \n";
    std::cout << "learning rate : " << learning_rate << "\n";
    std::cout << "true slope: " << slope;
    std::cout << " regression: " << a.item() << "\n";

    std::cout << "relative error (slope): " << relative_slope_error << "\n";
    std::cout << "absolute error (slope): " << slope_error << "\n";
    std::cout << "true intercept: " << intercept;
    std::cout << " regression: " << b.item() << "\n";
    std::cout << "absolute error (intercept): " << intercept_error << "\n";
    std::cout << "aelative error (intercept): " << relative_intercept_error << "\n";
    std::cout << "-------------------------------" << std::endl;
}
