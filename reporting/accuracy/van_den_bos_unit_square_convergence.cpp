/*
 * Convergence trace driver for van_den_bos_unit_square.
 *
 * Build:
 *   c++ -O3 -DNDEBUG -DBOOST_MATH_INSTRUMENT -std=c++17 -Iinclude \
 *     reporting/accuracy/van_den_bos_unit_square_convergence.cpp \
 *     -o vdb_unit_square_convergence
 *
 * Run:
 *   ./vdb_unit_square_convergence > van_den_bos_unit_square_convergence.csv
 */

#include <boost/math/quadrature/van_den_bos_unit_square.hpp>

#include <array>
#include <cmath>
#include <complex>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <limits>

namespace {

using point2 = std::array<double, 2>;
using boost::math::quadrature::van_den_bos_unit_square;

template <class F, class Exact>
void run_case(char const* name, F const& f, Exact const& exact, double tol = 1e-12)
{
    double error = 0;
    double l1 = 0;
    std::size_t evaluations = 0;

    std::cout << std::setprecision(std::numeric_limits<double>::max_digits10)
              << "VDB_CASE,name=" << name << ",exact=" << exact << '\n';
    (void)van_den_bos_unit_square(
        f, tol, 5, &error, &l1, &evaluations);
    std::cout << "VDB_END,name=" << name
              << ",error_estimate=" << error
              << ",evaluations=" << evaluations << '\n';
}

std::complex<double> exact_exp_5ixy()
{
    std::complex<double> exact(0, 0);
    std::complex<double> power(1, 0);
    double factorial = 1;

    for (unsigned k = 0; k < 100; ++k)
    {
        if (k > 0)
        {
            power *= std::complex<double>(0, 5);
            factorial *= double(k);
        }
        const double kp1 = double(k + 1);
        const auto term = power / (factorial * kp1 * kp1);
        exact += term;
        if (std::abs(term) < 1e-18)
            break;
    }
    return exact;
}

} // namespace

int main()
{
    run_case(
        "quadratic",
        [](point2 const& p)
        {
            return p[0] * p[0] + p[1] * p[1] + p[0] * p[1];
        },
        11.0 / 12.0);

    run_case(
        "exp_x_plus_y",
        [](point2 const& p)
        {
            return std::exp(p[0] + p[1]);
        },
        std::pow(std::expm1(1.0), 2));

    run_case(
        "exp_xy",
        [](point2 const& p)
        {
            return std::exp(p[0] * p[1]);
        },
        [] {
            double s = 0;
            double factorial = 1;
            for (unsigned k = 0; k < 100; ++k)
            {
                if (k > 0)
                    factorial *= double(k);
                const double kp1 = double(k + 1);
                const double term = 1.0 / (factorial * kp1 * kp1);
                s += term;
                if (term < 1e-18)
                    break;
            }
            return s;
        }());

    run_case(
        "one_over_1_plus_x_plus_y",
        [](point2 const& p)
        {
            return 1.0 / (1.0 + p[0] + p[1]);
        },
        3.0 * std::log(3.0) - 4.0 * std::log(2.0));

    run_case(
        "complex_exp_5ixy",
        [](point2 const& p)
        {
            return std::exp(std::complex<double>(0, 5 * p[0] * p[1]));
        },
        exact_exp_5ixy());

    run_case(
        "cos_5xy",
        [](point2 const& p)
        {
            return std::cos(5.0 * p[0] * p[1]);
        },
        exact_exp_5ixy().real());

    run_case(
        "sqrt_x_plus_y",
        [](point2 const& p)
        {
            return std::sqrt(p[0] + p[1]);
        },
        (4.0 / 15.0) * (std::pow(2.0, 2.5) - 2.0));

    return 0;
}
