//  Copyright 2026 Nicholas Thompson and Matt Borland
//  Distributed under the Boost Software License, Version 1.0.

#include <boost/math/quadrature/van_den_bos_unit_square.hpp>
#include "math_unit_test.hpp"

#include <array>
#include <cmath>
#include <complex>
#include <cstddef>
#include <thread>
#include <vector>

namespace {

using point = std::array<double, 2>;
using boost::math::quadrature::van_den_bos_unit_square;

std::complex<double> exact_complex_exp_5ixy()
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
        {
            break;
        }
    }

    return exact;
}

double shifted_legendre(unsigned n, double x)
{
    const double t = 2 * x - 1;
    if (n == 0)
    {
        return 1;
    }
    if (n == 1)
    {
        return t;
    }

    double pnm2 = 1;
    double pnm1 = t;
    for (unsigned k = 2; k <= n; ++k)
    {
        const double pn =
            ((2 * k - 1) * t * pnm1 - (k - 1) * pnm2) / k;
        pnm2 = pnm1;
        pnm1 = pn;
    }
    return pnm1;
}

void test_table_properties()
{
    namespace detail = boost::math::quadrature::detail;

    std::size_t previous_size = 0;
    for (std::size_t level_index = 0;
         level_index <= detail::van_den_bos_unit_square_max_level;
         ++level_index)
    {
        const detail::van_den_bos_unit_square_level_view level =
            detail::van_den_bos_unit_square_level(level_index);
        const std::size_t moment_count =
            (level.degree + 1) * (level.degree + 2) / 2;

        const std::size_t expected_size = level_index == 0
            ? (level.order + 1) * (level.order + 1)
            : previous_size + moment_count;
        CHECK_EQUAL(level.size, expected_size);

        double weight_sum = 0;
        for (std::size_t j = 0; j < level.size; ++j)
        {
            CHECK_EQUAL(level.weights[j] > 0, true);
            CHECK_GE(level.points[j][0], 0.0);
            CHECK_LE(level.points[j][0], 1.0);
            CHECK_GE(level.points[j][1], 0.0);
            CHECK_LE(level.points[j][1], 1.0);
            weight_sum += level.weights[j];

            for (std::size_t i = 0; i < j; ++i)
            {
                CHECK_EQUAL(level.points[i] == level.points[j], false);
            }
        }
        CHECK_ABSOLUTE_ERROR(weight_sum, 1.0, 4e-15);

        if (level_index > 0)
        {
            const detail::van_den_bos_unit_square_level_view previous =
                detail::van_den_bos_unit_square_level(level_index - 1);
            CHECK_EQUAL(previous.size, previous_size);
            for (std::size_t j = 0; j < previous.size; ++j)
            {
                CHECK_EQUAL(level.points[j] == previous.points[j], true);
            }
        }

        for (unsigned total = 0; total <= level.degree; ++total)
        {
            for (unsigned px = 0; px <= total; ++px)
            {
                const unsigned py = total - px;
                double moment = 0;
                for (std::size_t j = 0; j < level.size; ++j)
                {
                    moment += level.weights[j] *
                        shifted_legendre(px, level.points[j][0]) *
                        shifted_legendre(py, level.points[j][1]);
                }
                const double exact = total == 0 ? 1.0 : 0.0;
                CHECK_ABSOLUTE_ERROR(moment, exact, 5e-13);
            }
        }

        previous_size = level.size;
    }
}

void test_polynomial()
{
    double error = 0;
    std::size_t evaluations = 0;

    const double value = van_den_bos_unit_square(
        [](point const& p)
        {
            return p[0] * p[0] + p[1] * p[1] + p[0] * p[1];
        },
        1e-12,
        5,
        &error,
        nullptr,
        &evaluations);

    CHECK_ABSOLUTE_ERROR(value, 11.0 / 12.0, 2e-14);
    CHECK_GE(evaluations, std::size_t(1));
}

void test_exp_x_plus_y()
{
    const double exact = std::pow(std::expm1(1.0), 2);

    double error = 0;
    const double value = van_den_bos_unit_square(
        [](point const& p)
        {
            return std::exp(p[0] + p[1]);
        },
        1e-10,
        5,
        &error);

    CHECK_ABSOLUTE_ERROR(value, exact, 2e-10);
    CHECK_LE(error, 2e-10);
}

void test_complex_integrand()
{
    const auto exact = exact_complex_exp_5ixy();

    double error = 0;
    double L1 = 0;
    std::size_t evaluations = 0;

    const auto computed = van_den_bos_unit_square(
        [](point const& p)
        {
            return std::exp(
                std::complex<double>(0, 5 * p[0] * p[1]));
        },
        1e-8,
        5,
        &error,
        &L1,
        &evaluations);

    CHECK_ABSOLUTE_ERROR(exact.real(), computed.real(), 2e-8);
    CHECK_ABSOLUTE_ERROR(exact.imag(), computed.imag(), 5e-8);
    CHECK_LE(error, 1e-8);
    CHECK_GE(L1, 0.0);
    CHECK_GE(evaluations, std::size_t(1));
}

void test_max_level()
{
    double error = 0;

    const double value = van_den_bos_unit_square(
        [](point const& p)
        {
            return std::sqrt(p[0] + p[1]);
        },
        1e-30,
        5,
        &error);

    const double exact =
        (4.0 / 15.0) * (std::pow(2.0, 2.5) - 2.0);

    CHECK_ABSOLUTE_ERROR(value, exact, 1e-4);
    CHECK_GE(error, 0.0);
}

void test_concurrent_integrations()
{
    constexpr std::size_t thread_count = 4;

    std::array<double, thread_count> results{};
    std::array<double, thread_count> errors{};
    std::array<std::size_t, thread_count> evaluations{};
    std::vector<std::thread> threads;
    threads.reserve(thread_count);

    for (std::size_t k = 0; k < thread_count; ++k)
    {
        threads.emplace_back([&, k]
        {
            const double scale = double(k + 1);

            results[k] = van_den_bos_unit_square(
                [=](point const& p)
                {
                    return std::exp(
                        scale * (p[0] + p[1]) / 8.0);
                },
                1e-9,
                5,
                &errors[k],
                nullptr,
                &evaluations[k]);
        });
    }

    for (auto& thread : threads)
    {
        thread.join();
    }

    for (std::size_t k = 0; k < thread_count; ++k)
    {
        const double a = double(k + 1) / 8.0;
        const double exact =
            std::pow(std::expm1(a) / a, 2);

        CHECK_ABSOLUTE_ERROR(exact, results[k], 2e-8);
        CHECK_LE(
            errors[k],
            2e-9 * (std::max)(1.0, std::abs(results[k])));
        CHECK_GE(evaluations[k], std::size_t(1));
    }
}

} // namespace

int main()
{
    test_table_properties();
    test_polynomial();
    test_exp_x_plus_y();
    test_complex_integrand();
    test_max_level();
    test_concurrent_integrations();

    return boost::math::test::report_errors();
}
