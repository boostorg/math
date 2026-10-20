//  Copyright 2026 Nicholas Thompson and Matt Borland
//  Distributed under the Boost Software License, Version 1.0.
//  https://www.boost.org/LICENSE_1_0.txt

#include <boost/math/quadrature/van_den_bos_unit_square.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/math/special_functions/ellint_1.hpp>
#include "math_unit_test.hpp"

#include <array>
#include <cmath>
#include <complex>
#include <cstddef>
#include <limits>
#include <thread>
#include <type_traits>
#include <vector>

namespace {

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

std::vector<std::array<double, 2>> d4_orbit(std::array<double, 2> const& p)
{
    const std::array<std::array<double, 2>, 8> candidates{{
        {{p[0], p[1]}}, {{p[1], p[0]}},
        {{1 - p[0], p[1]}}, {{p[1], 1 - p[0]}},
        {{p[0], 1 - p[1]}}, {{1 - p[1], p[0]}},
        {{1 - p[0], 1 - p[1]}}, {{1 - p[1], 1 - p[0]}}}};
    std::vector<std::array<double, 2>> result;
    for (auto const& candidate : candidates)
    {
        bool duplicate = false;
        for (auto const& old : result)
            if (candidate == old) { duplicate = true; break; }
        if (!duplicate) result.push_back(candidate);
    }
    return result;
}

void test_table_properties()
{
    namespace detail = boost::math::quadrature::detail;
    constexpr std::array<unsigned, 7> exact_degrees{{4, 7, 11, 15, 21, 31, 45}};

    CHECK_EQUAL(detail::van_den_bos_unit_square_level_count,
        exact_degrees.size());

    std::size_t previous_size = 0;
    for (std::size_t level_index = 0;
         level_index < detail::van_den_bos_unit_square_level_count;
         ++level_index)
    {
        const detail::van_den_bos_unit_square_level_view level =
            detail::van_den_bos_unit_square_level(level_index);
        CHECK_EQUAL(level.size > previous_size, true);

        double weight_sum = 0;
        for (std::size_t j = 0; j < level.size; ++j)
        {
            CHECK_EQUAL(level.weights[j] > 0, true);
            const auto point = level.point(j);
            CHECK_EQUAL(point[0] > 0.0, true);
            CHECK_LE(point[0], 0.5);
            CHECK_EQUAL(point[1] > 0.0, true);
            CHECK_LE(point[1], point[0]);
            for (auto const& image : d4_orbit(point))
            {
                CHECK_EQUAL(image[0] > 0.0, true);
                CHECK_EQUAL(image[0] < 1.0, true);
                CHECK_EQUAL(image[1] > 0.0, true);
                CHECK_EQUAL(image[1] < 1.0, true);
            }
            weight_sum += level.weights[j] * d4_orbit(point).size();

            for (std::size_t i = 0; i < j; ++i)
            {
                CHECK_EQUAL(level.point(i) == point, false);
            }
        }
        CHECK_ABSOLUTE_ERROR(weight_sum, 1.0, 8e-15);

        if (level_index > 0)
        {
            const detail::van_den_bos_unit_square_level_view previous =
                detail::van_den_bos_unit_square_level(level_index - 1);
            CHECK_EQUAL(previous.size, previous_size);
            for (std::size_t j = 0; j < previous.size; ++j)
            {
                CHECK_EQUAL(level.point(j) == previous.point(j), true);
            }
        }

        for (unsigned total = 0; total <= exact_degrees[level_index]; ++total)
        {
            for (unsigned px = 0; px <= total; ++px)
            {
                const unsigned py = total - px;
                double moment = 0;
                for (std::size_t j = 0; j < level.size; ++j)
                {
                    for (auto const& point : d4_orbit(level.point(j)))
                        moment += level.weights[j] *
                            shifted_legendre(px, point[0]) *
                            shifted_legendre(py, point[1]);
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
        [](double x, double y)
        {
            return x * x + y * y + x * y;
        },
        1e-12,
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
        [](double x, double y)
        {
            return std::exp(x + y);
        },
        1e-10,
        &error);

    CHECK_ABSOLUTE_ERROR(value, exact, 2e-10);
    CHECK_LE(error, 2e-10);
}

void test_dilogarithm_integral()
{
    const double pi = boost::math::constants::pi<double>();
    const double log_two = std::log(2.0);
    const double exact = pi * pi / 6 - log_two * log_two;

    double error = 0;
    const double value = van_den_bos_unit_square(
        [](double x, double y)
        {
            return 1 / (1 - x * y / 2);
        },
        1e-10,
        &error);

    CHECK_ABSOLUTE_ERROR(value, exact, 2e-10);
    CHECK_LE(error, 2e-10);
}

void test_square_lattice_green_function()
{
    static constexpr double z = 3;
    const double pi = boost::math::constants::pi<double>();
    const double exact = 2 * boost::math::ellint_1(2 / z) / (pi * z);

    double error = 0;
    const double value = van_den_bos_unit_square(
        [](double x, double y)
        {
            const double two_pi =
                2 * boost::math::constants::pi<double>();
            return 1 / (z - std::cos(two_pi * x)
                          - std::cos(two_pi * y));
        },
        1e-6,
        &error);

    CHECK_ABSOLUTE_ERROR(value, exact, 5e-9);
    CHECK_LE(error, 1e-6);
}

void test_complex_integrand()
{
    const auto exact = exact_complex_exp_5ixy();

    double error = 0;
    double L1 = 0;
    std::size_t evaluations = 0;

    const auto computed = van_den_bos_unit_square(
        [](double x, double y)
        {
            return std::exp(
                std::complex<double>(0, 5 * x * y));
        },
        1e-8,
        &error,
        &L1,
        &evaluations);

    CHECK_ABSOLUTE_ERROR(exact.real(), computed.real(), 2e-8);
    CHECK_ABSOLUTE_ERROR(exact.imag(), computed.imag(), 5e-8);
    CHECK_LE(error, 1e-8);
    CHECK_GE(L1, 0.0);
    CHECK_GE(evaluations, std::size_t(1));
}

void test_final_precomputed_level()
{
    double error = 0;

    const double value = van_den_bos_unit_square(
        [](double x, double y)
        {
            return std::sqrt(x + y);
        },
        1e-30,
        &error);

    const double exact =
        (4.0 / 15.0) * (std::pow(2.0, 2.5) - 2.0);

    CHECK_ABSOLUTE_ERROR(value, exact, 1e-4);
    CHECK_GE(error, 0.0);
}

void test_oscillatory_convergence_guard()
{
    static constexpr double frequency = 9.41;
    const double exact = std::sin(frequency) / frequency;

    double error = 0;
    std::size_t evaluations = 0;
    const double value = van_den_bos_unit_square(
        [](double x, double)
        {
            return std::cos(frequency * x);
        },
        1e-6,
        &error,
        nullptr,
        &evaluations);

    CHECK_ABSOLUTE_ERROR(value, exact, 1e-12);
    CHECK_LE(error, 1e-6);
    CHECK_GE(evaluations, std::size_t(109));
}

void test_float_roundoff_floor()
{
    namespace detail = boost::math::quadrature::detail;

    double error = 0;
    std::size_t evaluations = 0;
    const auto value = van_den_bos_unit_square(
        [](double, double)
        {
            return float(1);
        },
        1e-10,
        &error,
        nullptr,
        &evaluations);

    static_assert(std::is_same<decltype(value), const float>::value,
        "The quadrature must preserve the integrand's floating-point type.");
    const float epsilon = (std::numeric_limits<float>::epsilon)();
    CHECK_ABSOLUTE_ERROR(value, 1.0f, 64 * epsilon);
    CHECK_GE(error, double(8 * epsilon));
    const auto final_level = detail::van_den_bos_unit_square_level(
        detail::van_den_bos_unit_square_level_count - 1);
    std::size_t expected_evaluations = 0;
    for (std::size_t i = 0; i < final_level.size; ++i)
        expected_evaluations += d4_orbit(final_level.point(i)).size();
    CHECK_EQUAL(evaluations, expected_evaluations);
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
                [=](double x, double y)
                {
                    return std::exp(
                        scale * (x + y) / 8.0);
                },
                1e-9,
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

void test_integrable_boundary_singularity()
{
    bool evaluated_boundary = false;
    double error_estimate = 0;
    double l1 = 0;
    std::size_t evaluations = 0;
    const double value = van_den_bos_unit_square(
        [&](double x, double y)
        {
            if (!(x > 0 && x < 1 && y > 0 && y < 1))
                evaluated_boundary = true;
            return 1 / std::sqrt(x * y);
        },
        1e-12,
        &error_estimate,
        &l1,
        &evaluations);

    CHECK_EQUAL(evaluated_boundary, false);
    CHECK_EQUAL(std::isfinite(value), true);
    CHECK_EQUAL(std::isfinite(error_estimate), true);
    CHECK_EQUAL(std::isfinite(l1), true);
    CHECK_ABSOLUTE_ERROR(value, 4.0, 0.1);
    CHECK_EQUAL(evaluations, std::size_t(2197));
}

} // namespace

int main()
{
    test_table_properties();
    test_polynomial();
    test_exp_x_plus_y();
    test_dilogarithm_integral();
    test_square_lattice_green_function();
    test_complex_integrand();
    test_final_precomputed_level();
    test_oscillatory_convergence_guard();
    test_float_roundoff_floor();
    test_integrable_boundary_singularity();
    test_concurrent_integrations();

    return boost::math::test::report_errors();
}
