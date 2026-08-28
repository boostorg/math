/*
 * Google Benchmark suite for van_den_bos_unit_square integrands.
 *
 * Example:
 *   c++ -O3 -DNDEBUG -std=c++17 -Iinclude \
 *     reporting/performance/van_den_bos_unit_square_performance.cpp \
 *     -lbenchmark -lpthread -o vdb_bench
 */

#include <boost/math/quadrature/van_den_bos_unit_square.hpp>
#include <benchmark/benchmark.h>

#include <array>
#include <cmath>
#include <complex>
#include <cstddef>

namespace {

using point2 = std::array<double, 2>;
using boost::math::quadrature::van_den_bos_unit_square;

template <class F>
void benchmark_integrand(benchmark::State& state, F const& f, double tol)
{
    for (auto _ : state)
    {
        const auto value = van_den_bos_unit_square(f, tol, 5);
        benchmark::DoNotOptimize(value);
    }
}

static void BM_quadratic(benchmark::State& state)
{
    benchmark_integrand(
        state,
        [](point2 const& p)
        {
            return p[0] * p[0] + p[1] * p[1] + p[0] * p[1];
        },
        1e-12);
}
BENCHMARK(BM_quadratic);

static void BM_exp_x_plus_y(benchmark::State& state)
{
    benchmark_integrand(
        state,
        [](point2 const& p)
        {
            return std::exp(p[0] + p[1]);
        },
        1e-10);
}
BENCHMARK(BM_exp_x_plus_y);

static void BM_exp_xy(benchmark::State& state)
{
    benchmark_integrand(
        state,
        [](point2 const& p)
        {
            return std::exp(p[0] * p[1]);
        },
        1e-10);
}
BENCHMARK(BM_exp_xy);

static void BM_one_over_1_plus_x_plus_y(benchmark::State& state)
{
    benchmark_integrand(
        state,
        [](point2 const& p)
        {
            return 1.0 / (1.0 + p[0] + p[1]);
        },
        1e-10);
}
BENCHMARK(BM_one_over_1_plus_x_plus_y);

static void BM_cos_5xy(benchmark::State& state)
{
    benchmark_integrand(
        state,
        [](point2 const& p)
        {
            return std::cos(5.0 * p[0] * p[1]);
        },
        1e-8);
}
BENCHMARK(BM_cos_5xy);

static void BM_complex_exp_5ixy(benchmark::State& state)
{
    benchmark_integrand(
        state,
        [](point2 const& p)
        {
            return std::exp(std::complex<double>(0, 5 * p[0] * p[1]));
        },
        1e-8);
}
BENCHMARK(BM_complex_exp_5ixy);

} // namespace

BENCHMARK_MAIN();
