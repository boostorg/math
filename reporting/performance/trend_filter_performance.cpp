// Copyright Nick Thompson 2026.
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0. (See accompanying file
// LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

/*
 * Google Benchmark suite for the trend filter: one fit at fixed k and lambda (k = 0: Condat's
 * algorithm; k = 1, 2: interior point), and the cross-validated filter, n = 2^9 .. 2^14.
 * The data is Doppler (Donoho & Johnstone) plus Gaussian noise at SNR 7.
 *
 * Example:
 *   c++ -O3 -DNDEBUG -std=c++14 -Iinclude \
 *     reporting/performance/trend_filter_performance.cpp \
 *     -lbenchmark -lpthread -o trend_filter_bench
 */

#include <boost/math/filters/trend_filter.hpp>
#include <benchmark/benchmark.h>

#include <cmath>
#include <random>
#include <vector>

namespace {

using boost::math::filters::trend_filter;

template <class Real>
std::vector<Real> noisy_doppler(std::size_t n)
{
    std::mt19937_64 gen(12345);
    std::normal_distribution<Real> noise(0, 1);
    std::vector<Real> y(n);
    for (std::size_t i = 0; i < n; ++i)
    {
        const Real t = Real(i) / Real(n);
        y[i] = 7 * 3.3 * std::sqrt(t * (1 - t)) * std::sin(2 * 3.14159265358979 * 1.05 / (t + 0.05)) + noise(gen);
    }
    return y;
}

template <class Real>
void BM_trend_filter_fixed(benchmark::State& state)
{
    const auto y = noisy_doppler<Real>(static_cast<std::size_t>(state.range(0)));
    const int k = static_cast<int>(state.range(1));
    // a lambda near the one cross-validation picks for this signal
    const Real lambda = k == 0 ? Real(10) : k == 1 ? Real(30) : Real(300);
    for (auto _ : state)
    {
        benchmark::DoNotOptimize(trend_filter(y, k, lambda).data());
    }
    state.SetComplexityN(state.range(0));
}
BENCHMARK_TEMPLATE(BM_trend_filter_fixed, double)->ArgsProduct({benchmark::CreateRange(512, 16384, 2), {0, 1, 2}})->Unit(benchmark::kMillisecond);

template <class Real>
void BM_trend_filter_cv(benchmark::State& state)
{
    const auto y = noisy_doppler<Real>(static_cast<std::size_t>(state.range(0)));
    for (auto _ : state)
    {
        benchmark::DoNotOptimize(trend_filter(y).data());
    }
    state.SetComplexityN(state.range(0));
}
BENCHMARK_TEMPLATE(BM_trend_filter_cv, float)->RangeMultiplier(2)->Range(512, 16384)->Unit(benchmark::kMillisecond)->Complexity();
BENCHMARK_TEMPLATE(BM_trend_filter_cv, double)->RangeMultiplier(2)->Range(512, 16384)->Unit(benchmark::kMillisecond)->Complexity();

} // namespace

BENCHMARK_MAIN();
