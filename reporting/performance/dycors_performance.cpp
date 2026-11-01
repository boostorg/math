// Copyright Nick Thompson 2026.
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0. (See accompanying file
// LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

/*
 * Google Benchmark suite for DYCORS: a whole run of 500 calls on the (cheap) Ackley function, so that the time is the
 * optimizer's overhead. The per_call counter is the time per function call.
 *
 * Example:
 *   c++ -O2 -DNDEBUG -std=c++17 -Iinclude \
 *     reporting/performance/dycors_performance.cpp \
 *     -lbenchmark -lpthread -o dycors_bench
 */

#include <boost/math/optimization/dycors.hpp>
#include <benchmark/benchmark.h>

#include <cmath>
#include <random>
#include <vector>

namespace {

double ackley(std::vector<double> const& x)
{
    double const two_pi = 6.283185307179586;
    double a = 0;
    double b = 0;
    for (double t : x)
    {
        a += t * t;
        b += std::cos(two_pi * t);
    }
    double const n = static_cast<double>(x.size());
    return -20 * std::exp(-0.2 * std::sqrt(a / n)) - std::exp(b / n) + 20 + std::exp(1.0);
}

void BM_dycors_ackley(benchmark::State& state)
{
    auto const d = static_cast<std::size_t>(state.range(0));
    std::size_t const calls = 500;
    boost::math::optimization::dycors_parameters<std::vector<double>> params;
    params.lower_bounds.assign(d, -15);
    params.upper_bounds.assign(d, 20);
    params.max_function_calls = calls;
    std::mt19937_64 gen(42);
    for (auto _ : state)
    {
        benchmark::DoNotOptimize(boost::math::optimization::dycors(ackley, params, gen).data());
    }
    state.counters["per_call"] = benchmark::Counter(static_cast<double>(calls),
        benchmark::Counter::kIsIterationInvariantRate | benchmark::Counter::kInvert);
}
BENCHMARK(BM_dycors_ackley)->Arg(10)->Arg(30)->Arg(100)->Unit(benchmark::kMillisecond)->MeasureProcessCPUTime()->UseRealTime();

} // namespace

BENCHMARK_MAIN();
