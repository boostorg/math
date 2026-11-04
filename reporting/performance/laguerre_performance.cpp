// Copyright Nick Thompson 2026.
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0. (See accompanying file
// LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

/*
 * Google Benchmark suite for the Laguerre polynomials L_n(x), the associated L_n^m(x) and the generalized
 * L_n^alpha(x), for degrees n = 8 to 512 (the benchmark's argument) at x = 3.5.
 *
 * Example:
 *   c++ -O2 -DNDEBUG -std=c++14 -Iinclude \
 *     reporting/performance/laguerre_performance.cpp \
 *     -lbenchmark -lpthread -o laguerre_bench
 */

#include <boost/math/special_functions/laguerre.hpp>
#include <benchmark/benchmark.h>

namespace {

void laguerre_plain(benchmark::State& state)
{
    const unsigned n = static_cast<unsigned>(state.range(0));
    double x = 3.5;
    for (auto _ : state)
    {
        benchmark::DoNotOptimize(x);
        benchmark::DoNotOptimize(boost::math::laguerre(n, x));
    }
}

void laguerre_associated(benchmark::State& state)
{
    const unsigned n = static_cast<unsigned>(state.range(0));
    double x = 3.5;
    for (auto _ : state)
    {
        benchmark::DoNotOptimize(x);
        benchmark::DoNotOptimize(boost::math::laguerre(n, 3u, x));
    }
}

void laguerre_generalized(benchmark::State& state)
{
    const unsigned n = static_cast<unsigned>(state.range(0));
    double x = 3.5;
    for (auto _ : state)
    {
        benchmark::DoNotOptimize(x);
        benchmark::DoNotOptimize(boost::math::laguerre(n, 0.5, x));
    }
}

} // namespace

BENCHMARK(laguerre_plain)->RangeMultiplier(4)->Range(8, 512);
BENCHMARK(laguerre_associated)->RangeMultiplier(4)->Range(8, 512);
BENCHMARK(laguerre_generalized)->RangeMultiplier(4)->Range(8, 512);

BENCHMARK_MAIN();
