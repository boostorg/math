// Copyright Nick Thompson 2026.
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0. (See accompanying file
// LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

/*
 * Google Benchmark suite for NSGA-II: 20 generations on one thread of the (cheap) ZDT1 and DTLZ2 problems in 12
 * variables, of 2 and 3 objectives, so that the time is the optimizer's overhead. The per_call counter is the time per
 * function call. The argument is the population size.
 *
 * Example:
 *   c++ -O2 -DNDEBUG -std=c++17 -Iinclude \
 *     reporting/performance/nsga2_performance.cpp \
 *     -lbenchmark -lpthread -o nsga2_bench
 */

#include <boost/math/optimization/nsga2.hpp>
#include <benchmark/benchmark.h>

#include "../../test/test_functions_for_optimization.hpp"

#include <random>
#include <vector>

namespace {

template <class Objectives>
void run(benchmark::State& state, Objectives objectives)
{
    auto const population_size = static_cast<std::size_t>(state.range(0));
    boost::math::optimization::nsga2_parameters<std::vector<double>> params;
    params.lower_bounds.assign(12, 0.0);
    params.upper_bounds.assign(12, 1.0);
    params.population_size = population_size;
    params.max_generations = 20;
    params.threads = 1;
    std::mt19937_64 gen(42);
    for (auto _ : state)
    {
        benchmark::DoNotOptimize(boost::math::optimization::nsga2(objectives, params, gen).data());
    }
    state.counters["per_call"] = benchmark::Counter(static_cast<double>(population_size * (params.max_generations + 1)),
        benchmark::Counter::kIsIterationInvariantRate | benchmark::Counter::kInvert);
}
void BM_nsga2_zdt1(benchmark::State& state) { run(state, [](std::vector<double> const& x) { return zdt(x, 1); }); }
void BM_nsga2_dtlz2(benchmark::State& state) { run(state, dtlz2<std::vector<double>>); }
BENCHMARK(BM_nsga2_zdt1)->Arg(100)->Arg(400)->Arg(1600)->Unit(benchmark::kMillisecond)->MeasureProcessCPUTime()->UseRealTime();
BENCHMARK(BM_nsga2_dtlz2)->Arg(100)->Arg(400)->Arg(1600)->Unit(benchmark::kMillisecond)->MeasureProcessCPUTime()->UseRealTime();

} // namespace

BENCHMARK_MAIN();
