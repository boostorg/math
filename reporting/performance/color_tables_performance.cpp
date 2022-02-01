//  (C) Copyright Matt Borland 2022.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <array>
#include <random>
#include <boost/math/tools/color_maps.hpp>
#include <benchmark/benchmark.h>

template <typename Table, typename Dist, typename Gen>
int helper(const Table& table, Dist d, Gen gen, std::int64_t size)
{
    for (std::int64_t i = 0; i < size; ++i)
    {
        table(d(gen));
    }

    return 0;
}

template <typename T>
void color_table_benchmark(benchmark::State& state)
{
    std::random_device rd;
    std::mt19937_64 gen(rd());
    std::uniform_real_distribution<T> dist(0, 1);
    constexpr boost::math::tools::smooth_cool_warm_color_map smooth_cool_warm;
    const std::int64_t size = state.range(0);

    for (auto _ : state)
    {
        benchmark::DoNotOptimize(helper(smooth_cool_warm, dist, gen, size));
    }
    state.SetComplexityN(state.range(0));
}

BENCHMARK_TEMPLATE(color_table_benchmark, float)->RangeMultiplier(2)->Range(1 << 6, 1 << 20)->Complexity(benchmark::oN)->UseRealTime();
BENCHMARK_TEMPLATE(color_table_benchmark, double)->RangeMultiplier(2)->Range(1 << 6, 1 << 20)->Complexity(benchmark::oN)->UseRealTime();
BENCHMARK_TEMPLATE(color_table_benchmark, long double)->RangeMultiplier(2)->Range(1 << 6, 1 << 20)->Complexity(benchmark::oN)->UseRealTime();

BENCHMARK_MAIN();
