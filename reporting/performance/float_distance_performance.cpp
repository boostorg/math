//  (C) Copyright Matt Borland 2022.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <boost/math/special_functions/next.hpp>
#include <benchmark/benchmark.h>

template <typename T>
void float_distance(benchmark::State& state)
{
    const auto difference = static_cast<int>(state.range(0));
    T left = 2;
    T right = boost::math::float_advance(left, difference);

    for (auto _ : state)
    {
        benchmark::DoNotOptimize(left);
        benchmark::DoNotOptimize(right);
        benchmark::DoNotOptimize(boost::math::float_distance(left, right));
    }
    state.SetComplexityN(state.range(0));
}

BENCHMARK_TEMPLATE(float_distance, float)->RangeMultiplier(2)->Range(1 << 1, 1 << 14)->Complexity()->UseRealTime();
BENCHMARK_TEMPLATE(float_distance, double)->RangeMultiplier(2)->Range(1 << 1, 1 << 14)->Complexity()->UseRealTime();

BENCHMARK_MAIN();

/*
Apple M-series, Homebrew clang 23, -O3 (Sept 2026), time per call:

                     generic code    bit patterns
float_distance<float>     14.8 ns         2.3 ns
float_distance<double>    14.5 ns         2.2 ns
*/
