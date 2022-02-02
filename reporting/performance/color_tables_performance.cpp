//  (C) Copyright Matt Borland 2022.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <array>
#include <random>
#include <boost/math/tools/color_maps.hpp>
#include <benchmark/benchmark.h>

template <typename Real>
void color_table_benchmark(benchmark::State& state)
{
    std::random_device rd;
    std::mt19937_64 gen(rd());
    std::uniform_real_distribution<Real> dist(0, 0.125);
    constexpr boost::math::tools::smooth_cool_warm_color_map smooth_cool_warm;
    Real x = dist(gen);
    for (auto _ : state)
    {
        benchmark::DoNotOptimize(smooth_cool_warm(x));
	x += std::numeric_limits<Real>::epsilon();
    }
}

BENCHMARK_TEMPLATE(color_table_benchmark, float);
BENCHMARK_TEMPLATE(color_table_benchmark, double);

BENCHMARK_MAIN();
