//  (C) Copyright Nick Thompson 2020.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <random>
#include <boost/math/statistics/univariate_statistics.hpp>
#include <boost/multiprecision/cpp_bin_float.hpp>
#include <benchmark/benchmark.h>

using boost::multiprecision::cpp_bin_float_50;

template <class Z>
void test_mode(benchmark::State& state)
{
    using boost::math::statistics::mode;
    
    std::random_device rd;
    std::mt19937_64 mt(rd());
    std::uniform_real_distribution<double> dist {1, 10};

    auto gen = [&dist, &mt](){return dist(mt);};

    std::vector<Z> v(100);
    std::generate(v.begin(), v.end(), gen);

    for (auto _ : state)
    {
        benchmark::DoNotOptimize(mode(v));
    }
}

BENCHMARK_TEMPLATE(test_mode, float);
BENCHMARK_TEMPLATE(test_mode, double);
BENCHMARK_TEMPLATE(test_mode, long double);
BENCHMARK_TEMPLATE(test_mode, int);
BENCHMARK_TEMPLATE(test_mode, cpp_bin_float_50);

BENCHMARK_MAIN();
