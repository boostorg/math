//  (C) Copyright Nick Thompson 2020.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <random>
#include <benchmark/benchmark.h>
#include <boost/math/tools/agm.hpp>

using boost::math::tools::agm;
template<class Real>
void AGM(benchmark::State& state)
{
    std::random_device rd;
    std::mt19937_64 mt(rd());
    std::uniform_real_distribution<Real> unif(1,100);

    Real x = unif(mt);
    for (auto _ : state)
    {
        benchmark::DoNotOptimize(agm(x,Real(1)));
    }
}

BENCHMARK_TEMPLATE(AGM, float);
BENCHMARK_TEMPLATE(AGM, double);
BENCHMARK_TEMPLATE(AGM, long double);




BENCHMARK_MAIN();