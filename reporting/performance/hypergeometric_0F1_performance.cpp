//  (C) Copyright Nick Thompson 2026.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <benchmark/benchmark.h>
#include <boost/math/special_functions/hypergeometric_0F1.hpp>

// Arguments: 2b, z.
template <class Real>
void Hypergeometric0F1(benchmark::State& state)
{
   Real b = static_cast<Real>(state.range(0)) / 2;
   Real z = static_cast<Real>(state.range(1));
   for (auto _ : state)
   {
      benchmark::DoNotOptimize(b);
      benchmark::DoNotOptimize(z);
      benchmark::DoNotOptimize(boost::math::hypergeometric_0F1(b, z));
   }
}

BENCHMARK_TEMPLATE(Hypergeometric0F1, double)->ArgsProduct({ { -9, 1, 3, 15, 41, 400 }, { -500, -170, -100, -30, -8, -3, -1, 5, 30 } });

BENCHMARK_MAIN();
