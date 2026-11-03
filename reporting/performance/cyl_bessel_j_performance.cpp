//  (C) Copyright Nick Thompson 2026.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <benchmark/benchmark.h>
#include <boost/math/special_functions/bessel.hpp>

// Arguments: n, 1000x.
template <class Real>
void CylBesselJ(benchmark::State& state)
{
   int n = static_cast<int>(state.range(0));
   Real x = static_cast<Real>(state.range(1)) / 1000;
   for (auto _ : state)
   {
      benchmark::DoNotOptimize(x);
      benchmark::DoNotOptimize(boost::math::cyl_bessel_j(n, x));
   }
}

BENCHMARK_TEMPLATE(CylBesselJ, double)->ArgsProduct({ { 0, 1, 7 }, { 2500, 4223, 7695, 12000 } });

BENCHMARK_MAIN();
