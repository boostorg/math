//  (C) Copyright Nick Thompson 2026.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <cmath>
#include <random>
#include <utility>
#include <vector>
#include <benchmark/benchmark.h>
#include <boost/math/special_functions/pow1p.hpp>
#ifdef BOOST_HAS_FLOAT128
#include <boost/multiprecision/float128.hpp>
#endif
#include <boost/multiprecision/cpp_bin_float.hpp>

#ifdef BOOST_HAS_FLOAT128
using boost::multiprecision::float128;
#endif
using boost::multiprecision::cpp_bin_float_50;

// Moderate arguments: x in (-1/2, 1/2), y in (-50, 50).
template<class Real>
std::vector<std::pair<Real, Real>> moderate_arguments()
{
    std::mt19937_64 mt(12345);
    std::uniform_real_distribution<double> ux(-0.5, 0.5);
    std::uniform_real_distribution<double> uy(-50, 50);
    std::vector<std::pair<Real, Real>> v(1024);
    for (auto& p : v)
    {
        p = {static_cast<Real>(ux(mt)), static_cast<Real>(uy(mt))};
    }
    return v;
}

// Tiny x and huge y with y*x of order 10: 1+x is inexact, so pow1p takes its double-T branch.
template<class Real>
std::vector<std::pair<Real, Real>> tiny_x_arguments()
{
    std::mt19937_64 mt(12345);
    std::uniform_real_distribution<double> ue(-18, -8);
    std::uniform_real_distribution<double> uw(-10, 10);
    std::vector<std::pair<Real, Real>> v(1024);
    for (auto& p : v)
    {
        double x = std::pow(10.0, ue(mt));
        p = {static_cast<Real>(x), static_cast<Real>(uw(mt) / x)};
    }
    return v;
}

template<class Real>
void Pow1p(benchmark::State& state)
{
    auto v = state.range(0) ? tiny_x_arguments<Real>() : moderate_arguments<Real>();
    std::size_t i = 0;
    for (auto _ : state)
    {
        benchmark::DoNotOptimize(boost::math::pow1p(v[i].first, v[i].second));
        i = (i + 1) & 1023;
    }
}

// The naive, inaccurate evaluation, for comparison.
template<class Real>
void PowOnePlusX(benchmark::State& state)
{
    using std::pow;
    auto v = state.range(0) ? tiny_x_arguments<Real>() : moderate_arguments<Real>();
    std::size_t i = 0;
    for (auto _ : state)
    {
        benchmark::DoNotOptimize(Real(pow(Real(1 + v[i].first), v[i].second)));
        i = (i + 1) & 1023;
    }
}

// Argument 0: moderate x and y.  Argument 1: tiny x, huge y.
BENCHMARK_TEMPLATE(Pow1p, float)->Arg(0)->Arg(1);
BENCHMARK_TEMPLATE(PowOnePlusX, float)->Arg(0)->Arg(1);
BENCHMARK_TEMPLATE(Pow1p, double)->Arg(0)->Arg(1);
BENCHMARK_TEMPLATE(PowOnePlusX, double)->Arg(0)->Arg(1);
BENCHMARK_TEMPLATE(Pow1p, long double)->Arg(0)->Arg(1);
BENCHMARK_TEMPLATE(PowOnePlusX, long double)->Arg(0)->Arg(1);
#ifdef BOOST_HAS_FLOAT128
BENCHMARK_TEMPLATE(Pow1p, float128)->Arg(0)->Arg(1);
BENCHMARK_TEMPLATE(PowOnePlusX, float128)->Arg(0)->Arg(1);
#endif
BENCHMARK_TEMPLATE(Pow1p, cpp_bin_float_50)->Arg(0)->Arg(1);
BENCHMARK_TEMPLATE(PowOnePlusX, cpp_bin_float_50)->Arg(0)->Arg(1);

BENCHMARK_MAIN();
