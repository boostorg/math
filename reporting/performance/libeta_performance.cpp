//  (C) Copyright Nick Thompson 2026.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

//
// Compares libeta(a, b, x) with log(ibeta(a, b, x)) where both are finite: a and b log-uniform
// in [0.5, 1000], x uniform in (0, 1). Then times libeta alone in the left tail of large a = b,
// where ibeta underflows and log(ibeta) is -infinity.
//
// Promotion is disabled so that each type is timed in its own arithmetic: by default,
// float is evaluated in double, and double in long double where that is wider.
//

#include <cmath>
#include <random>
#include <tuple>
#include <vector>
#include <benchmark/benchmark.h>
#include <boost/math/special_functions/beta.hpp>

using no_promotion = boost::math::policies::policy<boost::math::policies::promote_float<false>, boost::math::policies::promote_double<false>>;

template <typename Real>
std::vector<std::tuple<Real, Real, Real>> ordinary_arguments()
{
    std::mt19937_64 mt {42};
    std::uniform_real_distribution<double> log_ab(std::log10(0.5), 3);
    std::uniform_real_distribution<double> unit(0, 1);
    std::vector<std::tuple<Real, Real, Real>> v(1024);
    for (auto& t : v)
    {
        t = {static_cast<Real>(std::pow(10.0, log_ab(mt))), static_cast<Real>(std::pow(10.0, log_ab(mt))), static_cast<Real>(unit(mt))};
    }
    return v;
}

template <typename Real>
std::vector<std::tuple<Real, Real, Real>> underflow_arguments()
{
    std::mt19937_64 mt {42};
    std::uniform_real_distribution<double> log_a(4, 5);
    std::uniform_real_distribution<double> fraction(0.05, 0.4);
    std::vector<std::tuple<Real, Real, Real>> v(1024);
    for (auto& t : v)
    {
        Real a = static_cast<Real>(std::pow(10.0, log_a(mt)));
        t = {a, a, static_cast<Real>(fraction(mt))};
    }
    return v;
}

template <typename Real>
void libeta_performance(benchmark::State& state)
{
    auto v = ordinary_arguments<Real>();
    std::size_t i = 0;
    for (auto _ : state)
    {
        auto const & t = v[i++ & 1023];
        benchmark::DoNotOptimize(boost::math::libeta(std::get<0>(t), std::get<1>(t), std::get<2>(t), no_promotion()));
    }
}

template <typename Real>
void log_ibeta_performance(benchmark::State& state)
{
    using std::log;
    auto v = ordinary_arguments<Real>();
    std::size_t i = 0;
    for (auto _ : state)
    {
        auto const & t = v[i++ & 1023];
        benchmark::DoNotOptimize(log(boost::math::ibeta(std::get<0>(t), std::get<1>(t), std::get<2>(t), no_promotion())));
    }
}

template <typename Real>
void libeta_underflow_performance(benchmark::State& state)
{
    auto v = underflow_arguments<Real>();
    std::size_t i = 0;
    for (auto _ : state)
    {
        auto const & t = v[i++ & 1023];
        benchmark::DoNotOptimize(boost::math::libeta(std::get<0>(t), std::get<1>(t), std::get<2>(t), no_promotion()));
    }
}

BENCHMARK_TEMPLATE(libeta_performance, float);
BENCHMARK_TEMPLATE(libeta_performance, double);
BENCHMARK_TEMPLATE(libeta_performance, long double);
BENCHMARK_TEMPLATE(log_ibeta_performance, float);
BENCHMARK_TEMPLATE(log_ibeta_performance, double);
BENCHMARK_TEMPLATE(log_ibeta_performance, long double);
BENCHMARK_TEMPLATE(libeta_underflow_performance, float);
BENCHMARK_TEMPLATE(libeta_underflow_performance, double);
BENCHMARK_TEMPLATE(libeta_underflow_performance, long double);

BENCHMARK_MAIN();
