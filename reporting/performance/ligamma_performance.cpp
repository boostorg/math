//  (C) Copyright Nick Thompson 2026.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

//
// Compares ligamma(a, x) with log(tgamma(a, x)) where both are finite, even in float: a log-uniform in [0.5, 30],
// x/a log-uniform in [0.1, 10]. Then, for a log-uniform in [1e3, 1e5], where tgamma(a, x) overflows,
// compares ligamma with lgamma(a) + lgamma_q(a, x).
//
// Promotion is disabled so that each type is timed in its own arithmetic: by default,
// float is evaluated in double, and double in long double where that is wider.
//

#include <cmath>
#include <random>
#include <utility>
#include <vector>
#include <benchmark/benchmark.h>
#include <boost/math/special_functions/gamma.hpp>

using no_promotion = boost::math::policies::policy<boost::math::policies::promote_float<false>, boost::math::policies::promote_double<false>>;

template <typename Real>
std::vector<std::pair<Real, Real>> arguments(double log_a_min, double log_a_max)
{
    std::mt19937_64 mt {42};
    std::uniform_real_distribution<double> log_a(log_a_min, log_a_max);
    std::uniform_real_distribution<double> log_ratio(-1, 1);
    std::vector<std::pair<Real, Real>> v(1024);
    for (auto& p : v)
    {
        double a = std::pow(10.0, log_a(mt));
        p = {static_cast<Real>(a), static_cast<Real>(a * std::pow(10.0, log_ratio(mt)))};
    }
    return v;
}

template <typename Real>
void ligamma_performance(benchmark::State& state)
{
    auto v = arguments<Real>(std::log10(0.5), std::log10(30.0));
    std::size_t i = 0;
    for (auto _ : state)
    {
        auto const & p = v[i++ & 1023];
        benchmark::DoNotOptimize(boost::math::ligamma(p.first, p.second, no_promotion()));
    }
}

template <typename Real>
void log_tgamma_performance(benchmark::State& state)
{
    using std::log;
    auto v = arguments<Real>(std::log10(0.5), std::log10(30.0));
    std::size_t i = 0;
    for (auto _ : state)
    {
        auto const & p = v[i++ & 1023];
        benchmark::DoNotOptimize(log(boost::math::tgamma(p.first, p.second, no_promotion())));
    }
}

template <typename Real>
void ligamma_large_a_performance(benchmark::State& state)
{
    auto v = arguments<Real>(3, 5);
    std::size_t i = 0;
    for (auto _ : state)
    {
        auto const & p = v[i++ & 1023];
        benchmark::DoNotOptimize(boost::math::ligamma(p.first, p.second, no_promotion()));
    }
}

template <typename Real>
void lgamma_plus_lgamma_q_large_a_performance(benchmark::State& state)
{
    auto v = arguments<Real>(3, 5);
    std::size_t i = 0;
    for (auto _ : state)
    {
        auto const & p = v[i++ & 1023];
        benchmark::DoNotOptimize(boost::math::lgamma(p.first, no_promotion()) + boost::math::lgamma_q(p.first, p.second, no_promotion()));
    }
}

BENCHMARK_TEMPLATE(ligamma_performance, float);
BENCHMARK_TEMPLATE(ligamma_performance, double);
BENCHMARK_TEMPLATE(ligamma_performance, long double);
BENCHMARK_TEMPLATE(log_tgamma_performance, float);
BENCHMARK_TEMPLATE(log_tgamma_performance, double);
BENCHMARK_TEMPLATE(log_tgamma_performance, long double);
BENCHMARK_TEMPLATE(ligamma_large_a_performance, double);
BENCHMARK_TEMPLATE(lgamma_plus_lgamma_q_large_a_performance, double);

BENCHMARK_MAIN();
