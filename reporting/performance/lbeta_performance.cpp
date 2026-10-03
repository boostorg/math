//  (C) Copyright Nick Thompson 2026.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

//
// Compares lbeta(a, b) with the two obvious alternatives, log(beta(a, b)) and
// lgamma(a) + lgamma(b) - lgamma(a + b), with a and b log-uniform in [1e-3, 1e6]:
// the range where neither alternative fails outright for most pairs.
//
// Promotion is disabled so that each type is timed in its own arithmetic: by default,
// float is evaluated in double, and double in long double where that is wider.
//

#include <cmath>
#include <random>
#include <utility>
#include <vector>
#include <benchmark/benchmark.h>
#include <boost/math/special_functions/beta.hpp>
#include <boost/math/special_functions/gamma.hpp>

using no_promotion = boost::math::policies::policy<boost::math::policies::promote_float<false>, boost::math::policies::promote_double<false>>;

template <typename Real>
std::vector<std::pair<Real, Real>> arguments()
{
    std::mt19937_64 mt {42};
    std::uniform_real_distribution<double> dist(-3, 6);
    std::vector<std::pair<Real, Real>> v(1024);
    for (auto& p : v)
    {
        p = {static_cast<Real>(std::pow(10.0, dist(mt))), static_cast<Real>(std::pow(10.0, dist(mt)))};
    }
    return v;
}

template <typename Real>
void lbeta_performance(benchmark::State& state)
{
    auto v = arguments<Real>();
    std::size_t i = 0;
    for (auto _ : state)
    {
        auto const & p = v[i++ & 1023];
        benchmark::DoNotOptimize(boost::math::lbeta(p.first, p.second, no_promotion()));
    }
}

template <typename Real>
void log_beta_performance(benchmark::State& state)
{
    using std::log;
    auto v = arguments<Real>();
    std::size_t i = 0;
    for (auto _ : state)
    {
        auto const & p = v[i++ & 1023];
        benchmark::DoNotOptimize(log(boost::math::beta(p.first, p.second, no_promotion())));
    }
}

template <typename Real>
void lgamma_sum_performance(benchmark::State& state)
{
    using boost::math::lgamma;
    auto v = arguments<Real>();
    std::size_t i = 0;
    for (auto _ : state)
    {
        auto const & p = v[i++ & 1023];
        benchmark::DoNotOptimize(lgamma(p.first, no_promotion()) + lgamma(p.second, no_promotion()) - lgamma(p.first + p.second, no_promotion()));
    }
}

BENCHMARK_TEMPLATE(lbeta_performance, float);
BENCHMARK_TEMPLATE(lbeta_performance, double);
BENCHMARK_TEMPLATE(lbeta_performance, long double);
BENCHMARK_TEMPLATE(log_beta_performance, float);
BENCHMARK_TEMPLATE(log_beta_performance, double);
BENCHMARK_TEMPLATE(log_beta_performance, long double);
BENCHMARK_TEMPLATE(lgamma_sum_performance, float);
BENCHMARK_TEMPLATE(lgamma_sum_performance, double);
BENCHMARK_TEMPLATE(lgamma_sum_performance, long double);

BENCHMARK_MAIN();
