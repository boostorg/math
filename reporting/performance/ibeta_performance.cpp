//  (C) Copyright Nick Thompson 2026.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

//
// Times ibeta, ibetac and ibeta_invb over the regions whose evaluation was changed to meet the
// accuracy of the issue "Inverse incomplete gamma and beta functions are 10-35x less accurate than
// conditioning allows":
//
//   small:  a and b log-uniform in [0.01, 1], x uniform in (0, 1): the series, summed compensated
//   cf:     a log-uniform in [1, 15], b log-uniform in [1, 40], x within 10% of the mean a / (a + b): the
//           continued fraction, evaluated from the tail, in place of a sidestep to the series
//   mean:   a and b log-uniform in [15, 100], x within 5% of the mean: the power terms in the prefix
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
std::vector<std::tuple<Real, Real, Real>> arguments(double a_min, double a_max, double b_min, double b_max, double mean_jitter)
{
    std::mt19937_64 mt {42};
    std::uniform_real_distribution<double> u(0, 1);
    std::vector<std::tuple<Real, Real, Real>> v(1024);
    for (auto& t : v)
    {
        double a = a_min * std::pow(a_max / a_min, u(mt));
        double b = b_min * std::pow(b_max / b_min, u(mt));
        double x = mean_jitter > 0 ? a / (a + b) * (1 + mean_jitter * (2 * u(mt) - 1)) : u(mt);
        t = std::make_tuple(static_cast<Real>(a), static_cast<Real>(b), static_cast<Real>(x));
    }
    return v;
}

template <typename Real>
void ibeta_small(benchmark::State& state)
{
    auto v = arguments<Real>(0.01, 1, 0.01, 1, 0);
    std::size_t i = 0;
    for (auto _ : state)
    {
        auto const & t = v[i++ & 1023];
        benchmark::DoNotOptimize(boost::math::ibeta(std::get<0>(t), std::get<1>(t), std::get<2>(t), no_promotion()));
    }
}

template <typename Real>
void ibetac_small(benchmark::State& state)
{
    auto v = arguments<Real>(0.01, 1, 0.01, 1, 0);
    std::size_t i = 0;
    for (auto _ : state)
    {
        auto const & t = v[i++ & 1023];
        benchmark::DoNotOptimize(boost::math::ibetac(std::get<0>(t), std::get<1>(t), std::get<2>(t), no_promotion()));
    }
}

template <typename Real>
void ibeta_cf(benchmark::State& state)
{
    auto v = arguments<Real>(1, 15, 1, 40, 0.1);
    std::size_t i = 0;
    for (auto _ : state)
    {
        auto const & t = v[i++ & 1023];
        benchmark::DoNotOptimize(boost::math::ibeta(std::get<0>(t), std::get<1>(t), std::get<2>(t), no_promotion()));
    }
}

template <typename Real>
void ibeta_mean(benchmark::State& state)
{
    auto v = arguments<Real>(15, 100, 15, 100, 0.05);
    std::size_t i = 0;
    for (auto _ : state)
    {
        auto const & t = v[i++ & 1023];
        benchmark::DoNotOptimize(boost::math::ibeta(std::get<0>(t), std::get<1>(t), std::get<2>(t), no_promotion()));
    }
}

template <typename Real>
void ibeta_invb_small(benchmark::State& state)
{
    // b is found from a, x and p = ibeta(a, b, x) for b in (0.01, 1):
    auto v = arguments<Real>(0.01, 1, 0.01, 1, 0);
    for (auto& t : v)
        std::get<1>(t) = boost::math::ibeta(std::get<0>(t), std::get<1>(t), std::get<2>(t), no_promotion());
    std::size_t i = 0;
    for (auto _ : state)
    {
        auto const & t = v[i++ & 1023];
        benchmark::DoNotOptimize(boost::math::ibeta_invb(std::get<0>(t), std::get<2>(t), std::get<1>(t), no_promotion()));
    }
}

BENCHMARK_TEMPLATE(ibeta_small, float);
BENCHMARK_TEMPLATE(ibeta_small, double);
BENCHMARK_TEMPLATE(ibetac_small, double);
BENCHMARK_TEMPLATE(ibeta_cf, float);
BENCHMARK_TEMPLATE(ibeta_cf, double);
BENCHMARK_TEMPLATE(ibeta_mean, double);
BENCHMARK_TEMPLATE(ibeta_invb_small, double);

BENCHMARK_MAIN();
