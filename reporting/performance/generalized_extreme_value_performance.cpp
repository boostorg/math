//  (C) Copyright Nick Thompson 2026.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

//
// The generalized extreme value distribution's pdf, cdf and quantile at shape xi = 0.2, against the
// Gumbel distribution (extreme_value, which is xi = 0) as the baseline, and its skewness, whose
// moment series is the expensive part of the moments.
//
// Promotion is disabled so that each type is timed in its own arithmetic: by default,
// float is evaluated in double, and double in long double where that is wider.
//

#include <random>
#include <vector>
#include <benchmark/benchmark.h>
#include <boost/math/distributions/extreme_value.hpp>
#include <boost/math/distributions/generalized_extreme_value.hpp>

using no_promotion = boost::math::policies::policy<boost::math::policies::promote_float<false>, boost::math::policies::promote_double<false>>;

template <typename Real>
std::vector<Real> abscissas()
{
    std::mt19937_64 mt {42};
    std::uniform_real_distribution<double> dis(-2, 20);
    std::vector<Real> v(1024);
    for (auto& x : v)
    {
        x = static_cast<Real>(dis(mt));
    }
    return v;
}

template <typename Real>
std::vector<Real> probabilities()
{
    std::mt19937_64 mt {42};
    std::uniform_real_distribution<double> dis(0.001, 0.999);
    std::vector<Real> v(1024);
    for (auto& p : v)
    {
        p = static_cast<Real>(dis(mt));
    }
    return v;
}

template <typename Real>
auto gev()
{
    return boost::math::generalized_extreme_value_distribution<Real, no_promotion>(Real(0.5), Real(2), Real(0.2));
}

template <typename Real>
auto gumbel()
{
    return boost::math::extreme_value_distribution<Real, no_promotion>(Real(0.5), Real(2));
}

#define BENCH(name, make, fn, points)                     \
template <typename Real>                                  \
void name(benchmark::State& state)                        \
{                                                         \
    auto d = make<Real>();                                \
    auto v = points<Real>();                              \
    std::size_t i = 0;                                    \
    for (auto _ : state)                                  \
    {                                                     \
        benchmark::DoNotOptimize(fn(d, v[i++ & 1023]));   \
    }                                                     \
}

BENCH(gev_pdf, gev, boost::math::pdf, abscissas)
BENCH(gumbel_pdf, gumbel, boost::math::pdf, abscissas)
BENCH(gev_cdf, gev, boost::math::cdf, abscissas)
BENCH(gumbel_cdf, gumbel, boost::math::cdf, abscissas)
BENCH(gev_quantile, gev, boost::math::quantile, probabilities)
BENCH(gumbel_quantile, gumbel, boost::math::quantile, probabilities)

template <typename Real>
void gev_skewness(benchmark::State& state)
{
    std::mt19937_64 mt {42};
    std::uniform_real_distribution<double> dis(-0.5, 0.3);
    std::vector<boost::math::generalized_extreme_value_distribution<Real, no_promotion>> v;
    for (int i = 0; i < 1024; ++i)
    {
        v.emplace_back(Real(0), Real(1), static_cast<Real>(dis(mt)));
    }
    std::size_t i = 0;
    for (auto _ : state)
    {
        benchmark::DoNotOptimize(boost::math::skewness(v[i++ & 1023]));
    }
}

BENCHMARK_TEMPLATE(gev_pdf, float);
BENCHMARK_TEMPLATE(gev_pdf, double);
BENCHMARK_TEMPLATE(gumbel_pdf, double);
BENCHMARK_TEMPLATE(gev_cdf, float);
BENCHMARK_TEMPLATE(gev_cdf, double);
BENCHMARK_TEMPLATE(gumbel_cdf, double);
BENCHMARK_TEMPLATE(gev_quantile, float);
BENCHMARK_TEMPLATE(gev_quantile, double);
BENCHMARK_TEMPLATE(gumbel_quantile, double);
BENCHMARK_TEMPLATE(gev_skewness, double);

BENCHMARK_MAIN();
