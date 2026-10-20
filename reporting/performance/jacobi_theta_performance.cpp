//  (C) Copyright Nick Thompson 2020.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <random>
#include <benchmark/benchmark.h>
#include <boost/math/special_functions/jacobi_theta.hpp>
#ifdef BOOST_HAS_FLOAT128
#include <boost/multiprecision/float128.hpp>
#endif
#include <boost/multiprecision/mpfr.hpp>
#include <boost/multiprecision/cpp_bin_float.hpp>

using boost::multiprecision::number;
using boost::multiprecision::mpfr_float_backend;
#ifdef BOOST_HAS_FLOAT128
using boost::multiprecision::float128;
#endif
using boost::multiprecision::cpp_bin_float_50;
using boost::multiprecision::cpp_bin_float_100;
using boost::math::jacobi_theta1;
using boost::math::jacobi_theta1tau;

template<class Real>
void JacobiTheta1(benchmark::State& state)
{
    std::random_device rd;
    std::mt19937_64 mt(rd());
    std::uniform_real_distribution<long double> unif(0,0.01);

    Real x = static_cast<Real>(unif(mt));
    Real q = static_cast<Real>(unif(mt));
    for (auto _ : state)
    {
        benchmark::DoNotOptimize(jacobi_theta1(x, q));
        x += std::numeric_limits<Real>::epsilon();
    }
}

BENCHMARK_TEMPLATE(JacobiTheta1, float);
BENCHMARK_TEMPLATE(JacobiTheta1, double);
BENCHMARK_TEMPLATE(JacobiTheta1, long double);
#ifdef BOOST_HAS_FLOAT128
BENCHMARK_TEMPLATE(JacobiTheta1, float128);
#endif
BENCHMARK_TEMPLATE(JacobiTheta1, number<mpfr_float_backend<100>>);
BENCHMARK_TEMPLATE(JacobiTheta1, number<mpfr_float_backend<200>>);
BENCHMARK_TEMPLATE(JacobiTheta1, number<mpfr_float_backend<300>>);
BENCHMARK_TEMPLATE(JacobiTheta1, number<mpfr_float_backend<400>>);
BENCHMARK_TEMPLATE(JacobiTheta1, number<mpfr_float_backend<1000>>);
BENCHMARK_TEMPLATE(JacobiTheta1, cpp_bin_float_50);
BENCHMARK_TEMPLATE(JacobiTheta1, cpp_bin_float_100);

template<class Real>
void JacobiTheta1Tau(benchmark::State& state)
{
    std::random_device rd;
    std::mt19937_64 mt(rd());
    std::uniform_real_distribution<long double> unif(0,0.01);

    Real x = static_cast<Real>(unif(mt));
    Real q = static_cast<Real>(unif(mt));
    for (auto _ : state)
    {
        benchmark::DoNotOptimize(jacobi_theta1tau(x, q));
        x += std::numeric_limits<Real>::epsilon();
    }
}

BENCHMARK_TEMPLATE(JacobiTheta1Tau, float);
BENCHMARK_TEMPLATE(JacobiTheta1Tau, double);
BENCHMARK_TEMPLATE(JacobiTheta1Tau, long double);
#ifdef BOOST_HAS_FLOAT128
BENCHMARK_TEMPLATE(JacobiTheta1Tau, float128);
#endif
BENCHMARK_TEMPLATE(JacobiTheta1Tau, number<mpfr_float_backend<100>>);
BENCHMARK_TEMPLATE(JacobiTheta1Tau, number<mpfr_float_backend<200>>);
BENCHMARK_TEMPLATE(JacobiTheta1Tau, number<mpfr_float_backend<300>>);
BENCHMARK_TEMPLATE(JacobiTheta1Tau, number<mpfr_float_backend<400>>);
BENCHMARK_TEMPLATE(JacobiTheta1Tau, number<mpfr_float_backend<1000>>);
BENCHMARK_TEMPLATE(JacobiTheta1Tau, cpp_bin_float_50);
BENCHMARK_TEMPLATE(JacobiTheta1Tau, cpp_bin_float_100);

// The evaluation regimes, as (z, tau) pairs: the double-sided Gaussian sums
// used for tau < 1 and z != 0, the single-sided transformed series used for
// tau < 1 and z == 0, and the direct Fourier series used for tau >= 1. The
// argument index selects the pair.
static const double jacobi_theta_regimes[][2] = {
    { 0.5, 0.3 },  // Gaussian sums
    { 0.3, 0.05 }, // Gaussian sums, small tau
    { 0.0, 0.3 },  // z == 0 shortcut
    { 0.5, 3.0 },  // direct series
    { 0.0, 5.0 },  // direct series, one or two terms
};

struct theta1tau { template <class Real> Real operator()(Real z, Real tau) const { return boost::math::jacobi_theta1tau(z, tau); } };
struct theta2tau { template <class Real> Real operator()(Real z, Real tau) const { return boost::math::jacobi_theta2tau(z, tau); } };
struct theta3tau { template <class Real> Real operator()(Real z, Real tau) const { return boost::math::jacobi_theta3tau(z, tau); } };
struct theta4tau { template <class Real> Real operator()(Real z, Real tau) const { return boost::math::jacobi_theta4tau(z, tau); } };
struct theta4m1tau { template <class Real> Real operator()(Real z, Real tau) const { return boost::math::jacobi_theta4m1tau(z, tau); } };
// theta4(z, q) with q = exp(-pi tau), i.e. the q parameterization of the same point
struct theta4q { template <class Real> Real operator()(Real z, Real tau) const { return boost::math::jacobi_theta4(z, exp(-boost::math::constants::pi<Real>() * tau)); } };

template<class Real, class F>
void JacobiThetaRegime(benchmark::State& state)
{
    const double* regime = jacobi_theta_regimes[state.range(0)];
    Real z = static_cast<Real>(regime[0]);
    Real tau = static_cast<Real>(regime[1]);
    F f;
    for (auto _ : state)
    {
        benchmark::DoNotOptimize(f(z, tau));
        tau += std::numeric_limits<Real>::epsilon();
    }
}

#define JACOBI_THETA_REGIMES(F, Real) BENCHMARK_TEMPLATE(JacobiThetaRegime, Real, F)->DenseRange(0, 4)

JACOBI_THETA_REGIMES(theta1tau, float);
JACOBI_THETA_REGIMES(theta1tau, double);
JACOBI_THETA_REGIMES(theta2tau, double);
JACOBI_THETA_REGIMES(theta3tau, double);
JACOBI_THETA_REGIMES(theta4tau, float);
JACOBI_THETA_REGIMES(theta4tau, double);
JACOBI_THETA_REGIMES(theta4m1tau, double);
JACOBI_THETA_REGIMES(theta4q, double);
JACOBI_THETA_REGIMES(theta4tau, cpp_bin_float_50);

BENCHMARK_MAIN();
