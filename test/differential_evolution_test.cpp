/*
 * Copyright Nick Thompson, 2023
 * Use, modification and distribution are subject to the
 * Boost Software License, Version 1.0. (See accompanying file
 * LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
 */

#include "math_unit_test.hpp"
#include <random>
#include <boost/math/constants/constants.hpp>
#include <boost/math/tools/differential_evolution.hpp>
#ifdef BOOST_HAS_FLOAT128
#include <boost/multiprecision/float128.hpp>
using boost::multiprecision::float128;
#endif

using boost::math::tools::differential_evolution;
using boost::math::constants::two_pi;
using boost::math::constants::e;
using std::cbrt;
using std::sqrt;
using std::cos;
using std::exp;

// Taken from: https://en.wikipedia.org/wiki/Test_functions_for_optimization
template<typename Real>
Real ackley(std::array<Real, 2> const & v) {
    Real x = v[0];
    Real y = v[1];
    Real arg1 = -sqrt((x*x+y*y)/2)/5;
    Real arg2 = cos(two_pi<Real>()*x) + cos(two_pi<Real>()*y);
    return -20*exp(arg1) - exp(arg2/2) + 20 + e<Real>();
}


template<class Real>
void test_ackley()
{
    using ArgType = std::array<Real, 2>;
    std::vector<std::array<Real, 2>> bounds(2);
    bounds[0] = {-5, 5};
    bounds[1] = {-5, 5};
    auto de = differential_evolution(bounds, 0.9);
    std::random_device rd;
    std::mt19937_64 gen(rd());
    auto local_minima = de.template argmin<ArgType>(ackley<Real>, gen);
    CHECK_LE(std::abs(local_minima[0]), 10*std::numeric_limits<Real>::epsilon());
    CHECK_LE(std::abs(local_minima[1]), 10*std::numeric_limits<Real>::epsilon());
   
    // Works with a lambda:
    auto ack = [](std::array<Real, 2> const & x) { return ackley<Real>(x); };
    local_minima = de.template argmin<ArgType>(ack, gen);
    CHECK_LE(std::abs(local_minima[0]), 10*std::numeric_limits<Real>::epsilon());
    CHECK_LE(std::abs(local_minima[1]), 10*std::numeric_limits<Real>::epsilon());

    // Test that if an intial guess is the exact solution, the returned solution is the exact solution:
    std::array<Real, 2> initial_guess{0, 0};
    local_minima = de.template argmin<ArgType>(ack, gen,  &initial_guess);
    CHECK_EQUAL(local_minima[0], Real(0));
    CHECK_EQUAL(local_minima[1], Real(0));
}

template<typename Real>
auto rosenbrock_saddle(std::array<Real, 2> const & v) {
    auto x = v[0];
    auto y = v[1];
    return 100*(x*x - y)*(x*x - y) + (1 - x)*(1-x);
}

template<class Real>
void test_rosenbrock_saddle()
{
    using ArgType = std::array<Real, 2>;
    std::vector<ArgType> bounds(2);
    bounds[0] = {-2.048, 2.048};
    bounds[1] = {-2.048, 2.048};
    auto de = differential_evolution(bounds, 0.9);
    std::random_device rd;
    std::mt19937_64 gen(rd());
    auto local_minima = de.template argmin<ArgType>(rosenbrock_saddle<Real>, gen);
    CHECK_ABSOLUTE_ERROR(local_minima[0], Real(1), 10*std::numeric_limits<Real>::epsilon());
    CHECK_ABSOLUTE_ERROR(local_minima[1], Real(1), 10*std::numeric_limits<Real>::epsilon());

    // Does cancellation work?
    std::atomic<bool> cancel = true;
    gen.seed(12345);
    local_minima = de.template argmin<ArgType>(rosenbrock_saddle<Real>, gen,
                                               /*initial_guess*/  nullptr,
                                               /*target_value*/ std::numeric_limits<Real>::quiet_NaN(),
                                               &cancel);
    CHECK_GE(std::abs(local_minima[0] - Real(1)), 0.1);
} 

template<class Real>
Real rastrigin(std::vector<Real> const & v) {
   Real A = 10;
   Real y = 10*v.size();
   for (auto x : v) {
      y += x*x - A*cos(two_pi<Real>()*x);
   }
   return y;
}

template<class Real>
void test_rastrigin()
{
    using ArgType = std::vector<Real>;
    std::vector<std::array<Real,2>> bounds(8, {-5.12, 5.12});
    auto de = differential_evolution(bounds, 0.9);
    std::random_device rd;
    std::mt19937_64 gen(rd());
    auto local_minima = de.template argmin<ArgType>(rastrigin<Real>, gen);
    for (auto x : local_minima) {
       CHECK_ABSOLUTE_ERROR(x, Real(0), Real(2e-4));
    }

    // By definition, the value of the function which a target value is provided must be <= target_value.
    Real target_value = 1e-3;
    local_minima = de.template argmin<ArgType>(rastrigin<Real>, gen,
                                               /*initial guess*/ nullptr,
                                               target_value);
    CHECK_LE(rastrigin(local_minima), target_value);
} 

double sphere(std::vector<float> const & v) {
    double r = 0.0;
    for (auto x : v) {
        double x_ = x;
       r += x_*x_;
    }
    if (r >= 1) {
       return std::numeric_limits<double>::quiet_NaN();
    }
    return r;
}

// Tests NaN return types and return type != input type:
void test_sphere()
{
    using ArgType = std::vector<float>;
    std::vector<std::array<float,2>> bounds(8, {-2, 2});
    auto de = differential_evolution(bounds, 0.9);
    std::random_device rd;
    std::mt19937_64 gen(rd());
    auto local_minima = de.template argmin<ArgType>(sphere, gen);
    for (auto x : local_minima) {
       CHECK_ABSOLUTE_ERROR(x, 0.0, 2e-4);
    }
}

int main()
{
    test_ackley<float>();
    test_ackley<double>();
    test_rosenbrock_saddle<double>();
    test_rastrigin<float>();
    return boost::math::test::report_errors();
}
