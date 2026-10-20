// Copyright Nick Thompson, 2026
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt
#include <complex>
#include <cmath>
#include <limits>
#include <Eigen/Dense>
#include <boost/math/quadrature/gauss_kronrod.hpp>
#include "math_unit_test.hpp"

using boost::math::quadrature::gauss_kronrod;
using matrix = Eigen::MatrixXcd;

// Both parities of the embedded Gauss rule, with rectangular dynamic matrices.
template<unsigned N>
void test_adaptive()
{
    const matrix zero = matrix::Zero(2, 3);
    matrix value(2, 3);
    value << std::complex<double>(1, 2), 2., -3., 4., 5., std::complex<double>(-2, 1);
    const auto norm = [](const matrix& m) { return m.norm(); };
    unsigned calls = 0;
    const auto f = [&](double x) -> matrix {
        ++calls;
        return value / (1 + 10000*(x-0.37)*(x-0.37));
    };
    const double integral = (std::atan(100*(1-0.37)) + std::atan(100*0.37))/100;
    const matrix expected = value * integral;
    double error = 0, l1 = 0;
    matrix coarse = gauss_kronrod<double, N>::integrate(f, 0., 1., zero, norm, 0, 1e-10, &error, &l1);
    CHECK_EQUAL(calls, static_cast<unsigned>(N));
    const double coarse_error = (coarse - expected).norm();
    calls = 0;
    matrix result = gauss_kronrod<double, N>::integrate(f, 0., 1., zero, norm, 15, 1e-10, &error, &l1);
    CHECK_TRUE(calls > N);
    const unsigned refined_calls = calls;
    CHECK_TRUE((result - expected).norm() < coarse_error);
    CHECK_ABSOLUTE_ERROR(0, (result - expected).norm(), 1e-11);
    CHECK_ULP_CLOSE(value.norm()*integral, l1, 4);
    CHECK_GE(error, 0.0);
    CHECK_TRUE(error < 1e-10*expected.norm());

    // Rescaling a norm rescales error/L1, but preserves refinement decisions.
    const auto scaled_norm = [](const matrix& m) { return 8*m.norm(); };
    double scaled_error = 0, scaled_l1 = 0;
    calls = 0;
    matrix scaled = gauss_kronrod<double, N>::integrate(f, 0., 1., zero, scaled_norm, 15, 1e-10, &scaled_error, &scaled_l1);
    CHECK_EQUAL(calls, static_cast<unsigned>(refined_calls));
    CHECK_ABSOLUTE_ERROR(0, (scaled-result).norm(), 1e-14);
    CHECK_ULP_CLOSE(8*error, scaled_error, 4);
    CHECK_ULP_CLOSE(8*l1, scaled_l1, 4);
    result = gauss_kronrod<double, N>::integrate(f, 1., 0., zero, norm, 15, 1e-10, &error, &l1);
    CHECK_ABSOLUTE_ERROR(0, (result + expected).norm(), 1e-11);
    CHECK_ULP_CLOSE(value.norm()*integral, l1, 4);
}

void adaptive_matrix_integral()
{
    test_adaptive<15>();
    test_adaptive<21>();
}

void bounds_and_error_policies()
{
    using integrator = gauss_kronrod<double, 15>;
    const matrix zero = matrix::Zero(2, 3);
    const matrix value = matrix::Ones(2, 3);
    const auto norm = [](const matrix& m) { return m.norm(); };
    const double inf = std::numeric_limits<double>::infinity();
    const auto right = [&](double x) -> matrix { return value / ((1+x)*(1+x)); };
    const auto left = [&](double x) -> matrix { return value / ((1-x)*(1-x)); };
    const auto whole = [&](double x) -> matrix { return value / (boost::math::constants::pi<double>()*(1+x*x)); };
    double error = 0, l1 = 0;
    matrix result = integrator::integrate(right, 0., inf, zero, norm, 15, 1e-10, &error, &l1);
    CHECK_ABSOLUTE_ERROR(0, (result-value).norm(), 1e-12);
    CHECK_ULP_CLOSE(value.norm(), l1, 4);
    result = integrator::integrate(left, -inf, 0., zero, norm);
    CHECK_ABSOLUTE_ERROR(0, (result-value).norm(), 1e-12);
    result = integrator::integrate(whole, -inf, inf, zero, norm, 15, 1e-10);
    CHECK_ABSOLUTE_ERROR(0, (result-value).norm(), 1e-10);
    unsigned calls = 0;
    const auto counted = [&](double) -> matrix { ++calls; return value; };
    result = integrator::integrate(counted, 1., 1., zero, norm, 15, 1e-10, &error, &l1);
    CHECK_EQUAL(result.rows(), static_cast<Eigen::Index>(2));
    CHECK_EQUAL(result.cols(), static_cast<Eigen::Index>(3));
    CHECK_EQUAL(result.norm(), 0.0);
    CHECK_EQUAL(error, 0.0);
    CHECK_EQUAL(l1, 0.0);
    CHECK_EQUAL(calls, static_cast<unsigned>(0));
    const double nan = std::numeric_limits<double>::quiet_NaN();
    CHECK_THROW(integrator::integrate(counted, nan, 1., zero, norm), std::domain_error);
    using ignore = boost::math::policies::policy<boost::math::policies::domain_error<boost::math::policies::ignore_error>>;
    result = gauss_kronrod<double, 15, ignore>::integrate(counted, nan, 1., zero, norm, 15, 1e-10, &error, &l1);
    CHECK_EQUAL(result.rows(), static_cast<Eigen::Index>(2));
    CHECK_EQUAL(result.cols(), static_cast<Eigen::Index>(3));
    CHECK_NAN(error);
    CHECK_NAN(l1);
    for (Eigen::Index i = 0; i < result.size(); ++i)
        CHECK_NAN(result.data()[i].real());
    CHECK_EQUAL(calls, static_cast<unsigned>(0));
}

void scalar_compatibility_and_interval_scaling()
{
    using integrator = gauss_kronrod<double, 15>;
    const auto f = [](double x) { return std::exp(x); };
    const auto norm = [](double x) { return std::abs(x); };
    using function = decltype(f);
    double (*original)(function, double, double, unsigned, double, double*, double*) = &integrator::integrate<function>;
    CHECK_EQUAL(original(f, 0., 1., 0, 1e-10, nullptr, nullptr), integrator::integrate(f, 0., 1., 0., norm, 0, 1e-10));
    // An affine change of variables must scale the error estimate and L1.
    const auto stretched = [&](double x) { return f(x/4); };
    double error1 = 0, error2 = 0, l1 = 0, l2 = 0;
    double q1 = integrator::integrate(f, 0., 1., 0., norm, 0, 1e-10, &error1, &l1);
    double q2 = integrator::integrate(stretched, 0., 4., 0., norm, 0, 1e-10, &error2, &l2);
    CHECK_EQUAL(q2, 4*q1);
    CHECK_EQUAL(error2, 4*error1);
    CHECK_EQUAL(l2, 4*l1);
}

int main()
{
    adaptive_matrix_integral();
    bounds_and_error_policies();
    scalar_compatibility_and_interval_scaling();
    return boost::math::test::report_errors();
}
