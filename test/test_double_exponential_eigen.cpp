// Copyright Nick Thompson, 2026
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt
#define BOOST_TEST_MODULE double_exponential_eigen_test
#include <complex>
#include <cmath>
#include <limits>
#include <Eigen/Dense>
#include <boost/test/included/unit_test.hpp>
#include <boost/math/quadrature/exp_sinh.hpp>
#include <boost/math/quadrature/sinh_sinh.hpp>
#include <boost/math/quadrature/tanh_sinh.hpp>
#include <boost/math/quadrature/trapezoidal.hpp>

namespace q = boost::math::quadrature;
using matrix = Eigen::MatrixXcd;
const matrix zero = matrix::Zero(2, 3);
const matrix value = (matrix(2, 3) << std::complex<double>(1, 2), 2., -3., 4., 0., std::complex<double>(-2, 1)).finished();
const auto norm = [](const matrix& m) { return m.stableNorm(); };
const auto scaled_norm = [](const matrix& m) { return 8*m.stableNorm(); };
const double inf = std::numeric_limits<double>::infinity();
const double pi = boost::math::constants::pi<double>();

void check(const matrix& actual, const matrix& expected, double tolerance = 1e-9)
{
    BOOST_CHECK_EQUAL(actual.rows(), expected.rows());
    BOOST_CHECK_EQUAL(actual.cols(), expected.cols());
    BOOST_CHECK_SMALL((actual-expected).norm(), tolerance);
}

BOOST_AUTO_TEST_CASE(exp_sinh_matrices)
{
    q::exp_sinh<double> integrator;
    unsigned calls = 0;
    const auto f = [&](double x) -> matrix { ++calls; return value*std::exp(-std::sqrt(x)); };
    double error = 0, l1 = 0;
    std::size_t levels = 0;
    matrix coarse = integrator.integrate(f, zero, norm, 1e-3);
    const unsigned coarse_calls = calls;
    calls = 0;
    matrix result = integrator.integrate(f, zero, norm, 1e-11, &error, &l1, &levels);
    check(result, value*2);
    BOOST_CHECK_GT(calls, coarse_calls);
    BOOST_CHECK_GE(levels, 2);
    BOOST_CHECK_LT(error, 1e-9);
    BOOST_CHECK_CLOSE_FRACTION(l1, norm(value)*2, 1e-10);
    const unsigned refined_calls = calls;
    calls = 0;
    double scaled_error = 0, scaled_l1 = 0;
    check(integrator.integrate(f, zero, scaled_norm, 1e-11, &scaled_error, &scaled_l1), result);
    BOOST_CHECK_EQUAL(calls, refined_calls);
    BOOST_CHECK_CLOSE_FRACTION(scaled_error, 8*error, 1e-10);
    BOOST_CHECK_CLOSE_FRACTION(scaled_l1, 8*l1, 1e-10);
    check(integrator.integrate(f, 0., inf, zero, norm), value*2);
    const auto shifted = [&](double x) -> matrix { return value*std::exp(-(x-2)); };
    check(integrator.integrate(shifted, 2., inf, zero, norm), value);
    const auto left = [&](double x) -> matrix { return value*std::exp(x-2); };
    check(integrator.integrate(left, -inf, 2., zero, norm), value);
    // Distinct entrywise decay rates have independently known integrals.
    const auto mixed = [](double x) -> matrix {
        matrix m(2, 3);
        for (Eigen::Index i = 0; i < m.rows(); ++i)
            for (Eigen::Index j = 0; j < m.cols(); ++j)
                m(i, j) = value(i, j)*std::exp(-(1+i+2*j)*x);
        return m;
    };
    matrix expected(2, 3);
    for (Eigen::Index i = 0; i < expected.rows(); ++i)
        for (Eigen::Index j = 0; j < expected.cols(); ++j)
            expected(i, j) = value(i, j)/double(1+i+2*j);
    check(integrator.integrate(mixed, zero, norm), expected);
    BOOST_CHECK_THROW(integrator.integrate(f, 0., 1., zero, norm), std::domain_error);
}

BOOST_AUTO_TEST_CASE(sinh_sinh_matrices)
{
    q::sinh_sinh<double> integrator;
    unsigned calls = 0;
    const auto f = [&](double x) -> matrix { ++calls; return value*std::exp(-x*x); };
    double error = 0, l1 = 0;
    std::size_t levels = 0;
    integrator.integrate(f, zero, norm, 1e-3);
    const unsigned coarse_calls = calls;
    calls = 0;
    matrix result = integrator.integrate(f, zero, norm, 1e-11, &error, &l1, &levels);
    check(result, value*std::sqrt(pi));
    BOOST_CHECK_GT(calls, coarse_calls);
    BOOST_CHECK_GE(levels, 2);
    BOOST_CHECK_LT(error, 1e-9);
    BOOST_CHECK_CLOSE_FRACTION(l1, norm(value)*std::sqrt(pi), 1e-10);
    const unsigned refined_calls = calls;
    calls = 0;
    double scaled_error = 0, scaled_l1 = 0;
    check(integrator.integrate(f, zero, scaled_norm, 1e-11, &scaled_error, &scaled_l1), result);
    BOOST_CHECK_EQUAL(calls, refined_calls);
    BOOST_CHECK_CLOSE_FRACTION(scaled_error, 8*error, 1e-10);
    BOOST_CHECK_CLOSE_FRACTION(scaled_l1, 8*l1, 1e-10);
    const auto constant = [&](double) -> matrix { return value; };
    BOOST_CHECK_THROW(integrator.integrate(constant, zero, norm), std::domain_error);
}

BOOST_AUTO_TEST_CASE(tanh_sinh_matrices)
{
    q::tanh_sinh<double> integrator;
    unsigned calls = 0;
    // A narrow analytic peak, with a known antiderivative, needs refinement.
    const auto f = [&](double x) -> matrix { ++calls; return value/(1+10000*(x-0.37)*(x-0.37)); };
    const double exact = (std::atan(63.)+std::atan(37.))/100;
    integrator.integrate(f, 0., 1., zero, norm, 1e-3);
    const unsigned coarse_calls = calls;
    calls = 0;
    double error = 0, l1 = 0;
    std::size_t levels = 0;
    matrix result = integrator.integrate(f, 0., 1., zero, norm, 1e-11, &error, &l1, &levels);
    check(result, value*exact);
    BOOST_CHECK_GT(calls, coarse_calls);
    BOOST_CHECK_GT(levels, 4);
    BOOST_CHECK_LT(error, 1e-9);
    BOOST_CHECK_CLOSE_FRACTION(l1, norm(value)*exact, 1e-10);
    const unsigned refined_calls = calls;
    calls = 0;
    double scaled_error = 0, scaled_l1 = 0;
    check(integrator.integrate(f, 0., 1., zero, scaled_norm, 1e-11, &scaled_error, &scaled_l1), result);
    BOOST_CHECK_EQUAL(calls, refined_calls);
    BOOST_CHECK_CLOSE_FRACTION(scaled_error, 8*error, 1e-10);
    BOOST_CHECK_CLOSE_FRACTION(scaled_l1, 8*l1, 1e-10);
    check(integrator.integrate(f, 1., 0., zero, norm, 1e-11), -result);
    const auto smooth = [&](double x) -> matrix { return value*std::exp(x); };
    check(integrator.integrate(smooth, zero, norm), value*(std::exp(1.)-std::exp(-1.)));
    const auto gaussian = [&](double x) -> matrix { return value*std::exp(-x*x); };
    check(integrator.integrate(gaussian, -inf, inf, zero, norm), value*std::sqrt(pi));
    const auto right = [&](double x) -> matrix { return value*std::exp(-x); };
    const auto left = [&](double x) -> matrix { return value*std::exp(x); };
    check(integrator.integrate(right, 0., inf, zero, norm), value);
    check(integrator.integrate(left, -inf, 0., zero, norm), value);
    // Endpoint-distance functor: avoids cancellation in 1-x*x near +/-1.
    const auto endpoint = [&](double, double complement) -> matrix {
        const double c = std::abs(complement);
        return value/std::sqrt(c*(2-c));
    };
    check(integrator.integrate(endpoint, zero, norm), value*pi);
    check(integrator.integrate(endpoint, -1., 1., zero, norm), value*pi);
    check(integrator.integrate(endpoint, 1., -1., zero, norm), -value*pi);
    const auto shifted_endpoint = [&](double x, double complement) -> matrix {
        BOOST_CHECK_GE(x, 2.);
        BOOST_CHECK_LE(x, 6.); // x can round to the endpoint; complement retains distance.
        const double c = std::abs(complement);
        return value/std::sqrt(c*(4-c));
    };
    check(integrator.integrate(shifted_endpoint, 2., 6., zero, norm, 1e-11, &error, &l1), value*pi);
    BOOST_CHECK_CLOSE_FRACTION(l1, norm(value)*pi, 1e-10);
    // An endpoint singularity in the ordinary one-argument interface.
    const auto logarithm = [&](double x) -> matrix { return value*std::log(x); };
    check(integrator.integrate(logarithm, 0., 1., zero, norm), -value);
    calls = 0;
    check(integrator.integrate(f, 2., 2., zero, norm, 1e-11, &error, &l1, &levels), zero);
    BOOST_CHECK_EQUAL(calls, 0);
    BOOST_CHECK_EQUAL(error, 0);
    BOOST_CHECK_EQUAL(l1, 0);
    BOOST_CHECK_EQUAL(levels, 0);
    check(integrator.integrate(endpoint, 2., 2., zero, norm), zero);
}

BOOST_AUTO_TEST_CASE(trapezoidal_matrices)
{
    unsigned calls = 0;
    // Poisson kernel: its integral over one period is 2*pi.
    const double r = 0.8;
    const auto f = [&](double x) -> matrix { ++calls; return value*((1-r*r)/(1-2*r*std::cos(x)+r*r)); };
    q::trapezoidal(f, 0., 2*pi, zero, norm, 1e-3);
    const unsigned coarse_calls = calls;
    calls = 0;
    double error = 0, l1 = 0;
    matrix result = q::trapezoidal(f, 0., 2*pi, zero, norm, 1e-11, 12, &error, &l1);
    check(result, value*(2*pi));
    BOOST_CHECK_GT(calls, coarse_calls);
    BOOST_CHECK_LT(error, 1e-8);
    BOOST_CHECK_CLOSE_FRACTION(l1, norm(value)*(2*pi), 1e-10);
    const unsigned refined_calls = calls;
    calls = 0;
    double scaled_error = 0, scaled_l1 = 0;
    check(q::trapezoidal(f, 0., 2*pi, zero, scaled_norm, 1e-11, 12, &scaled_error, &scaled_l1), result);
    BOOST_CHECK_EQUAL(calls, refined_calls);
    BOOST_CHECK_CLOSE_FRACTION(scaled_error, 8*error, 1e-10);
    BOOST_CHECK_CLOSE_FRACTION(scaled_l1, 8*l1, 1e-10);
    check(q::trapezoidal(f, 2*pi, 0., zero, norm, 1e-11), -result);
    calls = 0;
    check(q::trapezoidal(f, 1., 1., zero, norm, 1e-11, 12, &error, &l1), zero);
    BOOST_CHECK_EQUAL(calls, 0);
    BOOST_CHECK_EQUAL(error, 0);
    BOOST_CHECK_EQUAL(l1, 0);
}

BOOST_AUTO_TEST_CASE(nonthrowing_errors)
{
    using ignore = boost::math::policies::policy<boost::math::policies::domain_error<boost::math::policies::ignore_error>>;
    const auto f = [&](double) -> matrix { return value; };
    const auto invalid = [](const matrix& m) {
        BOOST_CHECK_EQUAL(m.rows(), 2);
        BOOST_CHECK_EQUAL(m.cols(), 3);
        for (Eigen::Index i = 0; i < m.size(); ++i)
            BOOST_CHECK(std::isnan(m.data()[i].real()));
    };
    const double nan = std::numeric_limits<double>::quiet_NaN();
    invalid(q::exp_sinh<double, ignore>().integrate(f, 0., 1., zero, norm));
    invalid(q::exp_sinh<double, ignore>().integrate(f, zero, norm, 2.));
    invalid(q::sinh_sinh<double, ignore>().integrate(f, zero, norm));
    invalid(q::tanh_sinh<double, ignore>().integrate(f, nan, 1., zero, norm));
    invalid(q::trapezoidal(f, 0., inf, zero, norm, 1e-10, 12, static_cast<double*>(nullptr), static_cast<double*>(nullptr), ignore()));
}

BOOST_AUTO_TEST_CASE(evaluation_errors_and_scalar_compatibility)
{
    const auto bad = [](double x) -> matrix {
        if (std::abs(x) > 1e100) return zero;
        matrix m = value;
        m(0, 0) = std::numeric_limits<double>::quiet_NaN();
        return m;
    };
    BOOST_CHECK_THROW(q::exp_sinh<double>().integrate(bad, zero, norm), boost::math::evaluation_error);
    BOOST_CHECK_THROW(q::sinh_sinh<double>().integrate(bad, zero, norm), boost::math::evaluation_error);
    BOOST_CHECK_THROW(q::tanh_sinh<double>().integrate(bad, zero, norm), boost::math::evaluation_error);
    using ignore = boost::math::policies::policy<boost::math::policies::evaluation_error<boost::math::policies::ignore_error>>;
    // Non-throwing evaluation policies retain the concrete estimate, as in the
    // existing interface, rather than trying to construct a matrix from a scalar.
    matrix result = q::exp_sinh<double, ignore>().integrate(bad, zero, norm);
    BOOST_CHECK_EQUAL(result.rows(), 2);
    BOOST_CHECK_EQUAL(result.cols(), 3);
    BOOST_CHECK(std::isnan(result(0, 0).real()));
    result = q::sinh_sinh<double, ignore>().integrate(bad, zero, norm);
    BOOST_CHECK(std::isnan(result(0, 0).real()));
    result = q::tanh_sinh<double, ignore>().integrate(bad, zero, norm);
    BOOST_CHECK(std::isnan(result(0, 0).real()));

    const auto f = [](double x) { return std::exp(-x*x); };
    const auto scalar_norm = [](double x) { return std::abs(x); };
    using function = decltype(f);
    using exp_integrator = q::exp_sinh<double>;
    using sinh_integrator = q::sinh_sinh<double>;
    using tanh_integrator = q::tanh_sinh<double>;
    double (exp_integrator::*exp_original)(const function&, double, double*, double*, std::size_t*) const = &exp_integrator::integrate<function>;
    double (sinh_integrator::*sinh_original)(function, double, double*, double*, std::size_t*) const = &sinh_integrator::integrate<function>;
    double (tanh_integrator::*tanh_original)(function, double, double*, double*, std::size_t*) const = &tanh_integrator::integrate<function>;
    exp_integrator e;
    sinh_integrator s;
    tanh_integrator t;
    BOOST_CHECK_EQUAL((e.*exp_original)(f, 1e-10, nullptr, nullptr, nullptr), e.integrate(f, 0., scalar_norm, 1e-10));
    BOOST_CHECK_EQUAL((s.*sinh_original)(f, 1e-10, nullptr, nullptr, nullptr), s.integrate(f, 0., scalar_norm, 1e-10));
    BOOST_CHECK_EQUAL((t.*tanh_original)(f, 1e-10, nullptr, nullptr, nullptr), t.integrate(f, 0., scalar_norm, 1e-10));
    double (*trap_original)(function, double, double, double, std::size_t, double*, double*) = &q::trapezoidal<function, double>;
    BOOST_CHECK_EQUAL(trap_original(f, -1., 1., 1e-10, 12, nullptr, nullptr), q::trapezoidal(f, -1., 1., 0., scalar_norm, 1e-10));
}
