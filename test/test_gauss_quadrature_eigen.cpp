// Copyright Nick Thompson, 2026
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <complex>
#include <limits>
#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>
#include <boost/math/quadrature/gauss.hpp>
#include "math_unit_test.hpp"

using boost::math::quadrature::gauss;
using matrix = Eigen::Matrix<std::complex<double>, 2, 2>;

// Regression for https://github.com/boostorg/math/issues/337.
// Return concrete Eigen matrices so the test does not rely on the lifetime
// of expression-template temporaries. Supply a shaped zero and a scalar norm.
template<unsigned N>
void test_matrix()
{
    const matrix zero = matrix::Zero();
    const auto norm = [](const matrix& m) { return m.norm(); };
    const auto f = [](double x) -> matrix {
        matrix result;
        result(0, 0) = {1, x * x};
        result(0, 1) = {x, 2 * x};
        result(1, 0) = {-x * x, x};
        result(1, 1) = {x * x * x, -1};
        return result;
    };
    const auto primitive = [](double x) -> matrix {
        matrix result;
        result(0, 0) = {x, x * x * x / 3};
        result(0, 1) = {x * x / 2, x * x};
        result(1, 0) = {-x * x * x / 3, x * x / 2};
        result(1, 1) = {x * x * x * x / 4, -x};
        return result;
    };
    const auto check = [](const matrix& actual, const matrix& expected) {
        for (unsigned i = 0; i < 2; ++i)
            for (unsigned j = 0; j < 2; ++j)
                CHECK_ABSOLUTE_ERROR(0, std::abs(actual(i, j) - expected(i, j)), 100 * std::numeric_limits<double>::epsilon());
    };

    check(gauss<double, N>::integrate(f, zero, norm), primitive(1) + -primitive(-1));
    check(gauss<double, N>::integrate(f, 0., 1., zero, norm), primitive(1));
    check(gauss<double, N>::integrate(f, 1., 0., zero, norm), -primitive(1));
    check(gauss<double, N>::integrate(f, 1., 1., zero, norm), matrix::Zero());

    // A constant matrix has an exactly known integral and L1 norm.
    const matrix value = f(1);
    double l1 = 0;
    check(gauss<double, N>::integrate([&](double) { return value; }, 0., 3., zero, norm, &l1),
          value * 3);
    CHECK_ULP_CLOSE(3 * value.norm(), l1, 4);
}

void matrix_valued_quadrature()
{
    test_matrix<7>();
    test_matrix<10>();
}

void dynamic_matrices_and_infinite_bounds()
{
    using dynamic_matrix = Eigen::MatrixXcd;
    const dynamic_matrix zero = dynamic_matrix::Zero(2, 3);
    dynamic_matrix value(2, 3);
    value << 1., 2., 3., 4., 5., 6.;
    const auto norm = [](const dynamic_matrix& m) { return m.norm(); };
    const auto check = [&](const dynamic_matrix& actual) {
        CHECK_EQUAL(actual.rows(), static_cast<Eigen::Index>(2));
        CHECK_EQUAL(actual.cols(), static_cast<Eigen::Index>(3));
        CHECK_ABSOLUTE_ERROR(0, (actual - value).norm(), 1e-10);
    };
    const double inf = std::numeric_limits<double>::infinity();
    using integrator = gauss<double, 30>;
    // These functions become constants under the half-line substitution.
    const auto right = [&](double x) -> dynamic_matrix { return value / ((1+x)*(1+x)); };
    const auto left = [&](double x) -> dynamic_matrix { return value / ((1-x)*(1-x)); };
    double l1 = 0;
    check(integrator::integrate(right, 0., inf, zero, norm, &l1));
    CHECK_ULP_CLOSE(value.norm(), l1, 4);
    check(integrator::integrate(left, -inf, 0., zero, norm));
    const auto whole = [&](double x) -> dynamic_matrix {
        return value / (boost::math::constants::pi<double>() * (1+x*x));
    };
    check(integrator::integrate(whole, -inf, inf, zero, norm));

    unsigned calls = 0;
    const auto counted = [&](double) -> dynamic_matrix { ++calls; return value; };
    l1 = -1;
    dynamic_matrix empty = integrator::integrate(counted, 2., 2., zero, norm, &l1);
    CHECK_EQUAL(empty.rows(), static_cast<Eigen::Index>(2));
    CHECK_EQUAL(empty.cols(), static_cast<Eigen::Index>(3));
    CHECK_EQUAL(empty.norm(), 0.0);
    CHECK_EQUAL(l1, 0.0);
    CHECK_EQUAL(calls, static_cast<unsigned>(0));

    const double nan = std::numeric_limits<double>::quiet_NaN();
    CHECK_THROW(integrator::integrate(counted, nan, 1., zero, norm), std::domain_error);
    using ignore = boost::math::policies::policy<
        boost::math::policies::domain_error<boost::math::policies::ignore_error>>;
    dynamic_matrix invalid = gauss<double, 30, ignore>::integrate(counted, nan, 1., zero, norm, &l1);
    CHECK_EQUAL(invalid.rows(), static_cast<Eigen::Index>(2));
    CHECK_EQUAL(invalid.cols(), static_cast<Eigen::Index>(3));
    CHECK_NAN(l1);
    for (Eigen::Index i = 0; i < invalid.size(); ++i)
        CHECK_NAN(invalid.data()[i].real());
    CHECK_EQUAL(calls, static_cast<unsigned>(0));
}

void scalar_overload_compatibility()
{
    const auto f = [](double x) { return x*x; };
    const auto norm = [](double x) { return std::abs(x); };
    using integrator = gauss<double, 7>;
    double old_l1 = 0, new_l1 = 0;
    CHECK_EQUAL(integrator::integrate(f, &old_l1), integrator::integrate(f, 0., norm, &new_l1));
    CHECK_EQUAL(old_l1, new_l1);
    CHECK_EQUAL(integrator::integrate(f, 0., 1., &old_l1), integrator::integrate(f, 0., 1., 0., norm, &new_l1));
    CHECK_EQUAL(old_l1, new_l1);
    using function = decltype(f);
    double (*original)(function, double*) = &integrator::integrate<function>;
    CHECK_EQUAL(original(f, nullptr), integrator::integrate(f));
}

void controllability_gramian()
{
    // The dimensions and matrix-exponential expression from issue #337.
    // A*A = 0, so exp(A*t) = I + A*t and the integral is polynomial.
    using gramian_matrix = Eigen::Matrix<std::complex<double>, 12, 12>;
    Eigen::Matrix<double, 12, 12> a = Eigen::Matrix<double, 12, 12>::Zero();
    Eigen::Matrix<double, 12, 4> b = Eigen::Matrix<double, 12, 4>::Zero();
    a(0, 1) = 1;
    b(1, 0) = 1;
    const auto f = [&](double t) -> gramian_matrix {
        return (a*t).exp() * b * b.transpose() * (a.transpose()*t).exp();
    };
    const gramian_matrix zero = gramian_matrix::Zero();
    const auto norm = [](const gramian_matrix& m) { return m.norm(); };
    gramian_matrix expected = zero;
    expected(0, 0) = 1./3;
    expected(0, 1) = expected(1, 0) = 0.5;
    expected(1, 1) = 1;
    gramian_matrix result = gauss<double, 10>::integrate(f, 0., 1., zero, norm);
    CHECK_ABSOLUTE_ERROR(0, (result - expected).norm(), 1e-14);

    // Verify that L1 follows the supplied norm, including on the canonical
    // interval and when the integration bounds are reversed.
    const auto scaled_norm = [](const gramian_matrix& m) { return 2*m.norm(); };
    const auto constant = [&](double) -> gramian_matrix { return expected; };
    double l1 = 0;
    result = gauss<double, 7>::integrate(constant, zero, scaled_norm, &l1);
    CHECK_ABSOLUTE_ERROR(0, (result - 2*expected).norm(), 1e-14);
    CHECK_ULP_CLOSE(4*expected.norm(), l1, 4);
    result = gauss<double, 10>::integrate(constant, 3., 0., zero, scaled_norm, &l1);
    CHECK_ABSOLUTE_ERROR(0, (result + 3*expected).norm(), 1e-14);
    CHECK_ULP_CLOSE(6*expected.norm(), l1, 4);
}

int main()
{
    matrix_valued_quadrature();
    dynamic_matrices_and_infinite_bounds();
    scalar_overload_compatibility();
    controllability_gramian();
    return boost::math::test::report_errors();
}
