// Copyright Jacob Hass, 2026
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

#define BOOST_TEST_MODULE gauss_hermite_quadrature_eigen_test

#include <cmath>
#include <complex>
#include <limits>
#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>
#include <boost/test/included/unit_test.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/math/quadrature/hermite.hpp>
#include <boost/multiprecision/cpp_bin_float.hpp>

using boost::math::quadrature::hermite;
using matrix = Eigen::Matrix<std::complex<double>, 2, 2>;

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

    matrix expected;
    expected(0, 0) = {boost::math::constants::root_pi<double>(), boost::math::constants::root_pi<double>() * 0.5};
    expected(0, 1) = {0, 0};
    expected(1, 0) = {-boost::math::constants::root_pi<double>() * 0.5, 0};
    expected(1, 1) = {0, -boost::math::constants::root_pi<double>()};

    const auto check = [](const matrix& actual, const matrix& expected) {
        for (unsigned i = 0; i < 2; ++i)
            for (unsigned j = 0; j < 2; ++j)
                BOOST_CHECK_SMALL(std::abs(actual(i, j) - expected(i, j)),
                                  std::numeric_limits<double>::epsilon());
    };

    check(hermite<double, N>::integrate(f, zero, norm), expected);
}

// Check that hermite roots saved hermite.hpp are generated correctly with hermite_detail
template <typename Real, unsigned N>
void test_hermite_roots(double tol)
{
    using boost::math::quadrature::detail::hermite_detail;

    // Category 5 is not in hermite_detail. Thus, the weights and abscissa will be forced
    // to be regenerated
    hermite_detail<Real, 15, 5> integrator;
    const std::vector<Real> weights = integrator.weights();
    const std::vector<Real> vals = integrator.abscissa();

    hermite_detail<Real, 15, N> expected;
    const std::array<Real, 8> expected_vals = expected.abscissa();
    const std::array<Real, 8> expected_weights = expected.weights();

    for (unsigned i=0; i < vals.size(); i++)
    {
        BOOST_CHECK(abs(vals[i] - expected_vals[i]) < tol);
        BOOST_CHECK(abs(weights[i] - expected_weights[i]) < tol);
    }
}

template <typename Real>
void test_small_N(double tol)
{
    using boost::math::quadrature::detail::hermite_detail;

    hermite_detail<Real, 3, 5> integrator;
    const std::vector<Real> vals = integrator.abscissa();

    std::vector<Real> expected = {0, 0};
    expected[1] = sqrt(static_cast<Real>(3) / static_cast<Real>(2));

    for (unsigned i=0; i < vals.size(); i++)
    {
        std::cout << vals[i] << std::endl;
        std::cout << expected[i] << std::endl;
        BOOST_CHECK(abs(vals[i] - expected[i]) < tol);
    }
}

BOOST_AUTO_TEST_CASE(matrix_valued_quadrature)
{
    test_matrix<15>();

    typedef boost::multiprecision::number<boost::multiprecision::cpp_bin_float<250> > mp_type;
    test_hermite_roots<mp_type, 4>(1e-115);
    test_hermite_roots<double, 0>(1e-16);

    test_small_N<mp_type>(1e-115);
}
