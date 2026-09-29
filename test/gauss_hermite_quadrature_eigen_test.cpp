// Copyright Jacob Hass, 2026
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

#define BOOST_TEST_MODULE gauss_hermite_quadrature_eigen_test

#include <complex>
#include <limits>
#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>
#include <boost/test/included/unit_test.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/math/quadrature/hermite.hpp>

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

BOOST_AUTO_TEST_CASE(matrix_valued_quadrature)
{
    test_matrix<15>();
}
