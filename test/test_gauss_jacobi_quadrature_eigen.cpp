// Copyright Jacob Hass, 2026
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

#include "math_unit_test.hpp"
#include "gauss_quadrature_test_helpers.hpp"
#include <complex>
#include <limits>
#include <Eigen/Dense>
#include <boost/math/quadrature/gauss_jacobi.hpp>

using boost::math::quadrature::gauss_jacobi;

template <class Matrix, unsigned N>
void test_matrix(const Matrix& zero)
{
   // Against 1 -> 16/15, x -> 0, x^2 -> 16/105, odd powers are 0.
   const Matrix Q = integrate_matrix(gauss_jacobi<double>(N, 2.0, 2.0), zero);
   CHECK_ULP_CLOSE(16.0 / 15.0, Q(0, 0).real(), 4);
   CHECK_ULP_CLOSE(16.0 / 105.0, Q(0, 0).imag(), 4);
   CHECK_ABSOLUTE_ERROR(Q(0, 1).imag(), 0.0, std::numeric_limits<double>::epsilon());
   CHECK_ABSOLUTE_ERROR(Q(0, 1).real(), 0.0, std::numeric_limits<double>::epsilon());
   CHECK_ULP_CLOSE(-16.0 / 105.0, Q(1, 0).real(), 4);
   CHECK_ABSOLUTE_ERROR(Q(1, 0).imag(), 0.0, std::numeric_limits<double>::epsilon());
   CHECK_ABSOLUTE_ERROR(Q(1, 1).real(), 0.0, std::numeric_limits<double>::epsilon());
   CHECK_ULP_CLOSE(-16.0 / 15.0, Q(1, 1).imag(), 4);
}

int main()
{
   using fixed = Eigen::Matrix<std::complex<double>, 2, 2>;
   test_matrix<fixed, 7>(fixed::Zero());
   test_matrix<fixed, 20>(fixed::Zero());
   test_matrix<Eigen::MatrixXcd, 15>(Eigen::MatrixXcd::Zero(2, 2));
   return boost::math::test::report_errors();
}
