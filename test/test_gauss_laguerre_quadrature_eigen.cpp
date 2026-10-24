// Copyright Jacob Hass, 2026
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

#include "math_unit_test.hpp"
#include "gauss_quadrature_test_helpers.hpp"
#include <complex>
#include <Eigen/Dense>
#include <boost/math/quadrature/gauss_laguerre.hpp>

using boost::math::quadrature::gauss_laguerre;

template <class Matrix, unsigned N>
void test_matrix(const Matrix& zero)
{
   // Against exp(-x): 1 -> 1, e^(-x) * x -> 1, e^(-x) x^2 -> 2, e^(-x) x^3 -> 6.
   const Matrix Q = integrate_matrix(gauss_laguerre<double, N>(), zero);
   CHECK_ULP_CLOSE(1.0, Q(0, 0).real(), 4);
   CHECK_ULP_CLOSE(2.0, Q(0, 0).imag(), 32);
   CHECK_ULP_CLOSE(1.0, Q(0, 1).real(), 4);
   CHECK_ULP_CLOSE(2.0, Q(0, 1).imag(), 4);
   CHECK_ULP_CLOSE(-2.0, Q(1, 0).real(), 32);
   CHECK_ULP_CLOSE(Q(1, 0).imag(), 1.0, 4);
   CHECK_ULP_CLOSE(6.0, Q(1, 1).real(), 16);
   CHECK_ULP_CLOSE(-1.0, Q(1, 1).imag(), 4);
}

int main()
{
   using fixed = Eigen::Matrix<std::complex<double>, 2, 2>;
   test_matrix<fixed, 7>(fixed::Zero());
   test_matrix<fixed, 20>(fixed::Zero());
   test_matrix<Eigen::MatrixXcd, 15>(Eigen::MatrixXcd::Zero(2, 2));
   return boost::math::test::report_errors();
}
