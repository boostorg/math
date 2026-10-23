// Copyright Jacob Hass, 2026
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

#include "math_unit_test.hpp"
#include <complex>
#include <limits>
#include <Eigen/Dense>
#include <boost/math/quadrature/gauss_laguerre.hpp>

using boost::math::quadrature::gauss_laguerre;

// Return concrete Eigen matrices so the test does not rely on the lifetime
// of expression-template temporaries. Supply a shaped zero and a scalar norm.
template <class Matrix, unsigned N>
void test_matrix(const Matrix& zero)
{
   const auto norm = [](const Matrix& m) { return m.norm(); };
   const auto f = [&zero](double x) -> Matrix {
      Matrix result = zero;
      result(0, 0) = {1, x * x};
      result(0, 1) = {x, 2 * x};
      result(1, 0) = {-x * x, x};
      result(1, 1) = {x * x * x, -1};
      return result;
   };
   // Against exp(-x): 1 -> 1, e^(-x) * x -> 1, e^(-x) x^2 -> 2, e^(-x) x^3 -> 6.
   double L1;
   Matrix Q = gauss_laguerre<double, N>::integrate(f, zero, norm, nullptr, &L1);
   // The integrand is a polynomial that both rules integrate exactly, so the estimate is rounding.
   double error;
   gauss_laguerre<double, N>::integrate(f, zero, norm, &error, nullptr);
   CHECK_LE(error, 64 * std::numeric_limits<double>::epsilon() * L1);
   CHECK_EQUAL(Q.rows(), Eigen::Index(2));
   CHECK_EQUAL(Q.cols(), Eigen::Index(2));
   CHECK_ULP_CLOSE(1.0, Q(0, 0).real(), 4);
   CHECK_ULP_CLOSE(2.0, Q(0, 0).imag(), 16);
   CHECK_ULP_CLOSE(1.0, Q(0, 1).real(), 4);
   CHECK_ULP_CLOSE(2.0, Q(0, 1).imag(), 4);
   CHECK_ULP_CLOSE(-2.0, Q(1, 0).real(), 16);
   CHECK_ULP_CLOSE(Q(1, 0).imag(), 1.0, 4);
   CHECK_ULP_CLOSE(6.0, Q(1, 1).real(), 16);
   CHECK_ULP_CLOSE(-1.0, Q(1, 1).imag(), 4);
   CHECK_LE(norm(Q), L1);
}

int main()
{
   using fixed = Eigen::Matrix<std::complex<double>, 2, 2>;
   test_matrix<fixed, 7>(fixed::Zero());
   test_matrix<fixed, 20>(fixed::Zero());
   test_matrix<Eigen::MatrixXcd, 15>(Eigen::MatrixXcd::Zero(2, 2));
   return boost::math::test::report_errors();
}
