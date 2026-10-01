// Copyright Nick Thompson, 2026
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

#include "math_unit_test.hpp"
#include <complex>
#include <limits>
#include <Eigen/Dense>
#include <boost/math/constants/constants.hpp>
#include <boost/math/quadrature/gauss_hermite.hpp>

using boost::math::quadrature::gauss_hermite;

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
   // Against exp(-x^2): 1 -> sqrt(pi), x^2 -> sqrt(pi)/2, odd powers -> 0.
   const double root_pi = boost::math::constants::root_pi<double>();
   double L1;
   Matrix Q = gauss_hermite<double, N>::integrate(f, zero, norm, &L1);
   CHECK_EQUAL(Q.rows(), Eigen::Index(2));
   CHECK_EQUAL(Q.cols(), Eigen::Index(2));
   CHECK_ULP_CLOSE(root_pi, Q(0, 0).real(), 4);
   CHECK_ULP_CLOSE(root_pi / 2, Q(0, 0).imag(), 4);
   CHECK_EQUAL(Q(0, 1), std::complex<double>(0, 0));
   CHECK_ULP_CLOSE(-root_pi / 2, Q(1, 0).real(), 4);
   CHECK_EQUAL(Q(1, 0).imag(), 0.0);
   CHECK_EQUAL(Q(1, 1).real(), 0.0);
   CHECK_ULP_CLOSE(-root_pi, Q(1, 1).imag(), 4);
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
