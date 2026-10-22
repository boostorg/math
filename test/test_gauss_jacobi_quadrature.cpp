// Copyright Jacob Hass, 2026
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

#include "math_unit_test.hpp"
#include <cmath>
#include <complex>
#include <limits>
#include <valarray>
#include <boost/math/tools/test_value.hpp>
#include <boost/math/concepts/real_concept.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/math/quadrature/gauss_jacobi.hpp>
#include <boost/multiprecision/cpp_bin_float.hpp>
#include <boost/multiprecision/cpp_complex.hpp>
#include <boost/math/special_functions/factorials.hpp>

#ifdef BOOST_HAS_FLOAT128
#include <boost/multiprecision/float128.hpp>
#endif

using boost::math::quadrature::gauss_jacobi;

// Gauss-Jacobi integration is exact for polynomials of degree less than 2N - 1
// We use a polynomial of degree 10, so should be exact for N >= 6
template <class Real, unsigned N>
void test_moments(Real tol, const char* type_name)
{
   using std::pow;
   std::cout << "Testing Moments for N = " << N << " and type " << type_name << std::endl;
   auto f = [](const Real& x) { return pow(x, Real(10.0)); };
   Real L1;
   gauss_jacobi<Real> integrator = gauss_jacobi<Real>(N, -Real(0.5), -Real(0.25));
   Real Q = integrator.integrate(f, &L1);

   CHECK_ABSOLUTE_ERROR(Q, Real(0.604206393406616872500520394161000920557127865245), tol);
   // The integrand is nonnegative, so the L1 norm is the integral itself.
   CHECK_ABSOLUTE_ERROR(Q, L1, tol);
}

template <class Complex, unsigned N>
void test_complex(typename Complex::value_type tol)
{
   using Real = typename Complex::value_type;
   using std::cos;
   using std::exp;
   using std::sin;
   auto f = [](const Real& x) { return Complex(cos(x), sin(x)); };
   gauss_jacobi<Real> integrator = gauss_jacobi<Real>(N, 1.0, 1.0);
   Complex Q = integrator.integrate(f);
   Real real_expected = Real(4) * sin(Real(1)) - Real(4) * cos(Real(1));
   CHECK_ABSOLUTE_ERROR(real_expected, Q.real(), tol);
   CHECK_ABSOLUTE_ERROR(Real(0.0), Q.imag(), tol);
}

template <unsigned N>
void test_vector_valued()
{
   using vector = std::valarray<double>;
   const vector zero(0.0, 3);
   auto norm = [](const vector& v) { return std::sqrt((v * v).sum()); };
   auto f = [](double x) { return vector{1.0, x * x, std::cos(x)}; };
   double L1;
   gauss_jacobi<double> integrator = gauss_jacobi<double>(N, 1.0, 1.0);
   vector Q = integrator.integrate(f, zero, norm, &L1);
   CHECK_ULP_CLOSE(4.0 / 3.0, Q[0], 10);
   CHECK_ULP_CLOSE(4.0 / 15.0, Q[1], 20);
   CHECK_ULP_CLOSE(4.0 * sin(1.0) - 4.0 * cos(1.0), Q[2], 16);
   CHECK_LE(norm(Q), L1);
}

// No higher precision because it takes a long time because it has to
// compute the nodes/weights on the fly!
int main()
{
   test_moments<float, 6>(4 * std::numeric_limits<float>::epsilon(), "float");
   test_moments<double, 6>(4 * std::numeric_limits<double>::epsilon(), "double");
   test_moments<double, 10>(4 * std::numeric_limits<double>::epsilon(), "double");

#ifndef BOOST_MATH_NO_LONG_DOUBLE_MATH_FUNCTIONS
   test_moments<long double, 6>(4 * std::numeric_limits<double>::epsilon(), "long double");
#endif
#ifndef BOOST_MATH_NO_REAL_CONCEPT_TESTS
   test_moments<boost::math::concepts::real_concept, 6>(4e-16, "real concept");
#endif


   test_complex<std::complex<double>, 30>(4 * std::numeric_limits<double>::epsilon());

   test_vector_valued<30>();
   test_vector_valued<40>();

   return boost::math::test::report_errors();
}
