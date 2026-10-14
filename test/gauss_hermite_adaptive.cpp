// Copyright Jacob Hass, 2026
// Copyright Nick Thompson, 2026
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

#include "math_unit_test.hpp"
#include <cmath>
#include <complex>
#include <limits>
#include <valarray>
#include <boost/math/concepts/real_concept.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/math/quadrature/hermite_adaptive.hpp>
#include <boost/math/quadrature/gauss_hermite.hpp>
#include <boost/multiprecision/cpp_bin_float.hpp>
#include <boost/multiprecision/cpp_complex.hpp>

#ifdef BOOST_HAS_FLOAT128
#include <boost/multiprecision/float128.hpp>
#endif

using boost::math::quadrature::nested_hermite;
using boost::math::quadrature::gauss_hermite;

// Integral of exp(ix) exp(-x^2) over the real line is sqrt(pi) exp(-1/4).
template <class Complex, unsigned N>
void test_complex(typename Complex::value_type imaginary_tol)
{
   std::cout << "Testing complex valued function N = " << N << std::endl;
   using Real = typename Complex::value_type;
   using std::cos;
   using std::exp;
   using std::sin;
   auto f = [](const Real& x) { return Complex(cos(x), sin(x)); };
   Complex Q = nested_hermite<Real, N>::integrate(f);
   Real expected = boost::math::constants::root_pi<Real>() * exp(-Real(1) / 4);
   CHECK_ABSOLUTE_ERROR(expected, Q.real(), imaginary_tol);
   CHECK_ABSOLUTE_ERROR(Real(0), Q.imag(), imaginary_tol);

   Complex hermite = gauss_hermite<Real, (N-1)/2>::integrate(f);
   Real error = abs(hermite.real() - expected);
   Real nested_error = abs(Q.real() - expected);

   CHECK_LE(nested_error, error);
}

template <unsigned N>
void test_vector_valued(double tol)
{
   std::cout << "Testing vector valued N = " << N << std::endl;

   using vector = std::valarray<double>;
   const vector zero(0.0, 3);
   auto norm = [](const vector& v) { return std::sqrt((v * v).sum()); };
   auto f = [](double x) { return vector{1.0, x * x, std::cos(x)}; };
   double L1, error;
   vector Q = nested_hermite<double, N>::integrate(f, zero, norm, &error, &L1);
   const double root_pi = boost::math::constants::root_pi<double>();
   CHECK_ULP_CLOSE(root_pi, Q[0], 4);
   CHECK_ULP_CLOSE(root_pi / 2, Q[1], 4);
   CHECK_ABSOLUTE_ERROR(root_pi * std::exp(-0.25), Q[2], tol);
   CHECK_LE(norm(Q), L1);
}

int main()
{
   test_complex<std::complex<double>, 9>(1e-11);
   test_complex<std::complex<double>, 11>(10 * std::numeric_limits<double>::epsilon());
   test_complex<std::complex<double>, 35>(std::numeric_limits<double>::epsilon());

   test_complex<boost::multiprecision::cpp_complex_quad, 9>(1e-11);
   test_complex<boost::multiprecision::cpp_complex_quad, 11>(1e-14);
   test_complex<boost::multiprecision::cpp_complex_quad, 35>(1e-33);

   test_vector_valued<9>(1e-11);
   test_vector_valued<11>(10 * std::numeric_limits<double>::epsilon());
   test_vector_valued<35>(std::numeric_limits<double>::epsilon());

   return boost::math::test::report_errors();
}
