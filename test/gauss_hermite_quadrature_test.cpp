// Copyright Jacob Hass, 2026
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

#define BOOST_TEST_MODULE gauss_hermite_qudrature_test

#include <complex>

#include <boost/test/included/unit_test.hpp>
#include <boost/test/tools/floating_point_comparison.hpp>
#include <boost/math/tools/test_value.hpp>
#include <boost/multiprecision/cpp_bin_float.hpp>
#include <boost/multiprecision/cpp_complex.hpp>
#include <boost/math/quadrature/hermite.hpp>
#include <boost/multiprecision/complex128.hpp>
#include <boost/multiprecision/cpp_complex.hpp>

template<class Complex, unsigned Points, class RealType>
void test_complex_lambert_w(RealType tol)
{
   typedef typename Complex::value_type Real;
   auto lw = [](Real v)->Complex {
      using std::cos;
      using std::sin;
      Real sinv = sin(v);
      Real cosv = cos(v);

      Complex z{cosv, sinv};
      return z;
   };

   Complex Q = boost::math::quadrature::hermite<Real, Points>::integrate(lw);
   BOOST_CHECK_CLOSE_FRACTION(Q.real(), BOOST_MATH_TEST_VALUE(Real, 1.3803884470431429747734152467255912742707724655622107984502468507), tol);
   BOOST_CHECK_CLOSE_FRACTION(Q.imag(), BOOST_MATH_TEST_VALUE(Real, 0.0), tol);
}

BOOST_AUTO_TEST_CASE(gauss_hermite_qudrature_test)
{
   test_complex_lambert_w<std::complex<double>, 15>(2e-16);
   test_complex_lambert_w<boost::multiprecision::complex128, 15>(2e-25);
   test_complex_lambert_w<boost::multiprecision::cpp_complex_quad, 15>(2e-25);

   test_complex_lambert_w<std::complex<double>, 10>(2e-15);
   test_complex_lambert_w<boost::multiprecision::complex128, 10>(2e-15);
   test_complex_lambert_w<boost::multiprecision::cpp_complex_quad, 10>(2e-15);

   test_complex_lambert_w<std::complex<double>, 7>(1e-9);
   test_complex_lambert_w<boost::multiprecision::complex128, 7>(1e-9);
   test_complex_lambert_w<boost::multiprecision::cpp_complex_quad, 7>(1e-9);
}
