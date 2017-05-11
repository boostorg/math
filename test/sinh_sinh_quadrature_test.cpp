// Copyright Nick Thompson, 2017
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

#define BOOST_TEST_MODULE sinh_sinh_quadrature_test

#include <random>
#include <limits>
#include <functional>
#include <boost/type_index.hpp>
#include <boost/test/included/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/math/quadrature/sinh_sinh.hpp>
#include <boost/math/special_functions/sinc.hpp>
#include <boost/multiprecision/cpp_bin_float.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>

using std::expm1;
using std::exp;
using std::sin;
using std::cos;
using std::atan;
using std::tan;
using std::log;
using std::log1p;
using std::asinh;
using std::atanh;
using std::sqrt;
using std::isnormal;
using std::abs;
using std::sinh;
using std::tanh;
using std::cosh;
using std::pow;
using std::string;
using boost::multiprecision::cpp_bin_float_quad;
using boost::math::sinh_sinh;
using boost::math::constants::pi;
using boost::math::constants::pi_sqr;
using boost::math::constants::half_pi;
using boost::math::constants::two_div_pi;
using boost::math::constants::half;
using boost::math::constants::third;
using boost::math::constants::half;
using boost::math::constants::third;
using boost::math::constants::catalan;
using boost::math::constants::ln_two;
using boost::math::constants::root_two;
using boost::math::constants::root_two_pi;
using boost::math::constants::root_pi;


template<class Real>
void test_nr_examples()
{
    std::cout << "Testing type " << boost::typeindex::type_id<Real>().pretty_name() << "\n";
    Real tol = sqrt(std::numeric_limits<Real>::epsilon());
    std::cout << std::setprecision(std::numeric_limits<Real>::digits10);
    Real Q;
    Real Q_expected;
    Real L1;
    Real error;
    sinh_sinh<Real> integrator(tol, 10);

    auto f0 = [](Real x) { return (Real) 0; };
    Q = integrator.integrate(f0, &error, &L1);
    Q_expected = 0;
    BOOST_CHECK_CLOSE(Q, 0, 100*tol);
    BOOST_CHECK_CLOSE(L1, 0, 100*tol);

    // In spite of the poles at \pm i, we still get a doubling of the correct digits at each level of refinement.
    auto f1 = [](Real t) { return 1/(1+t*t); };
    Q = integrator.integrate(f1, &error, &L1);
    Q_expected = pi<Real>();
    BOOST_CHECK_CLOSE(Q, Q_expected, 100*tol);
    BOOST_CHECK_CLOSE(L1, Q_expected, 100*tol);

    auto f2 = [](Real x) { return exp(-x*x); };
    Q = integrator.integrate(f2, &error, &L1);
    Q_expected = root_pi<Real>();
    BOOST_CHECK_CLOSE(Q, Q_expected, 100*tol);
    BOOST_CHECK_CLOSE(L1, Q_expected, 100*tol);

    auto f5 = [](Real t) { return 1/cosh(t);};
    Q = integrator.integrate(f5, &error, &L1);
    Q_expected = pi<Real>();
    BOOST_CHECK_CLOSE(Q, Q_expected, 100*tol);
    BOOST_CHECK_CLOSE(L1, Q_expected, 100*tol);

    // This oscillatory integral has rapid convergence because the oscillations get swamped by the exponential growth of the denominator.
    auto f8 = [](Real t) { return cos(t)/cosh(t);};
    Q = integrator.integrate(f8, &error, &L1);
    Q_expected = pi<Real>()/cosh(half_pi<Real>());
    BOOST_CHECK_CLOSE(Q, Q_expected, 100*tol);

}

// Test formulas for in the CRC Handbook of Mathematical functions, 32nd edition.
template<class Real>
void test_crc()
{
    std::cout << "Testing CRC formulas on type " << boost::typeindex::type_id<Real>().pretty_name() << "\n";
    Real tol = sqrt(std::numeric_limits<Real>::epsilon());
    std::cout << std::setprecision(std::numeric_limits<Real>::digits10);
    Real Q;
    Real Q_expected;
    Real L1;
    Real error;
    sinh_sinh<Real> integrator(tol, 10);

    // CRC Definite integral 698:
    auto f0 = [](Real x) {
      if(x == 0) {
        return (Real) 1;
      }
      return x/sinh(x);
    };
    Q = integrator.integrate(f0, &error, &L1);
    Q_expected = pi_sqr<Real>()/2;
    BOOST_CHECK_CLOSE(Q, Q_expected, 100*tol);
    BOOST_CHECK_CLOSE(L1, Q_expected, 100*tol);


    // CRC Definite integral 695:
    auto f1 = [](Real x) {
      if(x == 0) {
        return (Real) 1;
      }
      return (Real) sin(x)/sinh(x);
    };
    Q = integrator.integrate(f1, &error, &L1);
    Q_expected = pi<Real>()*tanh(half_pi<Real>());
    BOOST_CHECK_CLOSE(Q, Q_expected, 100*tol);
}


BOOST_AUTO_TEST_CASE(sinh_sinh_quadrature_test)
{
    test_nr_examples<float>();
    test_nr_examples<double>();
    test_nr_examples<long double>();
    test_nr_examples<cpp_bin_float_quad>();

    test_crc<float>();
    test_crc<double>();
    test_crc<long double>();
    test_crc<cpp_bin_float_quad>();

}
