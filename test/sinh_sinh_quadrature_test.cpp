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
#ifdef __GNUC__
#ifndef __clang__
#include <boost/multiprecision/float128.hpp>
using boost::multiprecision::float128;
#endif
#endif

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
using boost::multiprecision::cpp_bin_float_50;
using boost::multiprecision::cpp_bin_float_100;
using boost::multiprecision::cpp_dec_float_50;
using boost::multiprecision::cpp_dec_float_100;
using boost::multiprecision::cpp_bin_float;
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


    // This test shows how oscillatory integrals with 1/t decay are approximated very poorly by this method.
    // In fact |f(x)| <= C/(1+x^(1+eps)) for large x for this method to converge.
    //Q = integrator.integrate(boost::math::sinc_pi<Real>);
    //Q_expected = pi<Real>();
    //BOOST_CHECK_CLOSE(Q, Q_expected, 100*tol);

    auto f5 = [](Real t) { return 1/cosh(t);};
    Q = integrator.integrate(f5, &error, &L1);
    Q_expected = pi<Real>();
    BOOST_CHECK_CLOSE(Q, Q_expected, 100*tol);
    BOOST_CHECK_CLOSE(L1, Q_expected, 100*tol);

    // This demonstrates that oscillatory integrals which are nonetheless L1 integrable converge, but so slowly as to be useless.
    //std::cout << "cos(t)/(1+t*t)\n";
    //auto f6 = [](Real t) { return cos(t)/(1+t*t);};
    //Q = integrator.integrate(f6, &error, &L1);
    //Q_expected = pi<Real>()/boost::math::constants::e<Real>();
    //BOOST_CHECK_CLOSE(Q, Q_expected, 100*tol);
    //std::cout << "Exact            " << Q_expected << std::endl;
    //std::cout << "L1               " << L1 << std::endl;

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
    #ifdef __GNUC__
    #ifndef __clang__
    test_nr_examples<float128>();
    #endif
    #endif

    test_crc<float>();
    test_crc<double>();
    test_crc<long double>();
    #ifdef __GNUC__
    #ifndef __clang__
    test_crc<float128>();
    #endif
    #endif


    // Unfortunately this fails with message
    // boost/multiprecision/detail/functions/trig.hpp(69): fatal error: in "void boost::multiprecision::default_ops::hyp0F1(T&, const T&, const T&)
    // [with T = boost::multiprecision::backends::cpp_bin_float<50u>]": boost::exception_detail::clone_impl<boost::exception_detail::error_info_injector<std::runtime_error> >:
    // H0F1 Failed to Converge
    //test_nr_examples<cpp_bin_float_50>();
}
