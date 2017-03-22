#define BOOST_TEST_MODULE exp_sinh_quadrature_test

#include <random>
#include <limits>
#include <functional>
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/type_index.hpp>
#include <boost/test/included/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/math/quadrature/exp_sinh.hpp>
#include <boost/math/special_functions/sinc.hpp>
#include <boost/multiprecision/cpp_bin_float.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/math/special_functions/next.hpp>
#include <boost/math/special_functions/gamma.hpp>
#ifdef __GNUC__
#ifndef __clang__
#include <boost/multiprecision/float128.hpp>
using boost::multiprecision::float128;
#endif
#endif

using std::expm1;
using std::exp;
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
using boost::math::exp_sinh;
using boost::math::constants::pi;
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
void test_right_limit_infinite()
{
    std::cout << "Testing right limit infinite for tanh_sinh in 'A Comparison of Three High Precision Quadrature Schemes' on type " << boost::typeindex::type_id<Real>().pretty_name() << "\n";
    Real tol = sqrt(std::numeric_limits<Real>::epsilon());
    Real Q;
    Real Q_expected;
    exp_sinh<Real> integrator(tol, 12);

    // Example 12
    auto f2 = [](Real t) { return exp(-t)/sqrt(t); };
    Q = integrator.integrate(f2, 0, std::numeric_limits<Real>::infinity());
    Q_expected = root_pi<Real>();
    BOOST_CHECK_CLOSE(Q, Q_expected, 100*tol);

    auto f3 = [](Real t) { return exp(-t)*cos(t); };
    Q = integrator.integrate(f3, 0, std::numeric_limits<Real>::infinity());
    Q_expected = half<Real>();
    BOOST_CHECK_CLOSE(Q, Q_expected, 100*tol);

    auto f4 = [](Real t) { return 1/(1+t*t); };
    Q = integrator.integrate(f4, 1, std::numeric_limits<Real>::infinity());
    Q_expected = pi<Real>()/4;
    BOOST_CHECK_CLOSE(Q, Q_expected, 100*tol);
}

template<class Real>
void test_left_limit_infinite()
{
    std::cout << "Testing left limit infinite for tanh_sinh in 'A Comparison of Three High Precision Quadrature Schemes' on type " << boost::typeindex::type_id<Real>().pretty_name() << "\n";
    Real tol = sqrt(std::numeric_limits<Real>::epsilon());
    Real Q;
    Real Q_expected;
    exp_sinh<Real> integrator(tol, 8);

    // Example 11:
    auto f1 = [](Real t) { return 1/(1+t*t);};
    Q = integrator.integrate(f1, -std::numeric_limits<Real>::infinity(), 0);
    Q_expected = half_pi<Real>();
    BOOST_CHECK_CLOSE(Q, Q_expected, 100*tol);
}


// Some examples of tough integrals from NR, section 4.5.4:
template<class Real>
void test_nr_examples()
{
    using std::sin;
    using std::pow;
    using std::exp;
    using std::sqrt;
    std::cout << "Testing type " << boost::typeindex::type_id<Real>().pretty_name() << "\n";
    Real tol = sqrt(std::numeric_limits<Real>::epsilon());
    std::cout << std::setprecision(std::numeric_limits<Real>::digits10);
    Real Q;
    Real Q_expected;
    exp_sinh<Real> integrator(tol, 12);

    auto f0 = [](Real x) { return (Real) 0; };
    Q = integrator.integrate(f0);
    Q_expected = 0;
    BOOST_CHECK_CLOSE(Q, 0, 100*tol);

    auto f = [](Real x) { return 1/(1+x*x); };
    Q = integrator.integrate(f);
    Q_expected = half_pi<Real>();
    BOOST_CHECK_CLOSE(Q, Q_expected, 100*tol);

    auto f1 = [](Real x) { return sin(x*half<Real>())*pow(x, -3*half<Real>())*exp(-x); };
    Q = integrator.integrate(f1, 0, std::numeric_limits<Real>::infinity());
    Q_expected = sqrt(pi<Real>()*(sqrt((Real) 5) - 2));
    BOOST_CHECK_CLOSE(Q, Q_expected, 100*tol);

    auto f2 = [](Real x) { return pow(x, -(Real) 2/(Real) 7)*exp(-x*x); };
    Q = integrator.integrate(f2, 0, std::numeric_limits<Real>::infinity());
    Q_expected = half<Real>()*boost::math::tgamma((Real) 5/ (Real) 14);
    BOOST_CHECK_CLOSE(Q, Q_expected, 100*tol);

    auto f3 = [](Real x) { return (Real) 1/ (sqrt(x)*(1+x)); };
    Q = integrator.integrate(f3, 0, std::numeric_limits<Real>::infinity());
    Q_expected = pi<Real>();
    BOOST_CHECK_CLOSE(Q, Q_expected, 100*std::numeric_limits<float>::epsilon());

    auto f4 = [](Real t) { return exp(-t*t*half<Real>()); };
    Q = integrator.integrate(f4, 0, std::numeric_limits<Real>::infinity());
    Q_expected = root_two_pi<Real>()/2;
    BOOST_CHECK_CLOSE(Q, Q_expected, 100*tol);

    // This test shows how oscillatory integrals with 1/t decay are approximated very poorly by this method:
    //Q = integrator.integrate(boost::math::sinc_pi<Real>, 0, std::numeric_limits<Real>::infinity());
    //Q_expected = half_pi<Real>();
    //BOOST_CHECK_CLOSE(Q, Q_expected, 100*tol);

    auto f5 = [](Real t) { return 1/cosh(t);};
    Q = integrator.integrate(f5, 0, std::numeric_limits<Real>::infinity());
    Q_expected = half_pi<Real>();
    BOOST_CHECK_CLOSE(Q, Q_expected, 100*tol);
}


BOOST_AUTO_TEST_CASE(exp_sinh_quadrature_test)
{
    test_nr_examples<float>();
    test_nr_examples<double>();
    test_nr_examples<long double>();
    #ifdef __GNUC__
    #ifndef __clang__
    test_nr_examples<float128>();
    #endif
    #endif

    test_left_limit_infinite<double>();

    test_right_limit_infinite<float>();
    test_right_limit_infinite<double>();
    test_right_limit_infinite<long double>();
    #ifdef __GNUC__
    #ifndef __clang__
    test_right_limit_infinite<float128>();
    #endif
    #endif
}
