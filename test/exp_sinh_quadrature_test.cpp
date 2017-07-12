// Copyright Nick Thompson, 2017
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

#define BOOST_TEST_MODULE exp_sinh_quadrature_test

#include <random>
#include <limits>
#include <functional>
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/math/concepts/real_concept.hpp>
#include <boost/type_index.hpp>
#include <boost/test/included/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/math/quadrature/exp_sinh.hpp>
#include <boost/math/special_functions/sinc.hpp>
#include <boost/math/special_functions/bessel.hpp>
#include <boost/multiprecision/cpp_bin_float.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/math/special_functions/next.hpp>
#include <boost/math/special_functions/gamma.hpp>

using std::exp;
using std::cos;
using std::tan;
using std::log;
using std::sqrt;
using std::abs;
using std::sinh;
using std::cosh;
using std::pow;
using boost::multiprecision::cpp_bin_float_50;
using boost::multiprecision::cpp_bin_float_100;
using boost::multiprecision::cpp_bin_float_quad;
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
using boost::math::quadrature::exp_sinh;

#if !defined(TEST1) && !defined(TEST2) && !defined(TEST3) && !defined(TEST4) && !defined(TEST5)
#  define TEST1
#  define TEST2
#  define TEST3
#  define TEST4
#  define TEST5
#endif

template<class Real>
void test_right_limit_infinite()
{
    std::cout << "Testing right limit infinite for tanh_sinh in 'A Comparison of Three High Precision Quadrature Schemes' on type " << boost::typeindex::type_id<Real>().pretty_name() << "\n";
    Real tol = sqrt(boost::math::tools::epsilon<Real>());
    Real Q;
    Real Q_expected;
    Real error;
    Real L1;
    exp_sinh<Real> integrator(tol, 12);

    // Example 12
    const auto f2 = [](const Real& t) { return exp(-t)/sqrt(t); };
    Q = integrator.integrate(f2, 0, std::numeric_limits<Real>::has_infinity ? std::numeric_limits<Real>::infinity() : boost::math::tools::max_value<Real>(), &error, &L1);
    Q_expected = root_pi<Real>();
    BOOST_CHECK_CLOSE_FRACTION(Q, Q_expected, tol);
    // The integrand is strictly positive, so it coincides with the value of the integral:
    BOOST_CHECK_CLOSE_FRACTION(L1, Q_expected, tol);

    auto f3 = [](Real t)->Real { Real z = exp(-t); if (z == 0) { return z; } return z*cos(t); };
    Q = integrator.integrate(f3, 0, std::numeric_limits<Real>::has_infinity ? std::numeric_limits<Real>::infinity() : boost::math::tools::max_value<Real>(), &error, &L1);
    Q_expected = half<Real>();
    BOOST_CHECK_CLOSE_FRACTION(Q, Q_expected, tol);

    auto f4 = [](Real t)->Real { return 1/(1+t*t); };
    Q = integrator.integrate(f4, 1, std::numeric_limits<Real>::has_infinity ? std::numeric_limits<Real>::infinity() : boost::math::tools::max_value<Real>(), &error, &L1);
    Q_expected = pi<Real>()/4;
    BOOST_CHECK_CLOSE_FRACTION(Q, Q_expected, tol);
    BOOST_CHECK_CLOSE_FRACTION(L1, Q_expected, tol);
}

template<class Real>
void test_left_limit_infinite()
{
    std::cout << "Testing left limit infinite for 1/(1+t^2) on type " << boost::typeindex::type_id<Real>().pretty_name() << "\n";
    Real tol = sqrt(boost::math::tools::epsilon<Real>());
    Real Q;
    Real Q_expected;
    Real error;
    Real L1;
    exp_sinh<Real> integrator(tol, 8);

    // Example 11:
    auto f1 = [](const Real& t) { return 1/(1+t*t);};
    Q = integrator.integrate(f1, std::numeric_limits<Real>::has_infinity ? -std::numeric_limits<Real>::infinity() : -boost::math::tools::max_value<Real>(), 0, &error, &L1);
    Q_expected = half_pi<Real>();
    BOOST_CHECK_CLOSE_FRACTION(Q, Q_expected, tol);
    BOOST_CHECK_CLOSE_FRACTION(L1, Q_expected, tol);
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
    Real tol = sqrt(boost::math::tools::epsilon<Real>());
    std::cout << std::setprecision(std::numeric_limits<Real>::digits10);
    Real Q;
    Real Q_expected;
    Real L1;
    Real error;
    exp_sinh<Real> integrator(tol, 12);

    auto f0 = [] (Real)->Real { return (Real) 0; };
    Q = integrator.integrate(f0, 0, std::numeric_limits<Real>::has_infinity ? std::numeric_limits<Real>::infinity() : boost::math::tools::max_value<Real>(), &error, &L1);
    Q_expected = 0;
    BOOST_CHECK_CLOSE_FRACTION(Q, 0.0f, tol);
    BOOST_CHECK_CLOSE_FRACTION(L1, 0.0f, tol);

    auto f = [](const Real& x) { return 1/(1+x*x); };
    Q = integrator.integrate(f, 0, std::numeric_limits<Real>::has_infinity ? std::numeric_limits<Real>::infinity() : boost::math::tools::max_value<Real>(), &error, &L1);
    Q_expected = half_pi<Real>();
    BOOST_CHECK_CLOSE_FRACTION(Q, Q_expected, tol);
    BOOST_CHECK_CLOSE_FRACTION(L1, Q_expected, tol);

    auto f1 = [](Real x)->Real {
        Real z1 = exp(-x);
        if (z1 == 0)
        {
            return (Real) 0;
        }
        Real z2 = pow(x, -3*half<Real>())*z1;
        if (z2 == 0)
        {
            return (Real) 0;
        }
        return sin(x*half<Real>())*z2;
    };

    Q = integrator.integrate(f1, 0, std::numeric_limits<Real>::has_infinity ? std::numeric_limits<Real>::infinity() : boost::math::tools::max_value<Real>(), &error, &L1);
    Q_expected = sqrt(pi<Real>()*(sqrt((Real) 5) - 2));

    // The integrand is oscillatory; the accuracy is low.
    BOOST_CHECK_CLOSE_FRACTION(Q, Q_expected, tol);

    auto f2 = [](Real x)->Real { return pow(x, -(Real) 2/(Real) 7)*exp(-x*x); };
    Q = integrator.integrate(f2, 0, std::numeric_limits<Real>::has_infinity ? std::numeric_limits<Real>::infinity() : boost::math::tools::max_value<Real>(), &error, &L1);
    Q_expected = half<Real>()*boost::math::tgamma((Real) 5/ (Real) 14);
    BOOST_CHECK_CLOSE_FRACTION(Q, Q_expected, tol);
    BOOST_CHECK_CLOSE_FRACTION(L1, Q_expected, tol);

    auto f3 = [](Real x)->Real { return (Real) 1/ (sqrt(x)*(1+x)); };
    Q = integrator.integrate(f3, 0, std::numeric_limits<Real>::has_infinity ? std::numeric_limits<Real>::infinity() : boost::math::tools::max_value<Real>(), &error, &L1);
    Q_expected = pi<Real>();

    BOOST_CHECK_CLOSE_FRACTION(Q, Q_expected, 10*boost::math::tools::epsilon<Real>());
    BOOST_CHECK_CLOSE_FRACTION(L1, Q_expected, 10*boost::math::tools::epsilon<Real>());

    auto f4 = [](const Real& t) { return exp(-t*t*half<Real>()); };
    Q = integrator.integrate(f4, 0, std::numeric_limits<Real>::has_infinity ? std::numeric_limits<Real>::infinity() : boost::math::tools::max_value<Real>(), &error, &L1);
    Q_expected = root_two_pi<Real>()/2;
    BOOST_CHECK_CLOSE_FRACTION(Q, Q_expected, tol);
    BOOST_CHECK_CLOSE_FRACTION(L1, Q_expected, tol);

    auto f5 = [](const Real& t) { return 1/cosh(t);};
    Q = integrator.integrate(f5, 0, std::numeric_limits<Real>::has_infinity ? std::numeric_limits<Real>::infinity() : boost::math::tools::max_value<Real>(), &error, &L1);
    Q_expected = half_pi<Real>();
    BOOST_CHECK_CLOSE_FRACTION(Q, Q_expected, tol);
    BOOST_CHECK_CLOSE_FRACTION(L1, Q_expected, tol);
}

// Definite integrals found in the CRC Handbook of Mathematical Formulas
template<class Real>
void test_crc()
{
    using std::sin;
    using std::pow;
    using std::exp;
    using std::sqrt;
    using std::log;
    std::cout << "Testing integral from CRC handbook on type " << boost::typeindex::type_id<Real>().pretty_name() << "\n";
    Real tol = sqrt(boost::math::tools::epsilon<Real>());
    std::cout << std::setprecision(std::numeric_limits<Real>::digits10);
    Real Q;
    Real Q_expected;
    Real L1;
    Real error;
    exp_sinh<Real> integrator(tol, 14);

    auto f0 = [](const Real& x) { return log(x)*exp(-x); };
    Q = integrator.integrate(f0, 0, std::numeric_limits<Real>::has_infinity ? std::numeric_limits<Real>::infinity() : boost::math::tools::max_value<Real>(), &error, &L1);
    Q_expected = -boost::math::constants::euler<Real>();
    BOOST_CHECK_CLOSE_FRACTION(Q, Q_expected, tol);

    // Test the integral representation of the gamma function:
    auto f1 = [](Real t)->Real { Real x = exp(-t);
        if(x == 0)
        {
            return (Real) 0;
        }
        return pow(t, (Real) 12 - 1)*x;
    };

    Q = integrator.integrate(f1, 0, std::numeric_limits<Real>::has_infinity ? std::numeric_limits<Real>::infinity() : boost::math::tools::max_value<Real>(), &error, &L1);
    Q_expected = boost::math::tgamma(12.0f);
    BOOST_CHECK_CLOSE_FRACTION(Q, Q_expected, tol);

    // Integral representation of the modified bessel function:
    // K_5(12)
    auto f2 = [](Real t)->Real {
        Real x = exp(-12*cosh(t));
        if (x == 0)
        {
            return (Real) 0;
        }
        return x*cosh(5*t);
    };
    Q = integrator.integrate(f2, 0, std::numeric_limits<Real>::has_infinity ? std::numeric_limits<Real>::infinity() : boost::math::tools::max_value<Real>(), &error, &L1);
    Q_expected = boost::math::cyl_bessel_k<int, Real>(5, 12);
    BOOST_CHECK_CLOSE_FRACTION(Q, Q_expected, tol);
    // Laplace transform of cos(at)
    Real a = 20;
    Real s = 1;
    auto f3 = [&](Real t)->Real {
        Real x = exp(-s*t);
        if (x == 0)
        {
            return (Real) 0;
        }
        return cos(a*t)*x;
    };

    // For high oscillation frequency, the quadrature sum is ill-conditioned.
    Q = integrator.integrate(f3, 0, std::numeric_limits<Real>::has_infinity ? std::numeric_limits<Real>::infinity() : boost::math::tools::max_value<Real>(), &error, &L1);
    Q_expected = s/(a*a+s*s);
    // Since the integrand is oscillatory, we increase the tolerance:
    BOOST_CHECK_CLOSE_FRACTION(Q, Q_expected, 100*tol);

    // Laplace transform of J_0(t):
    auto f4 = [&](Real t)->Real {
        Real x = exp(-s*t);
        if (x == 0)
        {
            return (Real) 0;
        }
        return boost::math::cyl_bessel_j(0, t)*x;
    };

    Q = integrator.integrate(f4, 0, std::numeric_limits<Real>::has_infinity ? std::numeric_limits<Real>::infinity() : boost::math::tools::max_value<Real>(), &error, &L1);
    Q_expected = 1/sqrt(1+s*s);
    BOOST_CHECK_CLOSE_FRACTION(Q, Q_expected, tol);

    auto f6 = [](const Real& t) { return exp(-t*t)*log(t);};
    Q = integrator.integrate(f6, 0, std::numeric_limits<Real>::has_infinity ? std::numeric_limits<Real>::infinity() : boost::math::tools::max_value<Real>(), &error, &L1);
    Q_expected = -boost::math::constants::root_pi<Real>()*(boost::math::constants::euler<Real>() + 2*ln_two<Real>())/4;
    BOOST_CHECK_CLOSE_FRACTION(Q, Q_expected, tol);
}


BOOST_AUTO_TEST_CASE(exp_sinh_quadrature_test)
{
#ifdef TEST1
    test_left_limit_infinite<float>();
    test_left_limit_infinite<double>();
    test_left_limit_infinite<long double>();
    test_right_limit_infinite<float>();
    test_right_limit_infinite<double>();
    test_right_limit_infinite<long double>();
    test_nr_examples<float>();
    test_nr_examples<double>();
    test_nr_examples<long double>();
    test_crc<float>();
    test_crc<double>();
    test_crc<long double>();
#endif
#ifdef TEST2
    test_left_limit_infinite<cpp_bin_float_quad>();
    test_right_limit_infinite<cpp_bin_float_quad>();
    test_nr_examples<cpp_bin_float_quad>();
    test_crc<cpp_bin_float_quad>();
#endif

#ifdef TEST3
    test_left_limit_infinite<boost::math::concepts::real_concept>();
    test_right_limit_infinite<boost::math::concepts::real_concept>();
    test_nr_examples<boost::math::concepts::real_concept>();
    test_crc<boost::math::concepts::real_concept>();
#endif

#ifdef TEST4
    test_left_limit_infinite<boost::multiprecision::cpp_dec_float_50>();
    test_right_limit_infinite<boost::multiprecision::cpp_dec_float_50>();
#endif
#ifdef TEST5
    test_nr_examples<boost::multiprecision::cpp_dec_float_50>();
    test_crc<boost::multiprecision::cpp_dec_float_50>();
#endif
}
