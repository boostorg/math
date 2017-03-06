#define BOOST_TEST_MODULE adaptive_trapezoidal

#include <boost/type_index.hpp>
#include <boost/test/included/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/math/quadrature/adaptive_trapezoidal.hpp>
#include <boost/multiprecision/cpp_bin_float.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>
#ifdef __GNUC__
#ifndef __clang__
#include <boost/multiprecision/float128.hpp>
#endif
#endif

using boost::multiprecision::cpp_bin_float_50;
using boost::multiprecision::cpp_bin_float_100;

template<class Real>
void test_constant()
{
    std::cout << "Testing constants are integrated correctly by the adaptive trapezoidal routine on type " << boost::typeindex::type_id<Real>().pretty_name()  << "\n";

    auto f = [](Real x) { return boost::math::constants::half<Real>(); };

    Real Q = boost::math::adaptive_trapezoidal<decltype(f), Real>(f, (Real) 0.0, (Real) 10.0);

    BOOST_CHECK_CLOSE(Q, 5.0, 100*std::numeric_limits<Real>::epsilon());
}


template<class Real>
void test_rational_periodic()
{
    using boost::math::constants::two_pi;
    using boost::math::constants::third;
    std::cout << "Testing that rational periodic functions are integrated correctly by trapezoidal rule on type " << boost::typeindex::type_id<Real>().pretty_name() << "\n";

    auto f = [](Real x) { return 1/(5 - 4*cos(x)); };

    Real tol = 100*std::numeric_limits<Real>::epsilon();
    Real Q = boost::math::adaptive_trapezoidal(f, (Real) 0.0, two_pi<Real>(), tol);

    BOOST_CHECK_CLOSE(Q, two_pi<Real>()*third<Real>(), 200*tol);
}

template<class Real>
void test_bump_function()
{
    std::cout << "Testing that bump functions are integrated correctly by trapezoidal rule on type " << boost::typeindex::type_id<Real>().pretty_name() << "\n";
    auto f = [](Real x) {
        if( x>= 1 || x <= -1)
        {
            return (Real) 0;
        }
        return (Real) exp(-(Real) 1/(1-x*x));
    };
    Real tol = 100*std::numeric_limits<Real>::epsilon();
    Real Q = boost::math::adaptive_trapezoidal(f, (Real) -1, (Real) 1, tol);
    // This is all the digits Mathematica gave me!
    BOOST_CHECK_CLOSE(Q, 0.443994, 1e-4);
}

template<class Real>
void test_zero_function()
{
    std::cout << "Testing that zero functions are integrated correctly by trapezoidal rule on type " << boost::typeindex::type_id<Real>().pretty_name() << "\n";
    auto f = [](Real x) { return (Real) 0;};
    Real tol = 100*std::numeric_limits<Real>::epsilon();
    Real Q = boost::math::adaptive_trapezoidal(f, (Real) -1, (Real) 1, tol);
    BOOST_CHECK_SMALL(Q, 100*tol);
}

template<class Real>
void test_sinsq()
{
    std::cout << "Testing that sin(x)^2 is integrated correctly by the trapezoidal rule on type " << boost::typeindex::type_id<Real>().pretty_name() << "\n";
    auto f = [](Real x) { return sin(10*x)*sin(10*x); };
    Real tol = 100*std::numeric_limits<Real>::epsilon();
    Real Q = boost::math::adaptive_trapezoidal(f, (Real) 0, (Real) boost::math::constants::pi<Real>(), tol);
    BOOST_CHECK_CLOSE(Q, boost::math::constants::half_pi<Real>(), 100*tol);

}

template<class Real>
void test_slowly_converging()
{
    std::cout << "Testing that non-periodic functions are correctly integrated by the trapezoidal rule, even if slowly, on type " << boost::typeindex::type_id<Real>().pretty_name() << "\n";
    // This function is not periodic, so it should not be fast to converge:
    auto f = [](Real x) { return sqrt(1 - x*x); };

    Real tol = sqrt(sqrt(std::numeric_limits<Real>::epsilon()));
    Real error_estimate;
    Real Q = boost::math::adaptive_trapezoidal(f, (Real) 0, (Real) 1, tol, 15, &error_estimate);
    BOOST_CHECK_CLOSE(Q, boost::math::constants::half_pi<Real>()/2, 1000*tol);
}

template<class Real>
void test_rational_sin()
{
    using std::pow;
    using boost::math::constants::two_pi;
    std::cout << "Testing that a rational sin function is integrated correctly by the trapezoidal rule on type " << boost::typeindex::type_id<Real>().pretty_name() << "\n";
    Real a = 5;
    auto f = [=](Real x) { Real t = a + sin(x); return 1.0/(t*t); };

    Real expected = two_pi<Real>()*a/pow(a*a - 1, 1.5);
    Real tol = 100*std::numeric_limits<Real>::epsilon();
    Real Q = boost::math::adaptive_trapezoidal(f, (Real) 0, (Real) boost::math::constants::two_pi<Real>(), tol);
    BOOST_CHECK_CLOSE(Q, expected, 100*tol);
}

BOOST_AUTO_TEST_CASE(adaptive_trapezoidal)
{
    test_constant<float>();
    test_constant<double>();
    test_constant<long double>();
    test_constant<cpp_bin_float_50>();
    test_constant<cpp_bin_float_100>();

    test_rational_periodic<float>();
    test_rational_periodic<double>();
    test_rational_periodic<long double>();
    test_rational_periodic<cpp_bin_float_50>();
    test_rational_periodic<cpp_bin_float_100>();

    test_bump_function<long double>();

    test_zero_function<float>();
    test_zero_function<double>();
    test_zero_function<long double>();
    test_zero_function<cpp_bin_float_50>();
    test_zero_function<cpp_bin_float_100>();

    test_sinsq<float>();
    test_sinsq<double>();
    test_sinsq<long double>();
    test_sinsq<cpp_bin_float_50>();
    test_sinsq<cpp_bin_float_100>();

    test_slowly_converging<float>();
    test_slowly_converging<double>();
    test_slowly_converging<long double>();

    test_rational_sin<float>();
    test_rational_sin<double>();
    test_rational_sin<long double>();

#ifdef __GNUC__
#ifndef __clang__
    test_constant<boost::multiprecision::float128>();
    test_rational_periodic<boost::multiprecision::float128>();
#endif
#endif
}
