#define BOOST_TEST_MODULE tanh_sinh_quadrature_test

#include <random>
#include <limits>
#include <functional>
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/type_index.hpp>
#include <boost/test/included/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/math/quadrature/detail/tanh_sinh_detail.hpp>
#include <boost/math/quadrature/tanh_sinh.hpp>
#include <boost/math/special_functions/sinc.hpp>
#include <boost/multiprecision/cpp_bin_float.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/math/special_functions/next.hpp>
#ifdef __GNUC__
#ifndef __clang__
#include <boost/multiprecision/float128.hpp>
using boost::multiprecision::float128;
#endif
#endif

using std::expm1;
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
using std::string;
using boost::multiprecision::cpp_bin_float_50;
using boost::multiprecision::cpp_bin_float_100;
using boost::multiprecision::cpp_dec_float_50;
using boost::multiprecision::cpp_dec_float_100;
using boost::multiprecision::cpp_bin_float;
using boost::math::tanh_sinh;
using boost::math::detail::tanh_sinh_detail;
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

template<typename Real>
Real g(Real t)
{
    return tanh(half_pi<Real>()*sinh(t));
}

template<typename Real>
Real g_prime(Real t)
{

    Real denom = cosh(half_pi<Real>()*sinh(t));
    return half_pi<Real>()*cosh(t)/(denom*denom);
}

void print_xt(cpp_dec_float_100 x, cpp_dec_float_100 t)
{
    std::cout << std::setprecision(std::numeric_limits<cpp_dec_float_100>::digits10 + 5) << std::fixed;
    std::cout << std::fixed << "            lexical_cast<Real>(\"" << x << "\"), // g(" << std::defaultfloat << t << ")\n";
}

void print_xt_prime(cpp_dec_float_100 x, cpp_dec_float_100 t)
{
    std::cout << std::setprecision(std::numeric_limits<cpp_dec_float_100>::digits10 + 5) << std::fixed;
    std::cout << std::fixed << "            lexical_cast<Real>(\"" << x << "\"), // g'(" << std::defaultfloat << t << ")\n";
}


void generate_constants()
{

    std::cout << "    static constexpr const std::vector<std::vector<Real>> abscissas{\n";
    std::cout << "        {\n";
    print_xt(g<cpp_dec_float_100>(0), 0);
    print_xt(g<cpp_dec_float_100>(1), 1);
    print_xt(g<cpp_dec_float_100>(2), 2);
    print_xt(g<cpp_dec_float_100>(3), 3);
    print_xt(g<cpp_dec_float_100>(4), 4);
    print_xt(g<cpp_dec_float_100>(5), 5);
    std::cout << "        },\n";
    size_t k = 1;
    while(k < 9)
    {
        cpp_dec_float_100 h =  (cpp_dec_float_100) 1/ (cpp_dec_float_100) (1 << k);

        std::cout << "        {\n";
        for(cpp_dec_float_100 t = h; t <= 5; t += 2*h)
        {
            auto x = g<cpp_dec_float_100>(t);
            print_xt(x, t);
        }
        std::cout << "        },\n";
        ++k;
    }
    std::cout << "};\n";


    std::cout << "    static constexpr const std::vector<std::vector<Real>> weights{\n";
    std::cout << "        {\n";
    print_xt_prime(g_prime<cpp_dec_float_100>(0), 0);
    print_xt_prime(g_prime<cpp_dec_float_100>(1), 1);
    print_xt_prime(g_prime<cpp_dec_float_100>(2), 2);
    print_xt_prime(g_prime<cpp_dec_float_100>(3), 3);
    print_xt_prime(g_prime<cpp_dec_float_100>(4), 4);
    print_xt_prime(g_prime<cpp_dec_float_100>(5), 5);
    std::cout << "        },\n";
    k = 1;
    while(k < 9)
    {
        cpp_dec_float_100 h =  (cpp_dec_float_100) 1/ (cpp_dec_float_100) (1 << k);

        std::cout << "        {\n";
        for(cpp_dec_float_100 t = h; t <= 5; t += 2*h)
        {
            auto x = g_prime<cpp_dec_float_100>(t);
            print_xt_prime(x, t);
        }
        std::cout << "        },\n";
        ++k;
    }
    std::cout << "    };\n";

}


template<class Real>
void test_detail()
{
    std::cout << "Testing tanh_sinh_detail on type " << boost::typeindex::type_id<Real>().pretty_name() << "\n";
    Real tol = 100*std::numeric_limits<Real>::epsilon();
    tanh_sinh_detail<Real> integrator(tol, 20);
    auto f = [](Real x) { return x*x; };
    Real Q = integrator.integrate(f);
    BOOST_CHECK_CLOSE(Q, 2*third<Real>(), 100*tol);
}


template<class Real>
void test_linear()
{
    std::cout << "Testing linear functions are integrated properly by tanh_sinh on type " << boost::typeindex::type_id<Real>().pretty_name() << "\n";
    Real tol = 100*std::numeric_limits<Real>::epsilon();
    tanh_sinh<Real> integrator(tol, 20);
    auto f = [](Real x) { return 5*x + 7; };
    Real error;
    Real Q = integrator.integrate(f, (Real) 0, (Real) 1, &error);
    BOOST_CHECK_CLOSE(Q, 9.5, 100*tol);
}


template<class Real>
void test_quadratic()
{
    std::cout << "Testing quadratic functions are integrated properly by tanh_sinh on type " << boost::typeindex::type_id<Real>().pretty_name() << "\n";
    Real tol = 100*std::numeric_limits<Real>::epsilon();
    tanh_sinh<Real> integrator(tol, 20);
    auto f = [](Real x) { return 5*x*x + 7*x + 12; };
    Real Q = integrator.integrate(f, (Real) 0, (Real) 1);
    BOOST_CHECK_CLOSE(Q, (Real) 17 + half<Real>()*third<Real>(), 100*tol);
}


template<class Real>
void test_singular()
{
    std::cout << "Testing integration of endpoint singularities on type " << boost::typeindex::type_id<Real>().pretty_name() << "\n";
    Real tol = 100*std::numeric_limits<Real>::epsilon();
    tanh_sinh<Real> integrator(tol, 20);
    auto f = [](Real x) { return log(x)*log(1-x); };
    Real Q = integrator.integrate(f, (Real) 0, (Real) 1);
    Real Q_expected = 2 - pi<Real>()*pi<Real>()*half<Real>()*third<Real>();

    BOOST_CHECK_CLOSE(Q, Q_expected, 100*tol);
}


// Examples taken from
//http://crd-legacy.lbl.gov/~dhbailey/dhbpapers/quadrature.pdf
template<class Real>
void test_ca()
{
    std::cout << "Testing integration of C(a) on type " << boost::typeindex::type_id<Real>().pretty_name() << "\n";
    Real tol = 100*std::numeric_limits<Real>::epsilon();
    Real error;

    tanh_sinh<Real> integrator(tol, 20);
    auto f1 = [](Real x) { return atan(x)/(x*(x*x + 1)) ; };
    Real Q = integrator.integrate(f1, (Real) 0, (Real) 1);
    Real Q_expected = pi<Real>()*ln_two<Real>()/8 + catalan<Real>()*half<Real>();
    BOOST_CHECK_CLOSE(Q, Q_expected, 100*tol);

    auto f2 = [](Real x) { Real t0 = x*x + 1; Real t1 = sqrt(t0); return atan(t1)/(t0*t1); };
    Q = integrator.integrate(f2, (Real) 0 , (Real) 1);
    Q_expected = pi<Real>()/4 - pi<Real>()/root_two<Real>() + 3*atan(root_two<Real>())/root_two<Real>();
    BOOST_CHECK_CLOSE(Q, Q_expected, 100*tol);

    auto f5 = [](Real t) { return t*t*log(t)/((t*t -1)*(t*t*t*t + 1)); };
    Q = integrator.integrate(f5, (Real) 0 , (Real) 1);
    Q_expected = pi<Real>()*pi<Real>()*(2 - root_two<Real>())/32;
    BOOST_CHECK_CLOSE(Q, Q_expected, 100*tol);


    // Oh it suffers on this one.
    auto f6 = [](Real t) { return log(t)*log(t); };
    Q = integrator.integrate(f6, (Real) 0 , (Real) 1);
    Q_expected = 2;
    BOOST_CHECK_CLOSE(Q, Q_expected, 5000*tol);


    // Although it doesn't get to the requested tolerance on this integral, the error bounds can be queried and are reasonable:
    auto f7 = [](Real t) { return sqrt(tan(t)); };
    Q = integrator.integrate(f7, (Real) 0 , (Real) half_pi<Real>(), &error);
    Q_expected = pi<Real>()/root_two<Real>();
    BOOST_CHECK_CLOSE(Q, Q_expected, 100*tol);

    auto f8 = [](Real t) { return log(cos(t)); };
    Q = integrator.integrate(f8, (Real) 0 , (Real) half_pi<Real>());
    Q_expected = -pi<Real>()*ln_two<Real>()*half<Real>();

    BOOST_CHECK_CLOSE(Q, Q_expected, 100*tol);
}


template<class Real>
void test_three_quadrature_schemes_examples()
{
    std::cout << "Testing integral in 'A Comparison of Three High Precision Quadrature Schemes' on type " << boost::typeindex::type_id<Real>().pretty_name() << "\n";
    Real tol = 100*std::numeric_limits<Real>::epsilon();
    Real Q;
    Real Q_expected;

    tanh_sinh<Real> integrator(tol, 20);
    // Example 1:
    auto f1 = [](Real t) { return t*log1p(t); };
    Q = integrator.integrate(f1, (Real) 0 , (Real) 1);
    Q_expected = half<Real>()*half<Real>();
    BOOST_CHECK_CLOSE(Q, Q_expected, 100*tol);


    // Example 2:
    auto f2 = [](Real t) { return t*t*atan(t); };
    Q = integrator.integrate(f2, (Real) 0 , (Real) 1);
    Q_expected = (pi<Real>() -2 + 2*ln_two<Real>())/12;
    BOOST_CHECK_CLOSE(Q, Q_expected, 100*tol);

    // Example 3:
    auto f3 = [](Real t) { return exp(t)*cos(t); };
    Q = integrator.integrate(f3, (Real) 0, half_pi<Real>());
    Q_expected = expm1(half_pi<Real>())*half<Real>();
    BOOST_CHECK_CLOSE(Q, Q_expected, 100*tol);

    // Example 4:
    auto f4 = [](Real x) { Real t0 = sqrt(x*x + 2); return atan(t0)/(t0*(x*x+1)); };
    Q = integrator.integrate(f4, (Real) 0 , (Real) 1);
    Q_expected = 5*pi<Real>()*pi<Real>()/96;
    BOOST_CHECK_CLOSE(Q, Q_expected, 100*tol);

    // Example 5:
    auto f5 = [](Real t) { return sqrt(t)*log(t); };
    Q = integrator.integrate(f5, (Real) 0 , (Real) 1);
    Q_expected = -4/ (Real) 9;
    BOOST_CHECK_CLOSE(Q, Q_expected, 100*tol);

    // Example 6:
    auto f6 = [](Real t) { return sqrt(1 - t*t); };
    Q = integrator.integrate(f6, (Real) 0 , (Real) 1);
    Q_expected = pi<Real>()/4;
    BOOST_CHECK_CLOSE(Q, Q_expected, 100*tol);
}


template<class Real>
void test_integration_over_real_line()
{
    std::cout << "Testing integrals over entire real line in 'A Comparison of Three High Precision Quadrature Schemes' on type " << boost::typeindex::type_id<Real>().pretty_name() << "\n";
    Real tol = 100*std::numeric_limits<Real>::epsilon();
    Real Q;
    Real Q_expected;
    Real error;
    tanh_sinh<Real> integrator(tol, 20);

    auto f1 = [](Real t) { return 1/(1+t*t);};
    Q = integrator.integrate(f1, -std::numeric_limits<Real>::infinity(), std::numeric_limits<Real>::infinity());
    Q_expected = pi<Real>();
    BOOST_CHECK_CLOSE(Q, Q_expected, 100*tol);

    auto f2 = [](Real t) { return exp(-t*t*half<Real>()); };
    Q = integrator.integrate(f2, -std::numeric_limits<Real>::infinity(), std::numeric_limits<Real>::infinity());
    Q_expected = root_two_pi<Real>();
    BOOST_CHECK_CLOSE(Q, Q_expected, 100*tol);

    // This test shows how oscillatory integrals are approximated very poorly by this method:
    std::cout << "Testing sinc integral: \n";
    Q = integrator.integrate(boost::math::sinc_pi<Real>, -std::numeric_limits<Real>::infinity(), std::numeric_limits<Real>::infinity(), &error);
    std::cout << "Error estimate of sinc integral: " << error << std::endl;
    Q_expected = pi<Real>();
    BOOST_CHECK_CLOSE(Q, Q_expected, 100*tol);

    auto f4 = [](Real t) { return 1/cosh(t);};
    Q = integrator.integrate(f4, -std::numeric_limits<Real>::infinity(), std::numeric_limits<Real>::infinity());
    Q_expected = pi<Real>();
    BOOST_CHECK_CLOSE(Q, Q_expected, 100*tol);

}

template<class Real>
void test_right_limit_infinite()
{
    std::cout << "Testing right limit infinite for tanh_sinh in 'A Comparison of Three High Precision Quadrature Schemes' on type " << boost::typeindex::type_id<Real>().pretty_name() << "\n";
    Real tol = sqrt(std::numeric_limits<Real>::epsilon());
    Real Q;
    Real Q_expected;
    tanh_sinh<Real> integrator;

    // Example 11:
    auto f1 = [](Real t) { return 1/(1+t*t);};
    Q = integrator.integrate(f1, 0, std::numeric_limits<Real>::infinity());
    Q_expected = half_pi<Real>();
    BOOST_CHECK_CLOSE(Q, Q_expected, 100*tol);

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
    std::cout << "Testing right limit infinite for tanh_sinh in 'A Comparison of Three High Precision Quadrature Schemes' on type " << boost::typeindex::type_id<Real>().pretty_name() << "\n";
    Real tol = sqrt(std::numeric_limits<Real>::epsilon());
    Real Q;
    Real Q_expected;
    tanh_sinh<Real> integrator;

    // Example 11:
    auto f1 = [](Real t) { return 1/(1+t*t);};
    Q = integrator.integrate(f1, -std::numeric_limits<Real>::infinity(), 0);
    Q_expected = half_pi<Real>();
    BOOST_CHECK_CLOSE(Q, Q_expected, 100*tol);
}


// A horrible function taken from
// http://www.chebfun.org/examples/quad/GaussClenCurt.html
template<class Real>
void test_horrible()
{
    Real tol = sqrt(std::numeric_limits<Real>::epsilon());
    Real Q;
    Real Q_expected;
    tanh_sinh<Real> integrator;

    auto f = [](Real x) { return x*sin(2*exp(2*sin(2*exp(2*x) ) ) ); };
    Q = integrator.integrate(f, -1, 1);
    // The example does not have more digits than this.
    // If this integration routine is correct, then all the digits here are correct.
    Q_expected = 0.336732834781728;
    BOOST_CHECK_CLOSE(Q, Q_expected, 100*std::numeric_limits<float>::epsilon());
}

BOOST_AUTO_TEST_CASE(tanh_sinh_quadrature_test)
{
    //generate_constants();

    test_right_limit_infinite<float>();
    test_right_limit_infinite<double>();
    test_right_limit_infinite<long double>();
    #ifdef __GNUC__
    #ifndef __clang__
    test_right_limit_infinite<float128>();
    #endif
    #endif


    test_left_limit_infinite<float>();
    test_left_limit_infinite<double>();
    test_left_limit_infinite<long double>();

    #ifdef __GNUC__
    #ifndef __clang__
    test_left_limit_infinite<float128>();
    #endif
    #endif


    test_detail<float>();
    test_detail<double>();
    test_detail<long double>();

    test_linear<float>();
    test_linear<double>();
    test_linear<long double>();

    #ifdef __GNUC__
    #ifndef __clang__
    test_linear<float128>();
    #endif
    #endif

    test_quadratic<float>();
    test_quadratic<double>();
    test_quadratic<long double>();

    #ifdef __GNUC__
    #ifndef __clang__
    test_quadratic<float128>();
    #endif
    #endif

    test_singular<float>();
    test_singular<double>();
    test_singular<long double>();

    test_ca<float>();
    test_ca<double>();
    test_ca<long double>();

    #ifdef __GNUC__
    #ifndef __clang__
    test_ca<float128>();
    #endif
    #endif

    test_three_quadrature_schemes_examples<float>();
    test_three_quadrature_schemes_examples<double>();
    test_three_quadrature_schemes_examples<long double>();

    #ifdef __GNUC__
    #ifndef __clang__
    test_three_quadrature_schemes_examples<float128>();
    #endif
    #endif


    test_integration_over_real_line<float>();
    test_integration_over_real_line<double>();
    test_integration_over_real_line<long double>();


    #ifdef __GNUC__
    #ifndef __clang__
    test_integration_over_real_line<float128>();
    #endif
    #endif

    test_horrible<float>();
    test_horrible<double>();
    test_horrible<long double>();

    #ifdef __GNUC__
    #ifndef __clang__
    test_horrible<float128>();
    #endif
    #endif


    test_singular<cpp_bin_float_50>();
    test_singular<cpp_bin_float_100>();

}
