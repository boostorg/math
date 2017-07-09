// Copyright Nick Thompson, 2017
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

#define BOOST_TEST_MODULE tanh_sinh_quadrature_test

#include <random>
#include <limits>
#include <functional>
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/math/concepts/real_concept.hpp>
#include <boost/type_index.hpp>
#include <boost/test/included/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/math/quadrature/detail/tanh_sinh_detail.hpp>
#include <boost/math/quadrature/tanh_sinh.hpp>
#include <boost/math/special_functions/sinc.hpp>
#include <boost/multiprecision/cpp_bin_float.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/math/special_functions/next.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/special_functions/beta.hpp>
#include <boost/math/special_functions/ellint_rc.hpp>
#include <boost/math/special_functions/ellint_rj.hpp>

#ifdef _MSC_VER
#pragma warning(disable:4127)  // Conditional expression is constant
#endif

#if !defined(TEST1) && !defined(TEST2) && !defined(TEST3) && !defined(TEST4) && !defined(TEST5) && !defined(TEST6) && !defined(TEST7)
#  define TEST1
#  define TEST2
#  define TEST3
#  define TEST4
#  define TEST5
#  define TEST6
#  define TEST7
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
using std::pow;
using std::string;
using boost::multiprecision::cpp_bin_float_50;
using boost::multiprecision::cpp_bin_float_100;
using boost::multiprecision::cpp_dec_float_50;
using boost::multiprecision::cpp_dec_float_100;
using boost::multiprecision::cpp_bin_float_quad;
using boost::math::quadrature::tanh_sinh;
using boost::math::quadrature::detail::tanh_sinh_detail;
using boost::math::constants::pi;
using boost::math::constants::half_pi;
using boost::math::constants::two_div_pi;
using boost::math::constants::two_pi;
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
    Real tol = sqrt(boost::math::tools::epsilon<Real>());
    tanh_sinh_detail<Real, boost::math::policies::policy<> > integrator(tol, 20);
    auto f = [](Real x, Real xc) { return x*x; };
    Real err;
    Real L1;
    Real Q = integrator.integrate(f, &err, &L1, 0);
    BOOST_CHECK_CLOSE_FRACTION(Q, 2*third<Real>(), tol);
    BOOST_CHECK_CLOSE_FRACTION(L1, 2*third<Real>(), tol);
}


template<class Real>
void test_linear()
{
    std::cout << "Testing linear functions are integrated properly by tanh_sinh on type " << boost::typeindex::type_id<Real>().pretty_name() << "\n";
    Real tol = 100*boost::math::tools::epsilon<Real>();
    tanh_sinh<Real> integrator(tol, 20);
    auto f = [](Real x) { return 5*x + 7; };
    Real error;
    Real L1;
    Real Q = integrator.integrate(f, (Real) 0, (Real) 1, &error, &L1);
    BOOST_CHECK_CLOSE_FRACTION(Q, 9.5, tol);
    BOOST_CHECK_CLOSE_FRACTION(L1, 9.5, tol);
}


template<class Real>
void test_quadratic()
{
    std::cout << "Testing quadratic functions are integrated properly by tanh_sinh on type " << boost::typeindex::type_id<Real>().pretty_name() << "\n";
    Real tol = 100*boost::math::tools::epsilon<Real>();
    tanh_sinh<Real> integrator(tol, 20);
    auto f = [](Real x) { return 5*x*x + 7*x + 12; };
    Real error;
    Real L1;
    Real Q = integrator.integrate(f, (Real) 0, (Real) 1, &error, &L1);
    BOOST_CHECK_CLOSE_FRACTION(Q, (Real) 17 + half<Real>()*third<Real>(), tol);
    BOOST_CHECK_CLOSE_FRACTION(L1, (Real) 17 + half<Real>()*third<Real>(), tol);
}


template<class Real>
void test_singular()
{
    std::cout << "Testing integration of endpoint singularities on type " << boost::typeindex::type_id<Real>().pretty_name() << "\n";
    Real tol = 100*boost::math::tools::epsilon<Real>();
    Real error;
    Real L1;
    tanh_sinh<Real> integrator(tol, 20);
    auto f = [](Real x) { return log(x)*log(1-x); };
    Real Q = integrator.integrate(f, (Real) 0, (Real) 1, &error, &L1);
    Real Q_expected = 2 - pi<Real>()*pi<Real>()*half<Real>()*third<Real>();

    BOOST_CHECK_CLOSE_FRACTION(Q, Q_expected, tol);
    BOOST_CHECK_CLOSE_FRACTION(L1, Q_expected, tol);
}


// Examples taken from
//http://crd-legacy.lbl.gov/~dhbailey/dhbpapers/quadrature.pdf
template<class Real>
void test_ca()
{
    std::cout << "Testing integration of C(a) on type " << boost::typeindex::type_id<Real>().pretty_name() << "\n";
    Real tol = sqrt(boost::math::tools::epsilon<Real>());
    Real error;
    Real L1;

    tanh_sinh<Real> integrator(tol, 20);
    auto f1 = [](Real x) { return atan(x)/(x*(x*x + 1)) ; };
    Real Q = integrator.integrate(f1, (Real) 0, (Real) 1, &error, &L1);
    Real Q_expected = pi<Real>()*ln_two<Real>()/8 + catalan<Real>()*half<Real>();
    BOOST_CHECK_CLOSE_FRACTION(Q, Q_expected, tol);
    BOOST_CHECK_CLOSE_FRACTION(L1, Q_expected, tol);

    auto f2 = [](Real x) { Real t0 = x*x + 1; Real t1 = sqrt(t0); return atan(t1)/(t0*t1); };
    Q = integrator.integrate(f2, (Real) 0 , (Real) 1, &error, &L1);
    Q_expected = pi<Real>()/4 - pi<Real>()/root_two<Real>() + 3*atan(root_two<Real>())/root_two<Real>();
    BOOST_CHECK_CLOSE_FRACTION(Q, Q_expected, tol);
    BOOST_CHECK_CLOSE_FRACTION(L1, Q_expected, tol);

    auto f5 = [](Real t) { return t*t*log(t)/((t*t - 1)*(t*t*t*t + 1)); };
    Q = integrator.integrate(f5, (Real) 0 , (Real) 1);
    Q_expected = pi<Real>()*pi<Real>()*(2 - root_two<Real>())/32;
    BOOST_CHECK_CLOSE_FRACTION(Q, Q_expected, tol);


    // Oh it suffers on this one.
    auto f6 = [](Real t) { Real x = log(t); return x*x; };
    Q = integrator.integrate(f6, (Real) 0 , (Real) 1, &error, &L1);
    Q_expected = 2;

    BOOST_CHECK_CLOSE_FRACTION(Q, Q_expected, 50*tol);
    BOOST_CHECK_CLOSE_FRACTION(L1, Q_expected, 50*tol);


    // Although it doesn't get to the requested tolerance on this integral, the error bounds can be queried and are reasonable:
    auto f7 = [](Real t) { return sqrt(tan(t)); };
    Q = integrator.integrate(f7, (Real) 0 , (Real) half_pi<Real>(), &error, &L1);
    Q_expected = pi<Real>()/root_two<Real>();
    BOOST_CHECK_CLOSE_FRACTION(Q, Q_expected, tol);
    BOOST_CHECK_CLOSE_FRACTION(L1, Q_expected, tol);

    auto f8 = [](Real t) { return log(cos(t)); };
    Q = integrator.integrate(f8, (Real) 0 , half_pi<Real>(), &error, &L1);
    Q_expected = -pi<Real>()*ln_two<Real>()*half<Real>();
    BOOST_CHECK_CLOSE_FRACTION(Q, Q_expected, tol);
    BOOST_CHECK_CLOSE_FRACTION(L1, -Q_expected, tol);
}


template<class Real>
void test_three_quadrature_schemes_examples()
{
    std::cout << "Testing integral in 'A Comparison of Three High Precision Quadrature Schemes' on type " << boost::typeindex::type_id<Real>().pretty_name() << "\n";
    Real tol = sqrt(boost::math::tools::epsilon<Real>());
    Real Q;
    Real Q_expected;

    tanh_sinh<Real> integrator(tol, 20);
    // Example 1:
    auto f1 = [](Real t) { return t*boost::math::log1p(t); };
    Q = integrator.integrate(f1, (Real) 0 , (Real) 1);
    Q_expected = half<Real>()*half<Real>();
    BOOST_CHECK_CLOSE_FRACTION(Q, Q_expected, tol);


    // Example 2:
    auto f2 = [](Real t) { return t*t*atan(t); };
    Q = integrator.integrate(f2, (Real) 0 , (Real) 1);
    Q_expected = (pi<Real>() -2 + 2*ln_two<Real>())/12;
    BOOST_CHECK_CLOSE_FRACTION(Q, Q_expected, tol);

    // Example 3:
    auto f3 = [](Real t) { return exp(t)*cos(t); };
    Q = integrator.integrate(f3, (Real) 0, half_pi<Real>());
    Q_expected = boost::math::expm1(half_pi<Real>())*half<Real>();
    BOOST_CHECK_CLOSE_FRACTION(Q, Q_expected, tol);

    // Example 4:
    auto f4 = [](Real x) { Real t0 = sqrt(x*x + 2); return atan(t0)/(t0*(x*x+1)); };
    Q = integrator.integrate(f4, (Real) 0 , (Real) 1);
    Q_expected = 5*pi<Real>()*pi<Real>()/96;
    BOOST_CHECK_CLOSE_FRACTION(Q, Q_expected, tol);

    // Example 5:
    auto f5 = [](Real t) { return sqrt(t)*log(t); };
    Q = integrator.integrate(f5, (Real) 0 , (Real) 1);
    Q_expected = -4/ (Real) 9;
    BOOST_CHECK_CLOSE_FRACTION(Q, Q_expected, tol);

    // Example 6:
    auto f6 = [](Real t) { return sqrt(1 - t*t); };
    Q = integrator.integrate(f6, (Real) 0 , (Real) 1);
    Q_expected = pi<Real>()/4;
    BOOST_CHECK_CLOSE_FRACTION(Q, Q_expected, tol);
}


template<class Real>
void test_integration_over_real_line()
{
    std::cout << "Testing integrals over entire real line in 'A Comparison of Three High Precision Quadrature Schemes' on type " << boost::typeindex::type_id<Real>().pretty_name() << "\n";
    Real tol = sqrt(boost::math::tools::epsilon<Real>());
    Real Q;
    Real Q_expected;
    Real error;
    Real L1;
    tanh_sinh<Real> integrator(tol, 20);

    auto f1 = [](Real t) { return 1/(1+t*t);};
    Q = integrator.integrate(f1, std::numeric_limits<Real>::has_infinity ? -std::numeric_limits<Real>::infinity() : -boost::math::tools::max_value<Real>(), std::numeric_limits<Real>::has_infinity ? std::numeric_limits<Real>::infinity() : boost::math::tools::max_value<Real>(), &error, &L1);
    Q_expected = pi<Real>();
    BOOST_CHECK_CLOSE_FRACTION(Q, Q_expected, tol);
    BOOST_CHECK_CLOSE_FRACTION(L1, Q_expected, tol);

    auto f2 = [](Real t) { return exp(-t*t*half<Real>()); };
    Q = integrator.integrate(f2, std::numeric_limits<Real>::has_infinity ? -std::numeric_limits<Real>::infinity() : -boost::math::tools::max_value<Real>(), std::numeric_limits<Real>::has_infinity ? std::numeric_limits<Real>::infinity() : boost::math::tools::max_value<Real>(), &error, &L1);
    Q_expected = root_two_pi<Real>();
    BOOST_CHECK_CLOSE_FRACTION(Q, Q_expected, tol);
    BOOST_CHECK_CLOSE_FRACTION(L1, Q_expected, tol);

    // This test shows how oscillatory integrals are approximated very poorly by this method:
    //std::cout << "Testing sinc integral: \n";
    //Q = integrator.integrate(boost::math::sinc_pi<Real>, -std::numeric_limits<Real>::infinity(), std::numeric_limits<Real>::infinity(), &error, &L1);
    //std::cout << "Error estimate of sinc integral: " << error << std::endl;
    //std::cout << "L1 norm of sinc integral " << L1 << std::endl;
    //Q_expected = pi<Real>();
    //BOOST_CHECK_CLOSE(Q, Q_expected, 100*tol);

    auto f4 = [](Real t) { return 1/cosh(t);};
    Q = integrator.integrate(f4, std::numeric_limits<Real>::has_infinity ? -std::numeric_limits<Real>::infinity() : -boost::math::tools::max_value<Real>(), std::numeric_limits<Real>::has_infinity ? std::numeric_limits<Real>::infinity() : boost::math::tools::max_value<Real>(), &error, &L1);
    Q_expected = pi<Real>();
    BOOST_CHECK_CLOSE_FRACTION(Q, Q_expected, tol);
    BOOST_CHECK_CLOSE_FRACTION(L1, Q_expected, tol);

}

template<class Real>
void test_right_limit_infinite()
{
    std::cout << "Testing right limit infinite for tanh_sinh in 'A Comparison of Three High Precision Quadrature Schemes' on type " << boost::typeindex::type_id<Real>().pretty_name() << "\n";
    Real tol = sqrt(boost::math::tools::epsilon<Real>());
    Real Q;
    Real Q_expected;
    Real error;
    Real L1;
    tanh_sinh<Real> integrator;

    // Example 11:
    auto f1 = [](Real t) { return 1/(1+t*t);};
    Q = integrator.integrate(f1, 0, std::numeric_limits<Real>::has_infinity ? std::numeric_limits<Real>::infinity() : boost::math::tools::max_value<Real>(), &error, &L1);
    Q_expected = half_pi<Real>();
    BOOST_CHECK_CLOSE(Q, Q_expected, 100*tol);

    // Example 12
    auto f2 = [](Real t) { return exp(-t)/sqrt(t); };
    Q = integrator.integrate(f2, 0, std::numeric_limits<Real>::has_infinity ? std::numeric_limits<Real>::infinity() : boost::math::tools::max_value<Real>(), &error, &L1);
    Q_expected = root_pi<Real>();
    BOOST_CHECK_CLOSE(Q, Q_expected, 1000*tol);

    auto f3 = [](Real t) { return exp(-t)*cos(t); };
    Q = integrator.integrate(f3, 0, std::numeric_limits<Real>::has_infinity ? std::numeric_limits<Real>::infinity() : boost::math::tools::max_value<Real>(), &error, &L1);
    Q_expected = half<Real>();
    BOOST_CHECK_CLOSE(Q, Q_expected, 100*tol);

    auto f4 = [](Real t) { return 1/(1+t*t); };
    Q = integrator.integrate(f4, 1, std::numeric_limits<Real>::has_infinity ? std::numeric_limits<Real>::infinity() : boost::math::tools::max_value<Real>(), &error, &L1);
    Q_expected = pi<Real>()/4;
    BOOST_CHECK_CLOSE(Q, Q_expected, 100*tol);
}

template<class Real>
void test_left_limit_infinite()
{
    std::cout << "Testing left limit infinite for tanh_sinh in 'A Comparison of Three High Precision Quadrature Schemes' on type " << boost::typeindex::type_id<Real>().pretty_name() << "\n";
    Real tol = sqrt(boost::math::tools::epsilon<Real>());
    Real Q;
    Real Q_expected;
    tanh_sinh<Real> integrator;

    // Example 11:
    auto f1 = [](Real t) { return 1/(1+t*t);};
    Q = integrator.integrate(f1, std::numeric_limits<Real>::has_infinity ? -std::numeric_limits<Real>::infinity() : -boost::math::tools::max_value<Real>(), 0);
    Q_expected = half_pi<Real>();
    BOOST_CHECK_CLOSE(Q, Q_expected, 100*tol);
}


// A horrible function taken from
// http://www.chebfun.org/examples/quad/GaussClenCurt.html
template<class Real>
void test_horrible()
{
    std::cout << "Testing Trefenthen's horrible integral on type " << boost::typeindex::type_id<Real>().pretty_name() << "\n";
    // We only know the integral to double precision, so requesting a higher tolerance doesn't make sense.
    Real tol = sqrt(std::numeric_limits<float>::epsilon());
    Real Q;
    Real Q_expected;
    Real error;
    Real L1;
    tanh_sinh<Real> integrator;

    auto f = [](Real x) { return x*sin(2*exp(2*sin(2*exp(2*x) ) ) ); };
    Q = integrator.integrate(f, (Real) -1, (Real) 1, &error, &L1);
    // NIntegrate[x*Sin[2*Exp[2*Sin[2*Exp[2*x]]]], {x, -1, 1}, WorkingPrecision -> 130, MaxRecursion -> 100]
    Q_expected = boost::lexical_cast<Real>("0.33673283478172753598559003181355241139806404130031017259552729882281");
    BOOST_CHECK_CLOSE_FRACTION(Q, Q_expected, tol);
}

// Some examples of tough integrals from NR, section 4.5.4:
template<class Real>
void test_nr_examples()
{
    std::cout << "Testing singular integrals from NR 4.5.4 on type " << boost::typeindex::type_id<Real>().pretty_name() << "\n";
    Real tol = sqrt(boost::math::tools::epsilon<Real>());
    Real Q;
    Real Q_expected;
    Real error;
    Real L1;
    tanh_sinh<Real> integrator(tol, 15);

    auto f1 = [](Real x) { return sin(x*half<Real>())*pow(x, -3*half<Real>())*exp(-x); };
    Q = integrator.integrate(f1, 0, std::numeric_limits<Real>::has_infinity ? std::numeric_limits<Real>::infinity() : boost::math::tools::max_value<Real>(), &error, &L1);
    Q_expected = sqrt(pi<Real>()*(sqrt((Real) 5) - 2));
    BOOST_CHECK_CLOSE_FRACTION(Q, Q_expected, 10*tol);

    auto f2 = [](Real x) { return pow(x, -(Real) 2/(Real) 7)*exp(-x*x); };
    Q = integrator.integrate(f2, 0, std::numeric_limits<Real>::has_infinity ? std::numeric_limits<Real>::infinity() : boost::math::tools::max_value<Real>());
    Q_expected = half<Real>()*boost::math::tgamma((Real) 5/ (Real) 14);
    BOOST_CHECK_CLOSE_FRACTION(Q, Q_expected, tol);

}

// Test integrand known to fool some termination schemes:
template<class Real>
void test_early_termination()
{
    std::cout << "Testing Clenshaw & Curtis's example of integrand which fools termination schemes on type " << boost::typeindex::type_id<Real>().pretty_name() << "\n";
    Real tol = sqrt(boost::math::tools::epsilon<Real>());
    Real Q;
    Real Q_expected;
    Real error;
    Real L1;
    tanh_sinh<Real> integrator;

    auto f1 = [](Real x) { return 23*cosh(x)/ (Real) 25 - cos(x) ; };
    Q = integrator.integrate(f1, (Real) -1, (Real) 1, &error, &L1);
    Q_expected = 46*sinh((Real) 1)/(Real) 25 - 2*sin((Real) 1);
    BOOST_CHECK_CLOSE_FRACTION(Q, Q_expected, tol);

}

// Test some definite integrals from the CRC handbook, 32nd edition:
template<class Real>
void test_crc()
{
    std::cout << "Testing CRC formulas on type " << boost::typeindex::type_id<Real>().pretty_name() << "\n";
    Real tol = sqrt(boost::math::tools::epsilon<Real>());
    Real Q;
    Real Q_expected;
    Real error;
    Real L1;
    tanh_sinh<Real> integrator;

    // CRC Definite integral 585
    auto f1 = [](Real x) { Real t = log(1/x); return x*x*t*t*t; };
    Q = integrator.integrate(f1, (Real) 0, (Real) 1, &error, &L1);
    Q_expected = (Real) 2/ (Real) 27;
    BOOST_CHECK_CLOSE_FRACTION(Q, Q_expected, tol);

    // CRC 636:
    //auto f2 = [](Real x) { return sqrt(cos(x)); };
    //Q = integrator.integrate(f2, (Real) 0, (Real) pi<Real>(), &error, &L1);
    //Q_expected = pow(two_pi<Real>(), 3*half<Real>())/pow(tgamma(0.25), 2);
    //BOOST_CHECK_CLOSE(Q, Q_expected, 100*tol);
}

template <class Real>
void test_sf()
{
   using std::sqrt;
   // Test some special functions that we already know how to evaluate:
   std::cout << "Testing special functions on type " << boost::typeindex::type_id<Real>().pretty_name() << "\n";
   Real tol = sqrt(boost::math::tools::epsilon<Real>());
   tanh_sinh<Real> integrator;
   // incomplete beta:
   if (std::numeric_limits<Real>::digits10 < 37) // Otherwise too slow
   {
      BOOST_CHECK_CLOSE_FRACTION(integrator.integrate([](Real x) { return boost::math::ibeta_derivative(100, 15, x); }, 0, Real(0.25)), boost::math::ibeta(100, 15, Real(0.25)), tol);
   }

   Real x = 2, y = 3, z = 0.5, p = 0.25;
   // This one has a higher than expected error rate, even though exp_sinh quadrature is just fine:
   //BOOST_CHECK_CLOSE_FRACTION(integrator.integrate([&](Real t) { return 1 / (sqrt(t + x) * (t + y)); }, 0, std::numeric_limits<Real>::has_infinity ? std::numeric_limits<Real>::infinity() : boost::math::tools::max_value<Real>()) / 2, boost::math::ellint_rc(x, y), tol);

   BOOST_CHECK_CLOSE_FRACTION(Real((Real(3) / 2) * integrator.integrate([&](Real t) { return 1 / (sqrt((t + x) * (t + y) * (t + z)) * (t + p)); }, 0, std::numeric_limits<Real>::has_infinity ? std::numeric_limits<Real>::infinity() : boost::math::tools::max_value<Real>())), boost::math::ellint_rj(x, y, z, p), tol);

   z = 5.5;
   BOOST_CHECK_CLOSE_FRACTION(integrator.integrate([&](Real t) { using std::pow; using std::exp; return t > 10000 ? Real(0) : Real(pow(t, z - 1) * exp(-t)); }, 0, std::numeric_limits<Real>::has_infinity ? std::numeric_limits<Real>::infinity() : boost::math::tools::max_value<Real>()),
      boost::math::tgamma(z), tol);
   BOOST_CHECK_CLOSE_FRACTION(integrator.integrate([](Real t) {  using std::exp; return exp(-t*t); }, std::numeric_limits<Real>::has_infinity ? -std::numeric_limits<Real>::infinity() : -boost::math::tools::max_value<Real>(), std::numeric_limits<Real>::has_infinity ? std::numeric_limits<Real>::infinity() : boost::math::tools::max_value<Real>()),
      boost::math::constants::root_pi<Real>(), tol);
}

BOOST_AUTO_TEST_CASE(tanh_sinh_quadrature_test)
{
    //generate_constants();

#ifdef TEST1

    test_detail<float>();
    test_right_limit_infinite<float>();
    test_left_limit_infinite<float>();
    test_linear<float>();
    test_quadratic<float>();
    test_singular<float>();
    test_ca<float>();
    test_three_quadrature_schemes_examples<float>();
    test_horrible<float>();
    test_integration_over_real_line<float>();
    test_nr_examples<float>();
    test_early_termination<float>();
    test_crc<float>();

#endif
#ifdef TEST2

    test_detail<double>();
    test_right_limit_infinite<double>();
    test_left_limit_infinite<double>();
    test_linear<double>();
    test_quadratic<double>();
    test_singular<double>();
    test_ca<double>();
    test_three_quadrature_schemes_examples<double>();
    test_horrible<double>();
    test_integration_over_real_line<double>();
    test_nr_examples<double>();
    test_early_termination<double>();
    test_crc<double>();
    test_sf<double>();

#endif

#ifdef TEST3

    test_detail<long double>();
    test_right_limit_infinite<long double>();
    test_left_limit_infinite<long double>();
    test_linear<long double>();
    test_quadratic<long double>();
    test_singular<long double>();
    test_ca<long double>();
    test_three_quadrature_schemes_examples<long double>();
    test_horrible<long double>();
    test_integration_over_real_line<long double>();
    test_nr_examples<long double>();
    test_early_termination<long double>();
    test_crc<long double>();
    test_sf<long double>();

#endif

#ifdef TEST4

    test_detail<cpp_bin_float_quad>();
    test_right_limit_infinite<cpp_bin_float_quad>();
    test_left_limit_infinite<cpp_bin_float_quad>();
    test_linear<cpp_bin_float_quad>();
    test_quadratic<cpp_bin_float_quad>();
    test_singular<cpp_bin_float_quad>();
    test_ca<cpp_bin_float_quad>();
    test_three_quadrature_schemes_examples<cpp_bin_float_quad>();
    test_horrible<cpp_bin_float_quad>();
    test_nr_examples<cpp_bin_float_quad>();
    test_early_termination<cpp_bin_float_quad>();
    test_crc<cpp_bin_float_quad>();
    test_sf<cpp_bin_float_quad>();

#endif
#ifdef TEST5

    test_sf<cpp_bin_float_50>();
    test_sf<cpp_bin_float_100>();
    test_sf<boost::multiprecision::number<boost::multiprecision::cpp_bin_float<150> > >();

#endif
#ifdef TEST6

    test_detail<boost::math::concepts::real_concept>();
    test_right_limit_infinite<boost::math::concepts::real_concept>();
    test_left_limit_infinite<boost::math::concepts::real_concept>();
    test_linear<boost::math::concepts::real_concept>();
    test_quadratic<boost::math::concepts::real_concept>();
    test_singular<boost::math::concepts::real_concept>();
    test_ca<boost::math::concepts::real_concept>();
    test_three_quadrature_schemes_examples<boost::math::concepts::real_concept>();
    test_horrible<boost::math::concepts::real_concept>();
    test_integration_over_real_line<boost::math::concepts::real_concept>();
    test_nr_examples<boost::math::concepts::real_concept>();
    test_early_termination<boost::math::concepts::real_concept>();
    test_crc<boost::math::concepts::real_concept>();
    test_sf<boost::math::concepts::real_concept>();

#endif
#ifdef TEST7
    test_sf<cpp_dec_float_50>();
#endif
}
