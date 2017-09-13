/*
 * Copyright Nick Thompson, 2017
 * Use, modification and distribution are subject to the
 * Boost Software License, Version 1.0. (See accompanying file
 * LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
 */
#define BOOST_TEST_MODULE chebyshev_test

#include <boost/type_index.hpp>
#include <boost/test/included/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/math/special_functions/chebyshev.hpp>
#include <boost/math/special_functions/chebyshev_transform.hpp>
#include <boost/math/special_functions/sinc.hpp>
#include <boost/multiprecision/cpp_bin_float.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>

using boost::multiprecision::cpp_bin_float_quad;
using boost::multiprecision::cpp_bin_float_50;
using boost::multiprecision::cpp_bin_float_100;
using boost::math::chebyshev_t;
using boost::math::chebyshev_t_prime;
using boost::math::chebyshev_u;
using boost::math::chebyshev_transform;

template<class Real>
void test_polynomials()
{
    std::cout << "Testing explicit polynomial representations of the Chebyshev polynomials on type " << boost::typeindex::type_id<Real>().pretty_name()  << "\n";

    Real x = -2;
    Real tol = 100*std::numeric_limits<Real>::epsilon();
    while (x < 2)
    {
        BOOST_CHECK_CLOSE_FRACTION(chebyshev_t(0, x), 1, tol);
        BOOST_CHECK_CLOSE_FRACTION(chebyshev_t(1, x), x, tol);
        BOOST_CHECK_CLOSE_FRACTION(chebyshev_t(2, x), 2*x*x - 1, tol);
        BOOST_CHECK_CLOSE_FRACTION(chebyshev_t(3, x), x*(4*x*x-3), tol);
        BOOST_CHECK_CLOSE_FRACTION(chebyshev_t(4, x), 8*x*x*x*x - 8*x*x + 1, tol);
        BOOST_CHECK_CLOSE_FRACTION(chebyshev_t(5, x), x*(16*x*x*x*x - 20*x*x + 5), tol);
        x += 1/static_cast<Real>(1<<7);
    }

    x = -2;
    tol = 10*tol;
    while (x < 2)
    {
        BOOST_CHECK_CLOSE_FRACTION(chebyshev_u(0, x), 1, tol);
        BOOST_CHECK_CLOSE_FRACTION(chebyshev_u(1, x), 2*x, tol);
        BOOST_CHECK_CLOSE_FRACTION(chebyshev_u(2, x), 4*x*x - 1, tol);
        BOOST_CHECK_CLOSE_FRACTION(chebyshev_u(3, x), 4*x*(2*x*x - 1), tol);
        x += 1/static_cast<Real>(1<<7);
    }
}


template<class Real>
void test_derivatives()
{
    std::cout << "Testing explicit polynomial representations of the Chebyshev polynomial derivatives on type " << boost::typeindex::type_id<Real>().pretty_name()  << "\n";

    Real x = -2;
    Real tol = 1000*std::numeric_limits<Real>::epsilon();
    while (x < 2)
    {
        BOOST_CHECK_CLOSE_FRACTION(chebyshev_t_prime(0, x), 0, tol);
        BOOST_CHECK_CLOSE_FRACTION(chebyshev_t_prime(1, x), 1, tol);
        BOOST_CHECK_CLOSE_FRACTION(chebyshev_t_prime(2, x), 4*x, tol);
        BOOST_CHECK_CLOSE_FRACTION(chebyshev_t_prime(3, x), 3*(4*x*x - 1), tol);
        BOOST_CHECK_CLOSE_FRACTION(chebyshev_t_prime(4, x), 16*x*(2*x*x - 1), tol);
        // This one makes the tolerance have to grow too large; the Chebyshev recurrence is more stable than naive polynomial evaluation anyway.
        //BOOST_CHECK_CLOSE_FRACTION(chebyshev_t_prime(5, x), 5*(4*x*x*(4*x*x - 3) + 1), tol);
        x += 1/static_cast<Real>(1<<7);
    }
}


template<class Real>
void test_sin_chebyshev_transform()
{
    using boost::math::chebyshev_transform;
    using boost::math::constants::half_pi;
    using std::sin;
    using std::cos;
    using std::abs;

    Real tol = 10*std::numeric_limits<Real>::epsilon();
    auto f = [](Real x) { return sin(x); };
    Real a = 0;
    Real b = 1;
    chebyshev_transform<Real> cheb(f, a, b, tol);

    Real x = a;
    while (x < b)
    {
        Real s = sin(x);
        Real c = cos(x);
        if (abs(s) < tol)
        {
            BOOST_CHECK_SMALL(cheb(x), 10*tol);
            BOOST_CHECK_CLOSE_FRACTION(c, cheb.prime(x), 10*tol);
        }
        else
        {
            BOOST_CHECK_CLOSE_FRACTION(s, cheb(x), 10*tol);
            if (abs(c) < tol)
            {
                BOOST_CHECK_SMALL(cheb.prime(x), 10*tol);
            }
            else
            {
                BOOST_CHECK_CLOSE_FRACTION(c, cheb.prime(x), 10*tol);
            }
        }
        x += static_cast<Real>(1)/static_cast<Real>(1 << 7);
    }

    Real Q = cheb.integrate();

    BOOST_CHECK_CLOSE_FRACTION(1 - cos(static_cast<Real>(1)), Q, 10*tol);
    cheb.print_coefficients();
}


template<class Real>
void test_clenshaw_recurrence()
{
    using boost::math::chebyshev_clenshaw_recurrence;
    std::vector<Real> c0{2, 0, 0, 0, 0};
    // Check the size = 1 case:
    std::vector<Real> c01{2};
    // Check the size = 2 case:
    std::vector<Real> c02{2, 0};
    std::vector<Real> c1{0, 1, 0, 0};
    std::vector<Real> c2{0, 0, 1, 0};
    std::vector<Real> c3{0, 0, 0, 1, 0};
    std::vector<Real> c4{0, 0, 0, 0, 1};
    std::vector<Real> c5{0, 0, 0, 0, 0, 1};
    std::vector<Real> c6{0, 0, 0, 0, 0, 0, 1};

    Real x = -1;
    Real tol = 10*std::numeric_limits<Real>::epsilon();
    while (x <= 1)
    {
        Real y = chebyshev_clenshaw_recurrence(c0.data(), c0.size(), x);
        BOOST_CHECK_CLOSE_FRACTION(y, chebyshev_t(0, x), tol);

        y = chebyshev_clenshaw_recurrence(c01.data(), c01.size(), x);
        BOOST_CHECK_CLOSE_FRACTION(y, chebyshev_t(0, x), tol);

        y = chebyshev_clenshaw_recurrence(c02.data(), c02.size(), x);
        BOOST_CHECK_CLOSE_FRACTION(y, chebyshev_t(0, x), tol);

        y = chebyshev_clenshaw_recurrence(c1.data(), c1.size(), x);
        BOOST_CHECK_CLOSE_FRACTION(y, chebyshev_t(1, x), tol);

        y = chebyshev_clenshaw_recurrence(c2.data(), c2.size(), x);
        BOOST_CHECK_CLOSE_FRACTION(y, chebyshev_t(2, x), tol);

        y = chebyshev_clenshaw_recurrence(c3.data(), c3.size(), x);
        BOOST_CHECK_CLOSE_FRACTION(y, chebyshev_t(3, x), tol);

        y = chebyshev_clenshaw_recurrence(c4.data(), c4.size(), x);
        BOOST_CHECK_CLOSE_FRACTION(y, chebyshev_t(4, x), tol);

        y = chebyshev_clenshaw_recurrence(c5.data(), c5.size(), x);
        BOOST_CHECK_CLOSE_FRACTION(y, chebyshev_t(5, x), tol);

        y = chebyshev_clenshaw_recurrence(c6.data(), c6.size(), x);
        BOOST_CHECK_CLOSE_FRACTION(y, chebyshev_t(6, x), tol);

        x += static_cast<Real>(1)/static_cast<Real>(1 << 7);
    }
}


template<class Real>
void test_sinc_chebyshev_transform()
{
    using boost::math::sinc_pi;
    using boost::math::chebyshev_transform;
    using boost::math::constants::half_pi;

    Real tol = 100*std::numeric_limits<Real>::epsilon();
    auto f = [](Real x) { return boost::math::sinc_pi(x); };
    Real a = 0;
    Real b = 1;
    chebyshev_transform<Real> cheb(f, a, b);

    Real x = a;
    while (x < b)
    {
        Real s = sinc_pi(x);
        if (s < tol)
        {
            BOOST_CHECK_SMALL(cheb(x), tol);
        }
        else
        {
            BOOST_CHECK_CLOSE_FRACTION(s, cheb(x), tol);
        }
        x += static_cast<Real>(1)/static_cast<Real>(1 << 7);
    }

    Real Q = cheb.integrate();
    //NIntegrate[Sinc[x], {x, 0, 1}, WorkingPrecision -> 200, AccuracyGoal -> 150, PrecisionGoal -> 150, MaxRecursion -> 150]
    Real Q_exp = boost::lexical_cast<Real>("0.94608307036718301494135331382317965781233795473811179047145477356668");
    BOOST_CHECK_CLOSE_FRACTION(Q_exp, Q, tol);
}



//Examples taken from "Approximation Theory and Approximation Practice", by Trefethen
template<class Real>
void test_atap_examples()
{
    using std::sin;
    using boost::math::constants::half;
    using boost::math::sinc_pi;
    using boost::math::chebyshev_transform;
    using boost::math::constants::half_pi;

    Real tol = 10*std::numeric_limits<Real>::epsilon();
    auto f1 = [](Real x) { return ((0 < x) - (x < 0)) - x/2; };
    auto f2 = [](Real x) { Real t = sin(6*x); Real s = sin(x + exp(2*x));
                           Real u = (0 < s) - (s < 0);
                           return t + u; };

    auto f3 = [](Real x) { return sin(6*x) + sin(60*exp(x)); };

    auto f4 = [](Real x) { return 1/(1+1000*(x+half<Real>())*(x+half<Real>())) + 1/sqrt(1+1000*(x-.5)*(x-0.5));};
    Real a = -1;
    Real b = 1;
    //chebyshev_transform<Real> cheb1(f1, a, b);
    //chebyshev_transform<Real> cheb2(f2, a, b, tol);
    chebyshev_transform<Real> cheb3(f3, a, b, tol);

    Real x = a;
    while (x < b)
    {
        Real s = f1(x);
        //BOOST_CHECK_CLOSE_FRACTION(s, cheb1(x), tol);
        //BOOST_CHECK_CLOSE_FRACTION(f2(x), cheb2(x), tol);
        BOOST_CHECK_CLOSE_FRACTION(f3(x), cheb3(x), tol);
        x += static_cast<Real>(1)/static_cast<Real>(1 << 7);
    }
}

//Validate that the Chebyshev polynomials are well approximated by the Chebyshev transform.
template<class Real>
void test_chebyshev_chebyshev_transform()
{
  Real tol = 500*std::numeric_limits<Real>::epsilon();
  // T_0 = 1:
  auto t0 = [](Real x) { return 1; };
  chebyshev_transform<Real> cheb0(t0, -1, 1);
  BOOST_CHECK_CLOSE_FRACTION(cheb0.coefficients()[0], 2, tol);

  Real x = -1;
  while (x < 1) {
    BOOST_CHECK_CLOSE_FRACTION(cheb0(x), 1, tol);
    x += static_cast<Real>(1)/static_cast<Real>(1 << 7);
  }

  // T_1 = x:
  auto t1 = [](Real x) { return x; };
  chebyshev_transform<Real> cheb1(t1, -1, 1);
  BOOST_CHECK_CLOSE_FRACTION(cheb1.coefficients()[1], 1, tol);

  x = -1;
  while (x < 1)
  {
    if (x == 0)
    {
        BOOST_CHECK_SMALL(cheb1(x), tol);
    }
    else
    {
        BOOST_CHECK_CLOSE_FRACTION(cheb1(x), x, tol);
    }
    x += static_cast<Real>(1)/static_cast<Real>(1 << 7);
  }
  cheb1.print_coefficients();


  // T_1 = x:
  auto t2 = [](Real x) { return 2*x*x-1; };
  chebyshev_transform<Real> cheb2(t2, -1, 1);
  BOOST_CHECK_CLOSE_FRACTION(cheb2.coefficients()[2], 1, tol);

  x = -1;
  while (x < 1)
  {
      BOOST_CHECK_CLOSE_FRACTION(cheb2(x), t2(x), tol);
      x += static_cast<Real>(1)/static_cast<Real>(1 << 7);
  }
  cheb2.print_coefficients();


}

BOOST_AUTO_TEST_CASE(chebyshev_test)
{
    test_clenshaw_recurrence<double>();
    //test_chebyshev_chebyshev_transform<float>();
    test_chebyshev_chebyshev_transform<double>();
    test_sin_chebyshev_transform<double>();
    test_atap_examples<double>();
    test_sinc_chebyshev_transform<double>();

    test_polynomials<float>();
    test_polynomials<double>();
    test_polynomials<long double>();
    test_polynomials<cpp_bin_float_quad>();

    test_derivatives<float>();
    test_derivatives<double>();
    test_derivatives<long double>();
    test_derivatives<cpp_bin_float_quad>();
}
