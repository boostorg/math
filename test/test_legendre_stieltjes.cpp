// Copyright Nick Thompson 2017.
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <iostream>
#include <iomanip>
#include <vector>
#include <limits>
#include "math_unit_test.hpp"
#include <boost/math/special_functions/legendre.hpp>
#include <boost/math/special_functions/legendre_stieltjes.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/multiprecision/cpp_bin_float.hpp>


using boost::math::legendre_stieltjes;
using boost::math::legendre_p;
using boost::multiprecision::cpp_bin_float_quad;


template<class Real>
void test_legendre_stieltjes()
{
    std::cout << std::setprecision(std::numeric_limits<Real>::digits10);
    using std::sqrt;
    using std::abs;
    using boost::math::constants::third;
    using boost::math::constants::half;

    Real tol = std::numeric_limits<Real>::epsilon();
    legendre_stieltjes<Real> ls1(1);
    legendre_stieltjes<Real> ls2(2);
    legendre_stieltjes<Real> ls3(3);
    legendre_stieltjes<Real> ls4(4);
    legendre_stieltjes<Real> ls5(5);
    legendre_stieltjes<Real> ls8(8);
    Real x = -1;
    while(x <= 1)
    {
        CHECK_ULP_CLOSE(x, ls1(x), 2);
        CHECK_ULP_CLOSE(1, ls1.prime(x), 2);

        Real p2 = legendre_p(2, x);
        CHECK_ULP_CLOSE(p2 - 2/static_cast<Real>(5), ls2(x), 2);
        CHECK_ULP_CLOSE(3*x, ls2.prime(x), 2);

        Real p3 = legendre_p(3, x);
        CHECK_ULP_CLOSE(p3 - 9*x/static_cast<Real>(14), ls3(x), 256);
        CHECK_ULP_CLOSE(15*x*x*half<Real>() -3*half<Real>()-9/static_cast<Real>(14), ls3.prime(x), 4);

        Real p4 = legendre_p(4, x);
        //-20P_2(x)/27 + 14P_0(x)/891
        Real E4 = p4 - 20*p2/static_cast<Real>(27) + 14/static_cast<Real>(891);
        CHECK_ULP_CLOSE(E4, ls4(x), 500);
        CHECK_ULP_CLOSE(35*x*(9*x*x -5)/static_cast<Real>(18), ls4.prime(x), 400);

        Real p5 = legendre_p(5, x);
        Real E5 = p5 - 35*p3/static_cast<Real>(44) + 135*x/static_cast<Real>(12584);
        CHECK_ULP_CLOSE(E5, ls5(x), 58000);
        Real E5prime = (315*(123 + 143*x*x*(11*x*x-9)))/static_cast<Real>(12584);
        CHECK_ULP_CLOSE(E5prime, ls5.prime(x), 316);
        x += 1/static_cast<Real>(1 << 9);
    }

    // Test norm:
    // E_1 = x
    Real expected_norm_sq = 2*third<Real>();
    CHECK_ULP_CLOSE(expected_norm_sq, ls1.norm_sq(), 2);

    // E_2 = P[sub 2](x) - 2P[sup 0](x)/5
    expected_norm_sq = 2/static_cast<Real>(5) + 8/static_cast<Real>(25);
    CHECK_ULP_CLOSE(expected_norm_sq, ls2.norm_sq(), 2);

    // E_3 = P[sub 3](x) - 9P[sub 1]/14
    expected_norm_sq = 2/static_cast<Real>(7) + 9*9*2*third<Real>()/static_cast<Real>(14*14);
    CHECK_ULP_CLOSE(expected_norm_sq, ls3.norm_sq(), 2);

    // E_4 = P[sub 4](x) -20P[sub 2](x)/27 + 14P[sub 0](x)/891
    expected_norm_sq = static_cast<Real>(2)/static_cast<Real>(9) + static_cast<Real>(20*20*2)/static_cast<Real>(27*27*5) + 14*14*2/static_cast<Real>(891*891);
    CHECK_ULP_CLOSE(expected_norm_sq, ls4.norm_sq(), 2);

    // E_5 = P[sub 5](x) - 35P[sub 3](x)/44 + 135P[sub 1](x)/12584
    expected_norm_sq = 2/static_cast<Real>(11) + (35*35/static_cast<Real>(44*44))*(2/static_cast<Real>(7)) + (135*135/static_cast<Real>(12584*12584))*2*third<Real>();
    CHECK_ULP_CLOSE(expected_norm_sq, ls5.norm_sq(), 2);

    // Only zero of E1 is 0:
    std::vector<Real> zeros = ls1.zeros();
    CHECK_EQUAL(zeros.size(), static_cast<std::size_t>(1));
    CHECK_ABSOLUTE_ERROR(0, zeros[0], tol);
    CHECK_ABSOLUTE_ERROR(0, ls1(zeros[0]), tol);

    zeros = ls2.zeros();
    CHECK_EQUAL(zeros.size(), static_cast<std::size_t>(1));
    CHECK_ULP_CLOSE(sqrt(3/static_cast<Real>(5)), zeros[0], 2);
    CHECK_ABSOLUTE_ERROR(0, ls2(zeros[0]), tol);

    zeros = ls3.zeros();
    CHECK_EQUAL(zeros.size(), static_cast<std::size_t>(2));
    CHECK_ABSOLUTE_ERROR(0, zeros[0], tol);
    CHECK_ULP_CLOSE(sqrt(6/static_cast<Real>(7)), zeros[1], 2);


    zeros = ls4.zeros();
    CHECK_EQUAL(zeros.size(), static_cast<std::size_t>(2));
    Real expected = sqrt( (55 - 2*sqrt(static_cast<Real>(330)))/static_cast<Real>(11) )/static_cast<Real>(3);
    CHECK_ULP_CLOSE(expected, zeros[0], 2);

    expected = sqrt( (55 + 2*sqrt(static_cast<Real>(330)))/static_cast<Real>(11) )/static_cast<Real>(3);
    CHECK_ULP_CLOSE(expected, zeros[1], 4);


    zeros = ls5.zeros();
    CHECK_EQUAL(zeros.size(), static_cast<std::size_t>(3));
    CHECK_ABSOLUTE_ERROR(0, zeros[0], tol);

    expected = sqrt( ( 195 - sqrt(static_cast<Real>(6045)) )/static_cast<Real>(286));
    CHECK_ULP_CLOSE(expected, zeros[1], 2);

    expected = sqrt( ( 195 + sqrt(static_cast<Real>(6045)) )/static_cast<Real>(286));
    CHECK_ULP_CLOSE(expected, zeros[2], 2);


    for (size_t i = 6; i < 50; ++i)
    {
        legendre_stieltjes<Real> En(i);
        zeros = En.zeros();
        for(auto const & zero : zeros)
        {
            CHECK_ABSOLUTE_ERROR(0, En(zero), 50*tol);
        }
    }

#ifndef BOOST_MATH_NO_EXCEPTIONS
    CHECK_THROW(legendre_stieltjes<Real>(0), std::domain_error);
#endif
}


int main()
{
    test_legendre_stieltjes<double>();
#ifndef BOOST_MATH_NO_LONG_DOUBLE_MATH_FUNCTIONS
    test_legendre_stieltjes<long double>();
#endif
    test_legendre_stieltjes<cpp_bin_float_quad>();
    //test_legendre_stieltjes<boost::multiprecision::cpp_bin_float_100>();
    return boost::math::test::report_errors();
}
