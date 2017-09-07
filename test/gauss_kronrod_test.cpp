// Copyright Nick Thompson 2017.
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)


#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE GaussKronrod

#include <boost/test/unit_test.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/multiprecision/cpp_bin_float.hpp>
#include <boost/math/quadrature/detail/gauss_kronrod_detail.hpp>

using boost::math::constants::half_pi;
using boost::math::quadrature::detail::g7_k15;
using boost::math::quadrature::detail::g10_k21;
using boost::multiprecision::cpp_bin_float_quad;


template<class Real>
void test_g7_k15()
{
    std::cout << std::setprecision(std::numeric_limits<Real>::digits10) << std::fixed;
    auto f = [](Real x) { return 1/(1+x*x); };
    Real error;
    Real Q = g7_k15<Real, decltype(f)>(f, &error);
    std::cout << Q << std::endl;
    std::cout << "Q - pi/2 = " << Q - half_pi<Real>() << std::endl;
    std::cout << "Error estimate: " << error << std::endl;
}

template<class Real>
void test_g10_k21()
{
    std::cout << std::setprecision(std::numeric_limits<Real>::digits10) << std::fixed;
    auto f = [](Real x) { return 1/(1+x*x); };
    Real error;
    Real Q = g10_k21<Real, decltype(f)>(f, &error);
    std::cout << Q << std::endl;
    std::cout << "Q - pi/2 = " << Q - half_pi<Real>() << std::endl;
    std::cout << "Error estimate: " << error << std::endl;
}



BOOST_AUTO_TEST_CASE(GaussKronrod)
{
    test_g7_k15<double>();
    test_g7_k15<long double>();
    test_g7_k15<cpp_bin_float_quad>();
    test_g7_k15<boost::multiprecision::cpp_bin_float_100>();

    test_g10_k21<double>();
    test_g10_k21<long double>();
    test_g10_k21<cpp_bin_float_quad>();
    test_g10_k21<boost::multiprecision::cpp_bin_float_100>();

}
