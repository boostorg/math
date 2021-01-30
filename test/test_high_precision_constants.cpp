//  Copyright Matt Borland 2021.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <limits>
#include <cmath>
#include <boost/math/constants/constants.hpp>
#include <boost/math/constants/high_precision_constants.hpp>
#include <boost/multiprecision/cpp_bin_float.hpp>
#include <boost/core/lightweight_test.hpp>

using boost::multiprecision::cpp_bin_float_50;
using boost::multiprecision::cpp_bin_float_100;
using std::abs;

template<typename T>
void test_pi()
{
    const T tol = std::numeric_limits<T>::epsilon();

    T old_pi = boost::math::constants::pi<T>();
    T new_pi = boost::math::numbers::pi_v<T>();

    BOOST_TEST(abs(old_pi-new_pi) < tol);
}

template<typename T>
void test_e()
{
    const T tol = std::numeric_limits<T>::epsilon();

    T old_e = boost::math::constants::e<T>();
    T new_e = boost::math::numbers::e_v<T>();

    BOOST_TEST(abs(old_e-new_e) < tol);
}

int main(void)
{
    test_pi<float>();
    test_pi<double>();
    test_pi<long double>();
    test_pi<cpp_bin_float_50>();
    test_pi<cpp_bin_float_100>();

    test_e<float>();
    test_e<double>();
    test_e<long double>();
    test_e<cpp_bin_float_50>();
    test_e<cpp_bin_float_100>();

    boost::report_errors();
}
