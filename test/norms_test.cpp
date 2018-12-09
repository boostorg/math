/*
 *  (C) Copyright Nick Thompson 2018.
 *  Use, modification and distribution are subject to the
 *  Boost Software License, Version 1.0. (See accompanying file
 *  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
 */

#include <vector>
#include <array>
#include <forward_list>
#include <algorithm>
#include <random>
#include <boost/core/lightweight_test.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/math/tools/norms.hpp>
#include <boost/multiprecision/cpp_bin_float.hpp>
#include <boost/multiprecision/cpp_complex.hpp>

using boost::multiprecision::cpp_bin_float_50;
using boost::multiprecision::cpp_complex_50;

/*
 * Test checklist:
 * 1) Does it work with multiprecision?
 * 2) Does it work with .cbegin()/.cend() if the data is not altered?
 * 3) Does it work with ublas and std::array? (Checking Eigen and Armadillo will make the CI system really unhappy.)
 * 4) Does it work with std::forward_list if a forward iterator is all that is required?
 * 5) Does it work with complex data if complex data is sensible?
 */


template<class Real>
void test_lp()
{
    Real tol = 50*std::numeric_limits<Real>::epsilon();

    std::array<Real, 3> u{1,0,0};
    Real l3 = boost::math::tools::lp_norm(u.begin(), u.end(), 3);
    BOOST_TEST(abs(l3 - 1) < tol);

    u[0] = -8;
    l3 = boost::math::tools::lp_norm(u.cbegin(), u.cend(), 3);
    BOOST_TEST(abs(l3 - 8) < tol);

    std::vector<Real> v(500);
    for (size_t i = 0; i < v.size(); ++i) {
        v[i] = 7;
    }
    Real l8 = boost::math::tools::lp_norm(v, 8);
    Real expected = 7*pow(v.size(), static_cast<Real>(1)/static_cast<Real>(8));
    BOOST_TEST(abs(l8 - expected) < tol*abs(expected));

    // Does it work with ublas vectors?
    // Does it handle the overflow of intermediates?
    boost::numeric::ublas::vector<Real> w(4);
    Real bignum = sqrt(std::numeric_limits<Real>::max())/256;
    for (size_t i = 0; i < w.size(); ++i)
    {
        w[i] = bignum;
    }
    Real l20 = boost::math::tools::lp_norm(w.cbegin(), w.cend(), 4);
    expected = bignum*pow(w.size(), static_cast<Real>(1)/static_cast<Real>(4));
    BOOST_TEST(abs(l20 - expected) < tol*expected);
}


template<class Complex>
void test_complex_lp()
{
    typedef typename Complex::value_type Real;
    Real tol = std::numeric_limits<Real>::epsilon();
    std::vector<Complex> v{{1,0}, {0,0}, {0,0}};
    Real l3 = boost::math::tools::lp_norm(v.cbegin(), v.cend(), 3);
    BOOST_TEST(abs(l3 - 1) < tol);

    l3 = boost::math::tools::lp_norm(v, 3);
    BOOST_TEST(abs(l3 - 1) < tol);

}

template<class Z>
void test_integer_total_variation()
{
    std::vector<Z> v{1,1};
    Z tv = boost::math::tools::total_variation(v);
    BOOST_TEST_EQ(tv,0);

    v[1] = 2;
    tv = boost::math::tools::total_variation(v.begin(), v.end());
    BOOST_TEST_EQ(tv,1);

    v.resize(50);
    for (size_t i = 0; i < v.size(); ++i) {
        v[i] = i;
    }

    tv = boost::math::tools::total_variation(v);
    BOOST_TEST_EQ(tv, v.size() -1);

    for (size_t i = 0; i < v.size(); ++i) {
        v[i] = i*i;
    }

    tv = boost::math::tools::total_variation(v);
    BOOST_TEST_EQ(tv, (v.size() -1)*(v.size()-1));

    // Work with std::array?
    std::array<Z, 2> w{1,1};
    tv = boost::math::tools::total_variation(w);
    BOOST_TEST_EQ(tv,0);
}

template<class Real>
void test_total_variation()
{
    Real tol = std::numeric_limits<Real>::epsilon();
    std::vector<Real> v{1,1};
    Real tv = boost::math::tools::total_variation(v.begin(), v.end());
    BOOST_TEST(tv >= 0 && abs(tv) < tol);

    tv = boost::math::tools::total_variation(v);
    BOOST_TEST(tv >= 0 && abs(tv) < tol);

    v[1] = 2;
    tv = boost::math::tools::total_variation(v.begin(), v.end());
    BOOST_TEST(abs(tv - 1) < tol);

    v.resize(50);
    for (size_t i = 0; i < v.size(); ++i) {
        v[i] = i;
    }

    tv = boost::math::tools::total_variation(v.begin(), v.end());
    BOOST_TEST(abs(tv - (v.size() -1)) < tol);

    for (size_t i = 0; i < v.size(); ++i) {
        v[i] = i*i;
    }

    tv = boost::math::tools::total_variation(v.begin(), v.end());
    BOOST_TEST(abs(tv - (v.size() -1)*(v.size()-1)) < tol);
}

template<class Real>
void test_sup_norm()
{
    Real tol = std::numeric_limits<Real>::epsilon();
    std::vector<Real> v{-2,1,0};
    Real s = boost::math::tools::sup_norm(v.begin(), v.end());
    BOOST_TEST(abs(s - 2) < tol);

    s = boost::math::tools::sup_norm(v);
    BOOST_TEST(abs(s - 2) < tol);

    // Work with std::array?
    std::array<Real, 3> w{-2,1,0};
    s = boost::math::tools::sup_norm(w);
    BOOST_TEST(abs(s - 2) < tol);

}

template<class Z>
void test_integer_sup_norm()
{
    std::vector<Z> v{-2,1,0};
    Z s = boost::math::tools::sup_norm(v.begin(), v.end());
    BOOST_TEST_EQ(s, 2);

    s = boost::math::tools::sup_norm(v);
    BOOST_TEST_EQ(s,2);
}

template<class Complex>
void test_complex_sup_norm()
{
    typedef typename Complex::value_type Real;
    Real tol = std::numeric_limits<Real>::epsilon();
    std::vector<Complex> w{{0,-8}, {1,1}, {3,2}};
    Real s = boost::math::tools::sup_norm(w.cbegin(), w.cend());
    BOOST_TEST(abs(s-8) < tol);

    s = boost::math::tools::sup_norm(w);
    BOOST_TEST(abs(s-8) < tol);
}

template<class Real>
void test_l0_pseudo_norm()
{
    std::vector<Real> v{0,0,1};
    size_t count = boost::math::tools::l0_pseudo_norm(v.begin(), v.end());
    BOOST_TEST_EQ(count, 1);

    // Compiles with cbegin()/cend()?
    count = boost::math::tools::l0_pseudo_norm(v.cbegin(), v.cend());
    BOOST_TEST_EQ(count, 1);

    count = boost::math::tools::l0_pseudo_norm(v);
    BOOST_TEST_EQ(count, 1);

    std::array<Real, 3> w{0,0,1};
    count = boost::math::tools::l0_pseudo_norm(w);
    BOOST_TEST_EQ(count, 1);
}

template<class Complex>
void test_complex_l0_pseudo_norm()
{
    std::vector<Complex> v{{0,0}, {0,0}, {1,0}};
    size_t count = boost::math::tools::l0_pseudo_norm(v.begin(), v.end());
    BOOST_TEST_EQ(count, 1);

    count = boost::math::tools::l0_pseudo_norm(v);
    BOOST_TEST_EQ(count, 1);

}

template<class Real>
void test_l1_norm()
{
    Real tol = std::numeric_limits<Real>::epsilon();
    std::vector<Real> v{1,1,1};
    Real l1 = boost::math::tools::l1_norm(v.begin(), v.end());
    BOOST_TEST(abs(l1 - 3) < tol);

    l1 = boost::math::tools::l1_norm(v);
    BOOST_TEST(abs(l1 - 3) < tol);

    std::array<Real, 3> w{1,1,1};
    l1 = boost::math::tools::l1_norm(w);
    BOOST_TEST(abs(l1 - 3) < tol);
}

template<class Z>
void test_integer_l1_norm()
{
    std::vector<Z> v{1,1,1};
    Z l1 = boost::math::tools::l1_norm(v.begin(), v.end());
    BOOST_TEST_EQ(l1, 3);
}

template<class Complex>
void test_complex_l1_norm()
{
    typedef typename Complex::value_type Real;
    Real tol = std::numeric_limits<Real>::epsilon();
    std::vector<Complex> v{{1,0}, {0,1},{0,-1}};
    Real l1 = boost::math::tools::l1_norm(v.begin(), v.end());
    BOOST_TEST(abs(l1 - 3) < tol);

    l1 = boost::math::tools::l1_norm(v);
    BOOST_TEST(abs(l1 - 3) < tol);

}

template<class Real>
void test_l2_norm()
{
    using std::sqrt;
    Real tol = std::numeric_limits<Real>::epsilon();
    std::vector<Real> v{1,1,1,1};
    Real l2 = boost::math::tools::l2_norm(v.begin(), v.end());
    BOOST_TEST(abs(l2 - 2) < tol);

    l2 = boost::math::tools::l2_norm(v);
    BOOST_TEST(abs(l2 - 2) < tol);

    std::array<Real, 4> w{1,1,1,1};
    l2 = boost::math::tools::l2_norm(w);
    BOOST_TEST(abs(l2 - 2) < tol);

    Real bignum = 4*sqrt(std::numeric_limits<Real>::max());
    v[0] = bignum;
    v[1] = 0;
    v[2] = 0;
    v[3] = 0;
    l2 = boost::math::tools::l2_norm(v.begin(), v.end());
    BOOST_TEST(abs(l2 - bignum) < tol*l2);
}

template<class Complex>
void test_complex_l2_norm()
{
    using std::sqrt;
    typedef typename Complex::value_type Real;
    Real tol = 100*std::numeric_limits<Real>::epsilon();
    std::vector<Complex> v{{1,0}, {0,1},{0,-1}, {1,0}};
    Real l2 = boost::math::tools::l2_norm(v.begin(), v.end());
    BOOST_TEST(abs(l2 - 2) < tol);

    l2 = boost::math::tools::l2_norm(v);
    BOOST_TEST(abs(l2 - 2) < tol);

}

int main()
{
    test_lp<float>();
    test_lp<double>();
    test_lp<long double>();
    test_lp<cpp_bin_float_50>();

    test_complex_lp<std::complex<float>>();
    test_complex_lp<std::complex<double>>();
    test_complex_lp<std::complex<long double>>();
    test_complex_lp<cpp_complex_50>();

    test_sup_norm<float>();
    test_sup_norm<double>();
    test_sup_norm<long double>();
    test_sup_norm<cpp_bin_float_50>();

    test_integer_sup_norm<int>();

    test_complex_sup_norm<std::complex<float>>();
    test_complex_sup_norm<std::complex<double>>();
    test_complex_sup_norm<std::complex<long double>>();
    test_complex_sup_norm<cpp_complex_50>();

    test_l0_pseudo_norm<int>();
    test_l0_pseudo_norm<float>();
    test_l0_pseudo_norm<double>();
    test_l0_pseudo_norm<long double>();
    test_l0_pseudo_norm<cpp_bin_float_50>();

    test_complex_l0_pseudo_norm<std::complex<float>>();
    test_complex_l0_pseudo_norm<std::complex<double>>();
    test_complex_l0_pseudo_norm<std::complex<long double>>();
    test_complex_l0_pseudo_norm<cpp_complex_50>();

    test_l1_norm<float>();
    test_l1_norm<double>();
    test_l1_norm<long double>();
    test_l1_norm<cpp_bin_float_50>();

    test_integer_l1_norm<int>();

    test_complex_l2_norm<std::complex<float>>();
    test_complex_l2_norm<std::complex<double>>();
    test_complex_l2_norm<std::complex<long double>>();
    test_complex_l2_norm<cpp_complex_50>();

    test_l2_norm<float>();
    test_l2_norm<double>();
    test_l2_norm<long double>();
    test_l2_norm<cpp_bin_float_50>();

    test_complex_l1_norm<std::complex<float>>();
    test_complex_l1_norm<std::complex<double>>();
    test_complex_l1_norm<std::complex<long double>>();
    test_complex_l1_norm<cpp_complex_50>();

    test_total_variation<float>();
    test_total_variation<double>();
    test_total_variation<long double>();
    test_total_variation<cpp_bin_float_50>();

    test_integer_total_variation<int>();

    return boost::report_errors();
}
