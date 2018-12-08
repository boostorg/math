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
#include <boost/math/tools/vector_functionals.hpp>
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
void test_mean()
{
    Real tol = std::numeric_limits<Real>::epsilon();
    std::vector<Real> v{1,2,3,4,5};
    Real mu = boost::math::tools::mean(v.begin(), v.end());
    BOOST_TEST(abs(mu - 3) < tol);

    // Can we successfully average only part of the vector?
    mu = boost::math::tools::mean(v.begin(), v.begin() + 3);
    BOOST_TEST(abs(mu - 2) < tol);

    // Does it work when we const qualify?
    mu = boost::math::tools::mean(v.cbegin(), v.cend());
    BOOST_TEST(abs(mu - 3) < tol);

    // Does it work for std::array?
    std::array<Real, 7> u{1,2,3,4,5,6,7};
    mu = boost::math::tools::mean(u.begin(), u.end());
    BOOST_TEST(abs(mu - 4) < tol);

    // Does it work for a forward iterator?
    std::forward_list<Real> l{1,2,3,4,5,6,7};
    mu = boost::math::tools::mean(l.begin(), l.end());
    BOOST_TEST(abs(mu - 4) < tol);

    // Does it work with ublas vectors?
    boost::numeric::ublas::vector<Real> w(7);
    for (size_t i = 0; i < w.size(); ++i)
    {
        w[i] = i+1;
    }
    mu = boost::math::tools::mean(w.cbegin(), w.cend());
    BOOST_TEST(abs(mu - 4) < tol);

}

template<class Complex>
void test_complex_mean()
{
    typedef typename Complex::value_type Real;
    Real tol = std::numeric_limits<Real>::epsilon();
    std::vector<Complex> v{{0,1},{0,2},{0,3},{0,4},{0,5}};
    auto mu = boost::math::tools::mean(v.begin(), v.end());
    BOOST_TEST(abs(mu.imag() - 3) < tol);
    BOOST_TEST(abs(mu.real()) < tol);
}

template<class Real>
void test_mean_and_population_variance()
{
    Real tol = std::numeric_limits<Real>::epsilon();
    std::vector<Real> v{1,1,1,1,1,1};
    auto [mu, sigma_sq] = boost::math::tools::mean_and_population_variance(v.begin(), v.end());
    BOOST_TEST(abs(mu - 1) < tol);
    BOOST_TEST(abs(sigma_sq) < tol);

    std::vector<Real> u{1};
    auto [mu1, sigma1_sq] = boost::math::tools::mean_and_population_variance(u.cbegin(), u.cend());
    BOOST_TEST(abs(mu1 - 1) < tol);
    BOOST_TEST(abs(sigma1_sq) < tol);

    std::array<Real, 8> w{0,1,0,1,0,1,0,1};
    auto [mu2, sigma2_sq] = boost::math::tools::mean_and_population_variance(w.begin(), w.end());
    BOOST_TEST(abs(mu2 - 1.0/2.0) < tol);
    BOOST_TEST(abs(sigma2_sq - 1.0/4.0) < tol);
}

template<class Real>
void test_median()
{
    std::mt19937 g(12);
    std::vector<Real> v{1,2,3,4,5,6,7};

    Real m = boost::math::tools::median(v.begin(), v.end());
    BOOST_TEST_EQ(m, 4);

    std::shuffle(v.begin(), v.end(), g);
    m = boost::math::tools::median(v.begin(), v.end());
    BOOST_TEST_EQ(m, 4);

    v = {1,2,3,3,4,5};
    m = boost::math::tools::median(v.begin(), v.end());
    BOOST_TEST_EQ(m, 3);
    std::shuffle(v.begin(), v.end(), g);
    m = boost::math::tools::median(v.begin(), v.end());
    BOOST_TEST_EQ(m, 3);

    v = {1};
    m = boost::math::tools::median(v.begin(), v.end());
    BOOST_TEST_EQ(m, 1);

    v = {1,1};
    m = boost::math::tools::median(v.begin(), v.end());
    BOOST_TEST_EQ(m, 1);

    v = {2,4};
    m = boost::math::tools::median(v.begin(), v.end());
    BOOST_TEST_EQ(m, 3);

    v = {1,1,1};
    m = boost::math::tools::median(v.begin(), v.end());
    BOOST_TEST_EQ(m, 1);

    v = {1,2,3};
    m = boost::math::tools::median(v.begin(), v.end());
    BOOST_TEST_EQ(m, 2);
    std::shuffle(v.begin(), v.end(), g);
    m = boost::math::tools::median(v.begin(), v.end());
    BOOST_TEST_EQ(m, 2);
}

template<class Real>
void test_absolute_median()
{
    std::mt19937 g(12);
    std::vector<Real> v{-1, 2, -3, 4, -5, 6, -7};

    Real m = boost::math::tools::absolute_median(v.begin(), v.end());
    BOOST_TEST_EQ(m, 4);

    std::shuffle(v.begin(), v.end(), g);
    m = boost::math::tools::absolute_median(v.begin(), v.end());
    BOOST_TEST_EQ(m, 4);

    v = {1, -2, -3, 3, -4, -5};
    m = boost::math::tools::absolute_median(v.begin(), v.end());
    BOOST_TEST_EQ(m, 3);
    std::shuffle(v.begin(), v.end(), g);
    m = boost::math::tools::absolute_median(v.begin(), v.end());
    BOOST_TEST_EQ(m, 3);

    v = {-1};
    m = boost::math::tools::absolute_median(v.begin(), v.end());
    BOOST_TEST_EQ(m, 1);

    v = {-1, 1};
    m = boost::math::tools::absolute_median(v.begin(), v.end());
    BOOST_TEST_EQ(m, 1);

    v = {2, -4};
    m = boost::math::tools::absolute_median(v.begin(), v.end());
    BOOST_TEST_EQ(m, 3);

    v = {1, -1, 1};
    m = boost::math::tools::absolute_median(v.begin(), v.end());
    BOOST_TEST_EQ(m, 1);

    v = {1, 2, -3};
    m = boost::math::tools::absolute_median(v.begin(), v.end());
    BOOST_TEST_EQ(m, 2);
    std::shuffle(v.begin(), v.end(), g);
    m = boost::math::tools::absolute_median(v.begin(), v.end());
    BOOST_TEST_EQ(m, 2);
}


template<class Complex>
void test_complex_absolute_median()
{
    typedef typename Complex::value_type Real;
    std::mt19937 g(18);
    std::vector<Complex> v{{0,1}, {0,-2},{0,3}, {0,-4}, {0,5}, {0,-6}, {0,7}};

    Real m = boost::math::tools::absolute_median(v.begin(), v.end());
    BOOST_TEST_EQ(m, 4);

    std::shuffle(v.begin(), v.end(), g);
    m = boost::math::tools::absolute_median(v.begin(), v.end());
    BOOST_TEST_EQ(m, 4);

    v = {{0,1}, {0,-2}, {0,-3}, {0,3}, {0,4}, {0,-5}};
    m = boost::math::tools::absolute_median(v.begin(), v.end());
    BOOST_TEST_EQ(m, 3);
    std::shuffle(v.begin(), v.end(), g);
    m = boost::math::tools::absolute_median(v.begin(), v.end());
    BOOST_TEST_EQ(m, 3);

    v = {{0, -1}};
    m = boost::math::tools::absolute_median(v.begin(), v.end());
    BOOST_TEST_EQ(m, 1);
}


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
    Real l8 = boost::math::tools::lp_norm(v.cbegin(), v.cend(), 8);
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

}


template<class Real>
void test_total_variation()
{
    Real tol = std::numeric_limits<Real>::epsilon();
    std::vector<Real> v{1,1};
    Real tv = boost::math::tools::total_variation(v.begin(), v.end());
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

}

template<class Complex>
void test_complex_sup_norm()
{
    typedef typename Complex::value_type Real;
    Real tol = std::numeric_limits<Real>::epsilon();
    std::vector<Complex> w{{0,-8}, {1,1}, {3,2}};
    Real s = boost::math::tools::sup_norm(w.cbegin(), w.cend());
    BOOST_TEST(abs(s-8) < tol);
}


template<class Real>
void test_gini_coefficient()
{
    Real tol = std::numeric_limits<Real>::epsilon();
    std::vector<Real> v{1,0,0};
    Real gini = boost::math::tools::gini_coefficient(v.begin(), v.end());
    BOOST_TEST(abs(gini - 1) < tol);

    v[0] = 1;
    v[1] = 1;
    v[2] = 1;
    gini = boost::math::tools::gini_coefficient(v.begin(), v.end());
    BOOST_TEST(abs(gini) < tol);

    v[0] = 0;
    v[1] = 0;
    v[2] = 0;
    gini = boost::math::tools::gini_coefficient(v.begin(), v.end());
    BOOST_TEST(abs(gini) < tol);
}

template<class Real>
void test_hoyer_sparsity()
{
    using std::sqrt;
    Real tol = 5*std::numeric_limits<Real>::epsilon();
    std::vector<Real> v{1,0,0};
    Real hs = boost::math::tools::hoyer_sparsity(v.begin(), v.end());
    BOOST_TEST(abs(hs - 1) < tol);

    // Does it work with constant iterators?
    hs = boost::math::tools::hoyer_sparsity(v.cbegin(), v.cend());
    BOOST_TEST(abs(hs - 1) < tol);

    v[0] = 1;
    v[1] = 1;
    v[2] = 1;
    hs = boost::math::tools::hoyer_sparsity(v.cbegin(), v.cend());
    BOOST_TEST(abs(hs) < tol);
}

template<class Complex>
void test_complex_hoyer_sparsity()
{
    typedef typename Complex::value_type Real;
    using std::sqrt;
    Real tol = 5*std::numeric_limits<Real>::epsilon();
    std::vector<Complex> v{{0,1}, {0, 0}, {0,0}};
    Real hs = boost::math::tools::hoyer_sparsity(v.begin(), v.end());
    BOOST_TEST(abs(hs - 1) < tol);

    // Does it work with constant iterators?
    hs = boost::math::tools::hoyer_sparsity(v.cbegin(), v.cend());
    BOOST_TEST(abs(hs - 1) < tol);

    // All are the same magnitude:
    v[0] = {0, 1};
    v[1] = {1, 0};
    v[2] = {0,-1};
    hs = boost::math::tools::hoyer_sparsity(v.cbegin(), v.cend());
    BOOST_TEST(abs(hs) < tol);
}


template<class Real>
void test_absolute_gini_coefficient()
{
    Real tol = std::numeric_limits<Real>::epsilon();
    std::vector<Real> v{-1,0,0};
    Real gini = boost::math::tools::absolute_gini_coefficient(v.begin(), v.end());
    BOOST_TEST(abs(gini - 1) < tol);

    v[0] = 1;
    v[1] = -1;
    v[2] = 1;
    gini = boost::math::tools::absolute_gini_coefficient(v.begin(), v.end());
    BOOST_TEST(abs(gini) < tol);

    std::vector<std::complex<Real>> w(128);
    std::complex<Real> i{0,1};
    for(size_t k = 0; k < w.size(); ++k)
    {
        w[k] = exp(i*static_cast<Real>(k)/static_cast<Real>(w.size()));
    }
    gini = boost::math::tools::absolute_gini_coefficient(w.begin(), w.end());
    BOOST_TEST(abs(gini) < tol);

    // The Gini index is invariant under "cloning": If w = v \oplus v, then G(w) = G(v).
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

}

template<class Complex>
void test_complex_l0_pseudo_norm()
{
    std::vector<Complex> v{{0,0}, {0,0}, {1,0}};
    size_t count = boost::math::tools::l0_pseudo_norm(v.begin(), v.end());
    BOOST_TEST_EQ(count, 1);
}

template<class Real>
void test_l1_norm()
{
    Real tol = std::numeric_limits<Real>::epsilon();
    std::vector<Real> v{1,1,1};
    Real l1 = boost::math::tools::l1_norm(v.begin(), v.end());
    BOOST_TEST(abs(l1 - 3) < tol);
}

template<class Complex>
void test_complex_l1_norm()
{
    typedef typename Complex::value_type Real;
    Real tol = std::numeric_limits<Real>::epsilon();
    std::vector<Complex> v{{1,0}, {0,1},{0,-1}};
    Real l1 = boost::math::tools::l1_norm(v.begin(), v.end());
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
}

template<class Real>
void test_shannon_entropy()
{
    Real tol = 100*std::numeric_limits<Real>::epsilon();
    using boost::math::constants::half;
    using boost::math::constants::ln_two;
    std::vector<Real> v(30, half<Real>());
    Real Hs = boost::math::tools::shannon_entropy(v.begin(), v.end());
    Real expected = v.size()*ln_two<Real>()/2;
    BOOST_TEST(abs(Hs - expected) < tol*expected);
}


int main()
{
    test_mean<float>();
    test_mean<double>();
    test_mean<long double>();
    test_mean<cpp_bin_float_50>();

    test_complex_mean<std::complex<float>>();
    test_complex_mean<cpp_complex_50>();

    test_mean_and_population_variance<float>();
    test_mean_and_population_variance<double>();
    test_mean_and_population_variance<long double>();
    test_mean_and_population_variance<cpp_bin_float_50>();

    test_median<float>();
    test_median<double>();
    test_median<long double>();
    test_median<cpp_bin_float_50>();

    test_absolute_median<float>();
    test_absolute_median<double>();
    test_absolute_median<long double>();
    test_absolute_median<cpp_bin_float_50>();

    test_complex_absolute_median<std::complex<float>>();
    test_complex_absolute_median<std::complex<double>>();
    test_complex_absolute_median<std::complex<long double>>();
    test_complex_absolute_median<cpp_complex_50>();

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

    test_complex_sup_norm<std::complex<float>>();
    test_complex_sup_norm<std::complex<double>>();
    test_complex_sup_norm<std::complex<long double>>();
    test_complex_sup_norm<cpp_complex_50>();

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

    test_gini_coefficient<float>();
    test_gini_coefficient<double>();
    test_gini_coefficient<long double>();
    test_gini_coefficient<cpp_bin_float_50>();

    test_absolute_gini_coefficient<float>();
    test_absolute_gini_coefficient<double>();
    test_absolute_gini_coefficient<long double>();
    test_absolute_gini_coefficient<cpp_bin_float_50>();

    test_hoyer_sparsity<float>();
    test_hoyer_sparsity<double>();
    test_hoyer_sparsity<long double>();
    test_hoyer_sparsity<cpp_bin_float_50>();

    test_shannon_entropy<float>();
    test_shannon_entropy<double>();
    test_shannon_entropy<long double>();

    test_complex_hoyer_sparsity<std::complex<float>>();
    test_complex_hoyer_sparsity<std::complex<double>>();
    test_complex_hoyer_sparsity<std::complex<long double>>();
    test_complex_hoyer_sparsity<cpp_complex_50>();

    return boost::report_errors();
}
