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
#include <boost/math/tools/descriptive_statistics.hpp>
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

template<class Z>
void test_integer_mean()
{
    double tol = std::numeric_limits<double>::epsilon();
    std::vector<Z> v{1,2,3,4,5};
    double mu = boost::math::tools::mean(v);
    BOOST_TEST(abs(mu - 3) < tol);

    // Work with std::array?
    std::array<Z, 5> w{1,2,3,4,5};
    mu = boost::math::tools::mean(w);
    BOOST_TEST(abs(mu - 3) < tol);
}

template<class Real>
void test_mean()
{
    Real tol = std::numeric_limits<Real>::epsilon();
    std::vector<Real> v{1,2,3,4,5};
    Real mu = boost::math::tools::mean(v.begin(), v.end());
    BOOST_TEST(abs(mu - 3) < tol);

    // Does range call work?
    mu = boost::math::tools::mean(v);
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

    // Does range work?
    mu = boost::math::tools::mean(v);
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

    auto [mu3, sigma3_sq] = boost::math::tools::mean_and_population_variance(w);
    BOOST_TEST(abs(mu3 - 1.0/2.0) < tol);
    BOOST_TEST(abs(sigma3_sq - 1.0/4.0) < tol);

}

template<class Z>
void test_integer_mean_and_population_variance()
{
    double tol = std::numeric_limits<double>::epsilon();
    std::vector<Z> v{1,1,1,1,1,1};
    auto [mu, sigma_sq] = boost::math::tools::mean_and_population_variance(v);
    BOOST_TEST(abs(mu - 1) < tol);
    BOOST_TEST(abs(sigma_sq) < tol);
}

template<class Real>
void test_median()
{
    std::mt19937 g(12);
    std::vector<Real> v{1,2,3,4,5,6,7};

    Real m = boost::math::tools::median(v.begin(), v.end());
    BOOST_TEST_EQ(m, 4);

    std::shuffle(v.begin(), v.end(), g);
    // Does range call work?
    m = boost::math::tools::median(v);
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

    // Does it work with std::array?
    std::array<Real, 3> w{1,2,3};
    m = boost::math::tools::median(w);
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
    m = boost::math::tools::absolute_median(v);
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

    std::array<Real, 3> w{1, 2, -3};
    m = boost::math::tools::absolute_median(w);
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
    m = boost::math::tools::absolute_median(v);
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
void test_gini_coefficient()
{
    Real tol = std::numeric_limits<Real>::epsilon();
    std::vector<Real> v{1,0,0};
    Real gini = boost::math::tools::gini_coefficient(v.begin(), v.end());
    BOOST_TEST(abs(gini - 1) < tol);

    gini = boost::math::tools::gini_coefficient(v);
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

    std::array<Real, 3> w{0,0,0};
    gini = boost::math::tools::gini_coefficient(w);
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

    hs = boost::math::tools::hoyer_sparsity(v);
    BOOST_TEST(abs(hs - 1) < tol);

    // Does it work with constant iterators?
    hs = boost::math::tools::hoyer_sparsity(v.cbegin(), v.cend());
    BOOST_TEST(abs(hs - 1) < tol);

    v[0] = 1;
    v[1] = 1;
    v[2] = 1;
    hs = boost::math::tools::hoyer_sparsity(v.cbegin(), v.cend());
    BOOST_TEST(abs(hs) < tol);

    std::array<Real, 3> w{1,1,1};
    hs = boost::math::tools::hoyer_sparsity(w);
    BOOST_TEST(abs(hs) < tol);
}

template<class Z>
void test_integer_hoyer_sparsity()
{
    using std::sqrt;
    double tol = 5*std::numeric_limits<double>::epsilon();
    std::vector<Z> v{1,0,0};
    double hs = boost::math::tools::hoyer_sparsity(v);
    BOOST_TEST(abs(hs - 1) < tol);

    v[0] = 1;
    v[1] = 1;
    v[2] = 1;
    hs = boost::math::tools::hoyer_sparsity(v);
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

    hs = boost::math::tools::hoyer_sparsity(v);
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

    gini = boost::math::tools::absolute_gini_coefficient(v);
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
void test_shannon_entropy()
{
    Real tol = 100*std::numeric_limits<Real>::epsilon();
    using boost::math::constants::half;
    using boost::math::constants::ln_two;
    std::vector<Real> v(30, half<Real>());
    Real Hs = boost::math::tools::shannon_entropy(v.begin(), v.end());
    Real expected = v.size()*ln_two<Real>()/2;
    BOOST_TEST(abs(Hs - expected) < tol*expected);

    Hs = boost::math::tools::shannon_entropy(v);
    BOOST_TEST(abs(Hs - expected) < tol*expected);

    std::array<Real, 3> w{half<Real>(), half<Real>(), half<Real>()};
    Hs = boost::math::tools::shannon_entropy(w);
    expected = 3*ln_two<Real>()/2;
    BOOST_TEST(abs(Hs - expected) < tol*expected);
}


int main()
{
    test_integer_mean<uint8_t>();
    test_integer_mean<int8_t>();
    test_integer_mean<int>();

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

    test_integer_mean_and_population_variance<int>();

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

    test_integer_hoyer_sparsity<int>();

    test_shannon_entropy<float>();
    test_shannon_entropy<double>();
    test_shannon_entropy<long double>();

    test_complex_hoyer_sparsity<std::complex<float>>();
    test_complex_hoyer_sparsity<std::complex<double>>();
    test_complex_hoyer_sparsity<std::complex<long double>>();
    test_complex_hoyer_sparsity<cpp_complex_50>();

    return boost::report_errors();
}
