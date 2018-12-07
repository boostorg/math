/*
 *  (C) Copyright Nick Thompson 2018.
 *  Use, modification and distribution are subject to the
 *  Boost Software License, Version 1.0. (See accompanying file
 *  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
 */

#include <vector>
#include <array>
#include <forward_list>
#include <boost/core/lightweight_test.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/math/tools/vector_functionals.hpp>
#include <boost/multiprecision/cpp_bin_float.hpp>
#include <boost/multiprecision/cpp_complex.hpp>

using boost::multiprecision::cpp_bin_float_50;
using boost::multiprecision::cpp_complex_50;

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
    boost::numeric::ublas::vector<double> w(7);
    for (size_t i = 0; i < w.size(); ++i)
    {
        w[i] = i+1;
    }
    mu = boost::math::tools::mean(w.cbegin(), w.cend());
    BOOST_TEST(abs(mu - 4) < tol);

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
    auto [mu1, sigma1_sq] = boost::math::tools::mean_and_population_variance(u.begin(), u.end());
    BOOST_TEST(abs(mu1 - 1) < tol);
    BOOST_TEST(abs(sigma1_sq) < tol);

    std::vector<Real> w{0,1,0,1,0,1,0,1};
    auto [mu2, sigma2_sq] = boost::math::tools::mean_and_population_variance(w.begin(), w.end());
    BOOST_TEST(abs(mu2 - 1.0/2.0) < tol);
    BOOST_TEST(abs(sigma2_sq - 1.0/4.0) < tol);
}

template<class Real>
void test_lp()
{
    Real tol = std::numeric_limits<Real>::epsilon();
    std::vector<std::complex<Real>> v{{1,0}, {0,0}, {0,0}};
    Real l3 = boost::math::tools::lp_norm(v.begin(), v.end(), 3);
    BOOST_TEST(abs(l3 - 1) < tol);

    std::vector<Real> u{1,0,0};
    l3 = boost::math::tools::lp_norm(u.begin(), u.end(), 3);
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

    std::vector<std::complex<Real>> w{{0,-8}, {1,1}, {3,2}};
    s = boost::math::tools::sup_norm(w.begin(), w.end());
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

}

template<class Real>
void test_l0_norm()
{
    std::vector<Real> v{0,0,1};
    size_t count = boost::math::tools::l0_norm(v.begin(), v.end());
    BOOST_TEST_EQ(count, 1);
}

int main()
{
    test_mean<float>();
    test_mean<double>();
    test_mean<long double>();
    test_mean<cpp_bin_float_50>();

    test_mean_and_population_variance<float>();
    test_mean_and_population_variance<double>();
    test_mean_and_population_variance<long double>();
    test_mean_and_population_variance<cpp_bin_float_50>();

    test_lp<float>();
    test_lp<double>();
    test_lp<long double>();

    test_total_variation<float>();
    test_total_variation<double>();
    test_total_variation<long double>();
    test_total_variation<cpp_bin_float_50>();

    test_sup_norm<float>();
    test_sup_norm<double>();
    test_sup_norm<long double>();

    test_gini_coefficient<float>();
    test_gini_coefficient<double>();
    test_gini_coefficient<long double>();

    test_absolute_gini_coefficient<float>();
    test_absolute_gini_coefficient<double>();
    test_absolute_gini_coefficient<long double>();

    test_l0_norm<float>();
    test_l0_norm<double>();
    test_l0_norm<long double>();

    return boost::report_errors();
}
