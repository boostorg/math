/*
 * Copyright Nick Thompson, 2019
 * Use, modification and distribution are subject to the
 * Boost Software License, Version 1.0. (See accompanying file
 * LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
 */

#include "math_unit_test.hpp"
#include <numeric>
#include <utility>
#include <random>
#include <boost/core/demangle.hpp>
#include <boost/math/distributions/anderson_darling.hpp>
#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/cauchy.hpp>
#ifdef BOOST_HAS_FLOAT128
#include <boost/multiprecision/float128.hpp>
using boost::multiprecision::float128;
#endif

using boost::math::distributions::anderson_darling;

template<class Real>
void test_uniform()
{
    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    Real a = -0.5;
    Real b =  0.5;
    std::uniform_real_distribution<Real> dis(a,b);
    std::vector<Real> v(100000);
    for (auto & x : v) {
        x = dis(gen);
    }

    std::sort(v.begin(), v.end());

    auto CDF = [&a, &b](Real x) {
        if(x < a) {
            return Real(0);
        }
        if(x > b) {
            return Real(1);
        }
        return (x-a)/(b-a);
    };

    Real Asq = anderson_darling(v, CDF);
    std::cout << "Uniform A^2 = " << Asq << "\n";
}


template<class Real>
void test_normal()
{
    using std::sqrt;
    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    Real mu = 0.0;
    Real sigma =  1.0;
    std::normal_distribution<Real> dis(mu, sqrt(sigma));

    std::vector<Real> v(100000);
    for (auto & x : v) {
        x = dis(gen);
    }

    std::sort(v.begin(), v.end());

    auto normal = boost::math::normal_distribution<Real>(mu, sqrt(sigma));

    auto CDF = [&normal](Real x) { return boost::math::cdf(normal, x); };
    Real Asq = anderson_darling(v, CDF);
    std::cout << "Normal A^2 = " << Asq << "\n";
}

template<class Real>
void test_cauchy()
{
    using std::sqrt;
    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    Real x0 = 0.0;
    Real gamma =  1.0;
    std::cauchy_distribution<Real> dis(x0, gamma);

    std::vector<Real> v(100000);
    for (auto & x : v) {
        x = dis(gen);
    }

    std::sort(v.begin(), v.end());

    auto cauchy = boost::math::cauchy_distribution<Real>(x0, gamma);

    auto CDF = [&cauchy](Real x) { return boost::math::cdf(cauchy, x); };
    Real Asq = anderson_darling(v, CDF);
    std::cout << "Cauchy A^2 = " << Asq << "\n";
}


int main()
{
    test_uniform<double>();

    test_normal<double>();

    test_cauchy<double>();
    return boost::math::test::report_errors();
}
