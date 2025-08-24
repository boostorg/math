/*
 * Copyright Nick Thompson, 2024
 * Use, modification and distribution are subject to the
 * Boost Software License, Version 1.0. (See accompanying file
 * LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
 */
#include <vector>
#include <iostream>
#include <list>
#include <random>
#include <cmath>
#include <complex>
#include <utility>
#include <limits>
#include <algorithm>
#include <boost/math/tools/polynomial.hpp>
#include "math_unit_test.hpp"
using boost::math::tools::polynomial;

#if __has_include(<Eigen/Eigenvalues>)

void test_random_coefficients() {
    std::random_device rd;
    uint32_t seed = rd(); 
    std::mt19937_64 mt(seed);
    std::uniform_real_distribution<double> unif(-1, 1);
    std::size_t n = seed % 3 + 3;
    std::vector<double> coeffs(n, std::numeric_limits<double>::quiet_NaN());
    for (std::size_t i = 0; i < coeffs.size(); ++i) {
        coeffs[i] = unif(mt);
    }
    coeffs[coeffs.size() -1] = 1.0;
    auto p = polynomial(std::move(coeffs));
    auto roots = p.roots();
    CHECK_EQUAL(roots.size(), p.size() - 1);
    std::complex<double> root_product = -1;
    std::complex<double> root_sum = 0.0;
    for (auto const & root : roots) {
        root_product *= static_cast<std::complex<double>>(root);
        root_sum += static_cast<std::complex<double>>(root);
    }
    if (p.size() & 1) {
       root_product *= -1;
    }
    CHECK_ULP_CLOSE(root_product.real(), p[0], 1000);
    CHECK_LE(root_product.imag(), 1e-6);

    CHECK_LE(root_sum.imag(), 1e-7);
    CHECK_ULP_CLOSE(root_sum.real(), -p[p.size() - 2], 1000);
}



void test_wilkinson_polynomial() {
    // CoefficientList[Expand[(x - 1)*(x - 2)*(x - 3)*(x - 4)*(x - 5)*(x - 6)*(x - 7)*(x - 8)*(x - 9)*(x - 10)], x]
    std::vector<float> coeffs{3628800.0, -10628640.0, 12753576.0, -8409500.0, 3416930.0, -902055.0, 157773.0, -18150.0, 1320.0, -55.0 ,1.0};
    auto p = polynomial(std::move(coeffs));
    auto roots = p.roots();
    CHECK_EQUAL(roots.size(), p.size() - 1);
    std::complex<double> root_product = -1;
    std::complex<double> root_sum = 0.0;
    for (auto const & root : roots) {
        root_product *= static_cast<std::complex<double>>(root);
        root_sum += static_cast<std::complex<double>>(root);
    }
    if (p.size() & 1) {
       root_product *= -1;
    }
    CHECK_ABSOLUTE_ERROR(root_product.real(), double(p[0]), double(1e-3*p[0]));
    CHECK_LE(root_product.imag(), 1e-8);

    CHECK_LE(root_sum.imag(), 1e-8);
    CHECK_ABSOLUTE_ERROR(root_sum.real(), -double(p[p.size() - 2]), 1e-5);

    std::complex<double> c = 0.0;
    for (std::size_t i = 0; i < roots.size(); ++i) {
        auto ri = static_cast<std::complex<double>>(roots[i]);
        for (std::size_t j = i + 1; j < roots.size(); ++j) {
            c += ri*static_cast<std::complex<double>>(roots[j]);
        }
    }
    CHECK_ULP_CLOSE(p[p.size()-3], static_cast<float>(c.real()), 10);
    CHECK_ABSOLUTE_ERROR(0.0, c.imag(), 1e-8);

}

template<typename T>
void test_singular_companion()
{
    std::vector<T> coeffs{0.0, 0.0, 1.0}; 
    auto p = polynomial(std::move(coeffs));
    auto roots = p.roots();
    CHECK_EQUAL(roots.size(), p.size() - 1);
    for (std::size_t i = 0; i < roots.size() - 1; ++i) {
        CHECK_ABSOLUTE_ERROR(T(0), roots[i].real(), std::numeric_limits<T>::epsilon());
        CHECK_ABSOLUTE_ERROR(T(0), roots[i].imag(), std::numeric_limits<T>::epsilon());
    }
}


int main()
{
    test_random_coefficients();
    test_singular_companion<float>();
    test_singular_companion<double>();
    test_wilkinson_polynomial();
    return boost::math::test::report_errors();
}
#endif
