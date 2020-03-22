/*
 * Copyright Nick Thompson, 2020
 * Use, modification and distribution are subject to the
 * Boost Software License, Version 1.0. (See accompanying file
 * LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
 */

#include "math_unit_test.hpp"
#include <numeric>
#include <utility>
#include <iomanip>
#include <iostream>
#include <cmath>
#include <boost/core/demangle.hpp>
#include <boost/hana/for_each.hpp>
#include <boost/hana/ext/std/integer_sequence.hpp>
#include <boost/math/quadrature/wavelet_transforms.hpp>
#include <boost/math/tools/minima.hpp>
#include <boost/math/quadrature/trapezoidal.hpp>

#ifdef BOOST_HAS_FLOAT128
#include <boost/multiprecision/float128.hpp>
using boost::multiprecision::float128;
#endif


using boost::math::constants::pi;
using boost::math::constants::root_two;
using boost::math::quadrature::daubechies_wavelet_transform;

template<typename Real, int p>
void test_wavelet_transform()
{
    std::cout << "Testing wavelet transform of " << p << " vanishing moment Daubechies wavelet on type " << boost::core::demangle(typeid(Real).name()) << "\n";
    auto psi = boost::math::daubechies_wavelet<Real, p>();

    auto abs_psi = [&psi](Real x) {
        return abs(psi(x));
    };
    auto [a, b] = psi.support();
    auto psil1 = boost::math::quadrature::trapezoidal(abs_psi, a, b);
    std::cout << "L1 norm of psi = " << psil1 << "\n";
    auto [x_min, psi_min] = boost::math::tools::brent_find_minima(psi, a, b, std::numeric_limits<Real>::digits10);
    std::cout << "Minimum value is " << psi_min <<  "  which occurs at x = " << x_min << "\n";
    auto neg_psi = [&](Real x) { return -psi(x); };
    auto [x_max, neg_psi_max] = boost::math::tools::brent_find_minima(neg_psi, a, b, std::numeric_limits<Real>::digits10);
    std::cout << "Maximum value of psi is " << -neg_psi_max << "\n";
    Real psi_sup_norm = std::max(abs(psi_min), std::abs(neg_psi_max));
    std::cout << "psi sup norm = " << psi_sup_norm << "\n";
    // An even function:
    auto f = [](Real x) {
        return std::exp(-abs(x))*std::cos(10*x);
    };
    Real fmax = 1;


    auto Wf = daubechies_wavelet_transform(f, psi);
    for (double s = 0; s < 10; s += 0.01)
    {
        Real w1 = Wf(s, 0.0);
        Real w2 = Wf(-s, 0.0);
        // Since f is an even function, we get w1 = w2:
        CHECK_ULP_CLOSE(w1, w2, 12);
        // Integral inequality:
        Real r1 = sqrt(abs(s))*fmax*psil1;
        //std::cout << "|w| = " << abs(w1) << ", ||f||_infty||psi||_1 = " << r1 << "\n";
        if(abs(w1) > r1)
        {
            std::cerr << " Integral inequality | int fg| <= ||f||_infty ||g||_infty is violated.\n";
        }
    }
    
}

int main()
{
    test_wavelet_transform<double, 8>();
    /*boost::hana::for_each(std::make_index_sequence<17>(), [&](auto i) {
        test_wavelet_transform<double, i+3>();
    });*/
    return boost::math::test::report_errors();
}
