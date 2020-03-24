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
using boost::math::quadrature::trapezoidal;

template<typename Real, int p>
void test_wavelet_transform()
{
    std::cout << "Testing wavelet transform of " << p << " vanishing moment Daubechies wavelet on type " << boost::core::demangle(typeid(Real).name()) << "\n";
    auto psi = boost::math::daubechies_wavelet<Real, p>();

    auto abs_psi = [&psi](Real x) {
        return abs(psi(x));
    };
    auto [a, b] = psi.support();
    auto psil1 = trapezoidal(abs_psi, a, b);
    Real psi_sup_norm = 0;
    for (double x = a; x < b; x += 0.0001)
    {
        Real y = psi(x);
        if (std::abs(y) > psi_sup_norm)
        {
            psi_sup_norm = std::abs(y);
        }
    }
    std::cout << "psi sup norm = " << psi_sup_norm << "\n";
    // An even function:
    auto f = [](Real x) {
        return std::exp(-abs(x));
    };
    Real fmax = 1;
    Real fl2 = 1;
    Real fl1 = 2;
    std::cout << "||f||_1 = " << fl1 << "\n";
    std::cout << "||f||_2 = " << fl2 << "\n";

    auto Wf = daubechies_wavelet_transform(f, psi);
    for (double s = 0; s < 10; s += 0.01)
    {
        Real w1 = Wf(s, 0.0);
        Real w2 = Wf(-s, 0.0);
        // Since f is an even function, we get w1 = w2:
        CHECK_ULP_CLOSE(w1, w2, 12);
    }

    // The wavelet transform with respect to Daubechies wavelets 
    for (double s = -10; s < 10; s += 0.1)
    {
        for (double t = -10; t < 10; t+= 0.1)
        {
            Real w = Wf(s, t);
            // Integral inequality:
            Real r1 = sqrt(abs(s))*fmax*psil1;
            if(abs(w) > r1)
            {
                std::cerr << " Integral inequality | int fg| <= ||f||_infty ||psi||_1 is violated.\n";
            }
            if (s != 0)
            {
                Real r2 = fl1*psi_sup_norm/sqrt(abs(s));
                if(abs(w) > r2)
                {
                    std::cerr << " Integral inequality | int fg| <= ||f||_1 ||psi||_infty/sqrt(|s|) is violated.\n";
                    std::cerr << " Violation: " << abs(w) << " !<= " << r2 << "\n";
                }
                Real r3 = fmax*psil1/sqrt(abs(s));
                if(abs(w) > r3)
                {
                    std::cerr << " Integral inequality | int fg| <= ||f||_infty ||psi||_1/sqrt(|s|) is violated.\n";
                    std::cerr << " Computed = " << abs(w) << ", expected " << r3 << "\n";
                }
            }
            if (abs(w) > fl2)
            {
                std::cerr << "  Integral inequality |f psi_s,t| <= ||f||_2 ||psi||2 violated.\n";
            }
            Real r4 = sqrt(abs(s))*fl1*psi_sup_norm;
            if (abs(w) > r4)
            {
                std::cerr << "  Integral inequality |W[f](s,t)| <= sqrt(|s|)||f||_1 ||psi||_infty is violated.\n";
            }
            Real r5 = sqrt(abs(s))*fmax*psil1;
            if (abs(w) > r5)
            {
                std::cerr << "  Integral inequality |W[f](s,t)| <= sqrt(|s|)||f||_infty ||psi||_1 is violated.\n";
            }

        }
    }
}

int main()
{
    boost::hana::for_each(std::make_index_sequence<17>(), [&](auto i) {
        test_wavelet_transform<double, i+3>();
    });
    return boost::math::test::report_errors();
}
