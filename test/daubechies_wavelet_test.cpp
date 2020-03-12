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
#include <random>
#include <boost/core/demangle.hpp>
#include <boost/hana/for_each.hpp>
#include <boost/hana/ext/std/integer_sequence.hpp>
#include <boost/math/tools/condition_numbers.hpp>
#include <boost/math/special_functions/daubechies_wavelet.hpp>
#include <boost/math/quadrature/trapezoidal.hpp>

#ifdef BOOST_HAS_FLOAT128
#include <boost/multiprecision/float128.hpp>
using boost::multiprecision::float128;
#endif


using boost::math::constants::pi;
using boost::math::constants::root_two;


template<class Real>
void test_wavelet_dyadic_grid()
{
    std::cout << "Testing wavelet dyadic grid on type " << boost::core::demangle(typeid(Real).name()) << "\n";
    auto f = [&](auto i)
    {
        auto phijk = boost::math::daubechies_scaling_dyadic_grid<Real, i+2, 0>(0);
        auto phik = boost::math::detail::daubechies_scaling_integer_grid<Real, i+2, 0>();
        assert(phik.size() == phijk.size());

        for (size_t k = 0; k < phik.size(); ++k)
        {
            CHECK_ULP_CLOSE(phik[k], phijk[k], 0);
        }

        for (int64_t j = 1; j < 10; ++j)
        {
            phijk = boost::math::daubechies_scaling_dyadic_grid<Real, i+2, 0>(j);
            phik = boost::math::detail::daubechies_scaling_integer_grid<Real, i+2, 0>();
            for (int64_t l = 0; l < static_cast<int64_t>(phik.size()); ++l)
            {
                CHECK_ULP_CLOSE(phik[l], phijk[l*(int64_t(1)<<j)], 0);
            }

            // This test is from Daubechies, Ten Lectures on Wavelets, Ch 7 "More About Compactly Supported Wavelets",
            // page 245: \forall y \in \mathbb{R}, \sum_{n \in \mathbb{Z}} \phi(y+n) = 1
            for (size_t k = 1; k < j; ++k)
            {
                auto cond = boost::math::tools::summation_condition_number<Real>(0);
                for (int64_t l = 0; l < static_cast<int64_t>(phik.size()); ++l)
                {
                    int64_t idx = l*(int64_t(1)<<j) + k;
                    if (idx < phijk.size())
                    {
                        cond += phijk[idx];
                    }
                }
                CHECK_MOLLIFIED_CLOSE(Real(1), cond.sum(), 10*cond()*std::numeric_limits<Real>::epsilon());
            }
        }
    };

    boost::hana::for_each(std::make_index_sequence<18>(), f);
}


template<typename Real, int p>
void test_quadratures()
{
    using boost::math::quadrature::trapezoidal;
    auto psi = boost::math::daubechies_wavelet<Real, p>();
        
    Real tol = std::numeric_limits<Real>::epsilon();
    Real error_estimate = std::numeric_limits<Real>::quiet_NaN();
    Real L1 = std::numeric_limits<Real>::quiet_NaN();
    auto [a, b] = psi.support();
    CHECK_ULP_CLOSE(Real(-p+1), a, 0);
    CHECK_ULP_CLOSE(Real(p), b, 0);
    // A wavelet is a function of zero average; ensure the quadrature over its support is zero.
    /*Real Q = trapezoidal(psi, a, b, tol, 15, &error_estimate, &L1);
    if (!CHECK_MOLLIFIED_CLOSE(Real(0), Q, Real(0.0001)))
    {
        std::cerr << "  Quadrature of " << p << " vanishing moment wavelet does not vanish.\n";
    }*/
    // psi is orthogonal to its integer translates: \int \psi(x-k) \psi(x) \, \mathrm{d}x = 0
    // psi has L2 norm 1:

    // g_n = 1/sqrt(2) <psi(t/2), phi(t-n)> (Mallat, 7.55)
}

int main()
{
    boost::hana::for_each(std::make_index_sequence<18>(), [&](auto i){
      test_quadratures<float, i+2>();
      test_quadratures<double, i+2>();
    });


    test_wavelet_dyadic_grid<float>();
    test_wavelet_dyadic_grid<double>();
    test_wavelet_dyadic_grid<long double>();
    #ifdef BOOST_HAS_FLOAT128
    test_wavelet_dyadic_grid<float128>();
    #endif

    return boost::math::test::report_errors();
}
