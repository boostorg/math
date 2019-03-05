/*
 * Copyright Nick Thompson, 2019
 * Use, modification and distribution are subject to the
 * Boost Software License, Version 1.0. (See accompanying file
 * LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
 */

#include "math_unit_test.hpp"
#include <numeric>
#include <boost/math/tools/condition_numbers.hpp>
#include <boost/math/special_functions/daubechies_scaling.hpp>
#include <boost/math/special_functions/daubechies_filters.hpp>
#include <boost/math/special_functions/detail/daubechies_scaling_integer_grid.hpp>
#include <boost/math/constants/constants.hpp>

#ifdef BOOST_HAS_FLOAT128
#include <boost/multiprecision/float128.hpp>
using boost::multiprecision::float128;
#endif


using boost::math::constants::pi;
using boost::math::constants::root_two;

// Mallat, Theorem 7.4, characterization number 3:
// A conjugate mirror filter has p vanishing moments iff h^{(n)}(pi) = 0 for 0 <= n < p.
template<class Real, unsigned p>
void test_daubechies_filters()
{
    Real tol = 3*std::numeric_limits<Real>::epsilon();
    using boost::math::daubechies_scaling_filter;
    using boost::math::daubechies_wavelet_filter;

    auto h = daubechies_scaling_filter<Real, p>();
    auto g = daubechies_wavelet_filter<Real, p>();

    auto inner = std::inner_product(h.begin(), h.end(), g.begin(), Real(0));
    CHECK_MOLLIFIED_CLOSE(0, inner, tol);

    // This is implied by Fourier transform of the two-scale dilatation equation;
    // If this doesn't hold, the infinite product for m_0 diverges.
    Real H0 = 0;
    for (size_t j = 0; j < h.size(); ++j)
    {
        H0 += h[j];
    }
    CHECK_MOLLIFIED_CLOSE(root_two<Real>(), H0, tol);

    // This is implied if we choose the scaling function to be an orthonormal basis of V0.
    Real scaling = 0;
    for (size_t j = 0; j < h.size(); ++j) {
      scaling += h[j]*h[j];
    }
    CHECK_MOLLIFIED_CLOSE(1, scaling, tol);

    using std::pow;
    // Daubechies wavelet of order p has p vanishing moments.
    // Unfortunately, the condition number of the sum is infinite.
    // Hence we must scale the tolerance by the summation condition number to ensure that we don't get spurious test failures.
    for (size_t k = 1; k < p; ++k)
    {
        Real hk = 0;
        Real abs_hk = 0;
        for (size_t n = 0; n < h.size(); ++n)
        {
            Real t = pow(n, k)*h[n];
            if (n & 1)
            {
                hk -= t;
            }
            else
            {
                hk += t;
            }
            abs_hk += abs(t);
        }
        // Multiply the tolerance by the condition number:
        Real cond = abs(hk) > 0 ? abs_hk/abs(hk) : 1/std::numeric_limits<Real>::epsilon();
        if (!CHECK_MOLLIFIED_CLOSE(0, hk, cond*tol))
        {
            std::cerr << "  The " << k << "th moment of the p = " << p << " filter did not vanish\n";
            std::cerr << "  Condition number = " << abs_hk/abs(hk) << "\n";
        }
    }

    // For the scaling function to be orthonormal to its integer translates,
    // sum h_k h_{k-2l} = \delta_{0,l}.
    // See Theoretical Numerical Analysis, Atkinson, Exercise 4.5.2.
    // This is the last condition we could test to ensure that the filters are correct,
    // but I'm not gonna bother because it's painful!
}

template<class Real1, class Real2, size_t p>
void test_filter_ulp_distance()
{
    using boost::math::daubechies_scaling_filter;
    auto h1 = daubechies_scaling_filter<Real1, p>();
    auto h2 = daubechies_scaling_filter<Real2, p>();

    for (size_t i = 0; i < h1.size(); ++i)
    {
        if(!CHECK_ULP_CLOSE(h1[i], h2[i], 0))
        {
            std::cerr << "  Index " << i << " at order " << p << " failed tolerance check\n";
        }
    }
}

template<class Real, unsigned p, unsigned order>
void test_integer_grid()
{
    using boost::math::detail::daubechies_scaling_integer_grid;
    using boost::math::tools::summation_condition_number;
    Real unit_roundoff = std::numeric_limits<Real>::epsilon()/2;
    auto grid = daubechies_scaling_integer_grid<Real, p, order>();

    if constexpr (order == 0) {
        auto cond = summation_condition_number<Real>(0);
        for (auto & x : grid) {
            cond += x;
        }
        CHECK_MOLLIFIED_CLOSE(1, cond.sum(), 2*cond.l1_norm()*unit_roundoff);
    }

    if constexpr (order == 1) {
        std::cout << "grid = {";

        auto cond = summation_condition_number<Real>(0);
        for (auto & x : grid) {
            std::cout << x << ", ";
            cond += x;
        }
        std::cout << "}\n";
        CHECK_MOLLIFIED_CLOSE(-1, cond.sum(), 2*cond.l1_norm()*unit_roundoff);
    }

}


int main()
{
    //auto phi = boost::math::daubechies_scaling<double, 4>();

    //std::cout << phi(2.3) << "\n";

    test_integer_grid<float, 2, 0>();
    test_integer_grid<float, 3, 0>();
    test_integer_grid<float, 4, 0>();
    test_integer_grid<float, 5, 0>();
    test_integer_grid<float, 6, 0>();
    test_integer_grid<float, 7, 0>();
    test_integer_grid<float, 8, 0>();
    test_integer_grid<float, 9, 0>();
    test_integer_grid<float, 10, 0>();
    test_integer_grid<float, 11, 0>();
    test_integer_grid<float, 12, 0>();
    test_integer_grid<float, 13, 0>();

    test_integer_grid<double, 2, 0>();
    test_integer_grid<double, 3, 0>();
    test_integer_grid<double, 4, 0>();
    test_integer_grid<double, 5, 0>();
    test_integer_grid<double, 6, 0>();
    test_integer_grid<double, 7, 0>();
    test_integer_grid<double, 8, 0>();
    test_integer_grid<double, 9, 0>();
    test_integer_grid<double, 10, 0>();
    test_integer_grid<double, 11, 0>();
    test_integer_grid<double, 12, 0>();
    test_integer_grid<double, 13, 0>();

    test_integer_grid<long double, 2, 0>();
    test_integer_grid<long double, 3, 0>();
    test_integer_grid<long double, 4, 0>();
    test_integer_grid<long double, 5, 0>();
    test_integer_grid<long double, 6, 0>();
    test_integer_grid<long double, 7, 0>();
    test_integer_grid<long double, 8, 0>();
    test_integer_grid<long double, 9, 0>();
    test_integer_grid<long double, 10, 0>();
    test_integer_grid<long double, 11, 0>();
    test_integer_grid<long double, 12, 0>();
    test_integer_grid<long double, 13, 0>();

    test_integer_grid<float128, 2, 0>();
    test_integer_grid<float128, 3, 0>();
    test_integer_grid<float128, 4, 0>();
    test_integer_grid<float128, 5, 0>();
    test_integer_grid<float128, 6, 0>();
    test_integer_grid<float128, 7, 0>();
    test_integer_grid<float128, 8, 0>();
    test_integer_grid<float128, 9, 0>();
    test_integer_grid<float128, 10, 0>();
    test_integer_grid<float128, 11, 0>();
    test_integer_grid<float128, 12, 0>();
    test_integer_grid<float128, 13, 0>();

    test_integer_grid<float, 3, 1>();
    /*test_integer_grid<float, 3, 1>();
    test_integer_grid<float, 4, 1>();
    test_integer_grid<float, 5, 1>();
    test_integer_grid<float, 6, 1>();
    test_integer_grid<float, 7, 1>();
    test_integer_grid<float, 8, 1>();
    test_integer_grid<float, 9, 1>();
    test_integer_grid<float, 10, 1>();
    test_integer_grid<float, 11, 1>();
    test_integer_grid<float, 12, 1>();
    test_integer_grid<float, 13, 1>();

    test_integer_grid<double, 2, 1>();
    test_integer_grid<double, 3, 1>();
    test_integer_grid<double, 4, 1>();
    test_integer_grid<double, 5, 1>();
    test_integer_grid<double, 6, 1>();
    test_integer_grid<double, 7, 1>();
    test_integer_grid<double, 8, 1>();
    test_integer_grid<double, 9, 1>();
    test_integer_grid<double, 10, 1>();
    test_integer_grid<double, 11, 1>();
    test_integer_grid<double, 12, 1>();
    test_integer_grid<double, 13, 1>();

    test_integer_grid<long double, 2, 1>();
    test_integer_grid<long double, 3, 1>();
    test_integer_grid<long double, 4, 1>();
    test_integer_grid<long double, 5, 1>();
    test_integer_grid<long double, 6, 1>();
    test_integer_grid<long double, 7, 1>();
    test_integer_grid<long double, 8, 1>();
    test_integer_grid<long double, 9, 1>();
    test_integer_grid<long double, 10, 1>();
    test_integer_grid<long double, 11, 1>();
    test_integer_grid<long double, 12, 1>();
    test_integer_grid<long double, 13, 1>();

    test_integer_grid<float128, 2, 1>();
    test_integer_grid<float128, 3, 1>();
    test_integer_grid<float128, 4, 1>();
    test_integer_grid<float128, 5, 1>();
    test_integer_grid<float128, 6, 1>();
    test_integer_grid<float128, 7, 1>();
    test_integer_grid<float128, 8, 1>();
    test_integer_grid<float128, 9, 1>();
    test_integer_grid<float128, 10, 1>();
    test_integer_grid<float128, 11, 1>();
    test_integer_grid<float128, 12, 1>();
    test_integer_grid<float128, 13, 1>();*/



    test_daubechies_filters<float, 1>();
    test_daubechies_filters<float, 2>();
    test_daubechies_filters<float, 3>();
    test_daubechies_filters<float, 4>();
    test_daubechies_filters<float, 5>();
    test_daubechies_filters<float, 6>();
    test_daubechies_filters<float, 7>();
    test_daubechies_filters<float, 8>();
    test_daubechies_filters<float, 9>();


    test_daubechies_filters<double, 1>();
    test_daubechies_filters<double, 2>();
    test_daubechies_filters<double, 3>();
    test_daubechies_filters<double, 4>();
    test_daubechies_filters<double, 5>();
    test_daubechies_filters<double, 6>();
    test_daubechies_filters<double, 7>();
    test_daubechies_filters<double, 8>();
    test_daubechies_filters<double, 9>();
    test_daubechies_filters<double, 10>();
    test_daubechies_filters<double, 11>();
    test_daubechies_filters<double, 12>();
    test_daubechies_filters<double, 13>();

    test_daubechies_filters<long double, 1>();
    test_daubechies_filters<long double, 2>();
    test_daubechies_filters<long double, 3>();
    test_daubechies_filters<long double, 4>();
    test_daubechies_filters<long double, 5>();
    test_daubechies_filters<long double, 6>();
    test_daubechies_filters<long double, 7>();
    test_daubechies_filters<long double, 8>();
    test_daubechies_filters<long double, 9>();
    test_daubechies_filters<long double, 10>();
    test_daubechies_filters<long double, 11>();
    test_daubechies_filters<long double, 12>();

    #ifdef BOOST_HAS_FLOAT128
    test_filter_ulp_distance<float128, long double, 1>();
    test_filter_ulp_distance<float128, long double, 2>();
    test_filter_ulp_distance<float128, long double, 3>();
    test_filter_ulp_distance<float128, long double, 4>();
    test_filter_ulp_distance<float128, long double, 5>();
    test_filter_ulp_distance<float128, long double, 6>();
    test_filter_ulp_distance<float128, long double, 7>();
    test_filter_ulp_distance<float128, long double, 8>();
    test_filter_ulp_distance<float128, long double, 9>();
    test_filter_ulp_distance<float128, long double, 10>();
    test_filter_ulp_distance<float128, long double, 11>();
    test_filter_ulp_distance<float128, long double, 12>();
    test_filter_ulp_distance<float128, long double, 13>();
    test_filter_ulp_distance<float128, long double, 14>();
    test_filter_ulp_distance<float128, long double, 15>();
    test_filter_ulp_distance<float128, long double, 16>();
    test_filter_ulp_distance<float128, long double, 17>();
    test_filter_ulp_distance<float128, long double, 18>();
    test_filter_ulp_distance<float128, long double, 19>();
    test_filter_ulp_distance<float128, long double, 20>();
    test_filter_ulp_distance<float128, long double, 21>();
    test_filter_ulp_distance<float128, long double, 22>();
    test_filter_ulp_distance<float128, long double, 23>();
    test_filter_ulp_distance<float128, long double, 24>();

    test_daubechies_filters<float128, 1>();
    test_daubechies_filters<float128, 2>();
    test_daubechies_filters<float128, 3>();
    test_daubechies_filters<float128, 4>();
    test_daubechies_filters<float128, 5>();
    test_daubechies_filters<float128, 6>();
    test_daubechies_filters<float128, 7>();
    test_daubechies_filters<float128, 8>();
    test_daubechies_filters<float128, 9>();
    test_daubechies_filters<float128, 10>();
    test_daubechies_filters<float128, 11>();
    test_daubechies_filters<float128, 12>();
    #endif

    return boost::math::test::report_errors();
}
