/*
 * Copyright Nick Thompson, 2026
 * Use, modification and distribution are subject to the
 * Boost Software License, Version 1.0. (See accompanying file
 * LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
 */

// Without threads, the grids of the Daubechies scaling functions and wavelets are computed one after another.
#define BOOST_MATH_DISABLE_THREADS
#include "math_unit_test.hpp"
#include <boost/math/special_functions/daubechies_scaling.hpp>
#include <boost/math/special_functions/daubechies_wavelet.hpp>
#include <chrono>
#include <future>

template<int p>
void test_sequential_construction()
{
    const int refinements = 4;
    const boost::math::daubechies_scaling<double, p> phi(refinements);
    const boost::math::daubechies_wavelet<double, p> psi(refinements);
    const auto phi_grid = boost::math::daubechies_scaling_dyadic_grid<double, p, 0>(refinements);
    const auto phi_prime_grid = boost::math::daubechies_scaling_dyadic_grid<double, p, 1>(refinements);
    const auto psi_grid = boost::math::daubechies_wavelet_dyadic_grid<double, p, 0>(refinements);
    const double h = 1.0 / (1 << refinements);
    // The last node is the end of the support, where the function is zero:
    for (size_t i = 0; i + 1 < phi_grid.size(); ++i)
    {
        // The interpolants pass through the grid, up to the rounding of their evaluation:
        CHECK_ABSOLUTE_ERROR(phi_grid[i], phi(i * h), 1e-12);
        CHECK_ABSOLUTE_ERROR(psi_grid[i], psi(psi.support().first + i * h), 1e-12);
        if constexpr (p > 2)
        {
            CHECK_ABSOLUTE_ERROR(phi_prime_grid[i], phi.prime(i * h), 1e-9);
        }
    }
}

void test_deferred()
{
    bool ran = false;
    auto f = boost::math::detail::daubechies_async([&ran]() { ran = true; return 1; });
    CHECK_EQUAL(true, f.wait_for(std::chrono::seconds(0)) == std::future_status::deferred);
    CHECK_EQUAL(false, ran);
    CHECK_EQUAL(1, f.get());
    CHECK_EQUAL(true, ran);
}

int main()
{
    test_deferred();
    test_sequential_construction<2>();
    test_sequential_construction<3>();
    test_sequential_construction<6>();
    test_sequential_construction<10>();
    return boost::math::test::report_errors();
}
