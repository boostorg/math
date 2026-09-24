/*
 *  (C) Copyright Nick Thompson 2018.
 *  (C) Copyright Matt Borland 2026.
 *  Use, modification and distribution are subject to the
 *  Boost Software License, Version 1.0. (See accompanying file
 *  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
 *
 *  Module build check for the boost::math::statistics signal-statistics entry
 *  points (absolute_gini_coefficient, sample_absolute_gini_coefficient,
 *  hoyer_sparsity, oracle_snr, oracle_snr_db). Building and evaluating them here
 *  forces the exported templates to be instantiated and resolved from a module
 *  consumer. Data and expected values are copied from signal_statistics_test.cpp;
 *  only the double/float, mathematically-certain results are checked (its
 *  multiprecision, complex, integer and long double paths are omitted).
 */

#ifndef BOOST_MATH_BUILD_MODULE
#include <boost/math/statistics/signal_statistics.hpp>
#else
import boost.math;
#endif

#ifndef BOOST_MATH_BUILD_MODULE
#include <vector>
#include <array>
#include <limits>
#include <cstddef>
#endif
#include "math_unit_test.hpp"

// The absolute Gini coefficient of {-1, 0, 0} is the population value 2/3, and
// the sample value n/(n-1) times that, which is 1 for n = 3. A vector whose
// entries all share a single magnitude ({1, -1, 1}) has coefficient 0.
template <class Real>
void test_absolute_gini_coefficient()
{
    using boost::math::statistics::absolute_gini_coefficient;
    using boost::math::statistics::sample_absolute_gini_coefficient;

    const Real tol = std::numeric_limits<Real>::epsilon();

    // Both forms sort the range in place; the final coefficient is unchanged by order.
    std::vector<Real> v{static_cast<Real>(-1), static_cast<Real>(0), static_cast<Real>(0)};
    Real gini = sample_absolute_gini_coefficient(v.begin(), v.end());
    CHECK_ULP_CLOSE(static_cast<Real>(1), gini, 4);

    gini = absolute_gini_coefficient(v);
    CHECK_ULP_CLOSE(static_cast<Real>(2) / static_cast<Real>(3), gini, 4);

    std::vector<Real> w{static_cast<Real>(1), static_cast<Real>(-1), static_cast<Real>(1)};
    gini = absolute_gini_coefficient(w.begin(), w.end());
    CHECK_ABSOLUTE_ERROR(static_cast<Real>(0), gini, tol);
    gini = sample_absolute_gini_coefficient(w.begin(), w.end());
    CHECK_ABSOLUTE_ERROR(static_cast<Real>(0), gini, tol);
}

// A maximally sparse vector (one nonzero of three) has Hoyer sparsity 1; a flat
// vector (all entries equal) has Hoyer sparsity 0.
template <class Real>
void test_hoyer_sparsity()
{
    using boost::math::statistics::hoyer_sparsity;

    const Real tol = 5 * std::numeric_limits<Real>::epsilon();

    std::vector<Real> v{static_cast<Real>(1), static_cast<Real>(0), static_cast<Real>(0)};
    CHECK_ULP_CLOSE(static_cast<Real>(1), hoyer_sparsity(v), 4);

    std::array<Real, 3> w{static_cast<Real>(1), static_cast<Real>(1), static_cast<Real>(1)};
    CHECK_ABSOLUTE_ERROR(static_cast<Real>(0), hoyer_sparsity(w), tol);
}

// A length-100 constant signal perturbed by 1 in a single sample has
// ||signal||^2 = 100 and ||noise||^2 = 1, so the SNR is exactly 100 and the
// SNR in decibels is 10 log10(100) = 20.
template <class Real>
void test_oracle_snr()
{
    using boost::math::statistics::oracle_snr;
    using boost::math::statistics::oracle_snr_db;

    const std::size_t length = 100;
    std::vector<Real> signal(length, static_cast<Real>(1));
    std::vector<Real> noisy_signal = signal;
    noisy_signal[0] += static_cast<Real>(1);

    Real snr = oracle_snr(signal, noisy_signal);
    CHECK_ULP_CLOSE(static_cast<Real>(length), snr, 4);

    Real snr_db = oracle_snr_db(signal, noisy_signal);
    CHECK_ULP_CLOSE(static_cast<Real>(20), snr_db, 16);
}

int main()
{
    test_absolute_gini_coefficient<double>();
    test_hoyer_sparsity<double>();
    test_oracle_snr<double>();

    test_absolute_gini_coefficient<float>();
    test_hoyer_sparsity<float>();
    test_oracle_snr<float>();

    return boost::math::test::report_errors();
}
