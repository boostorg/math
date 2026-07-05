//  (C) Copyright Matt Borland 2022.
//  (C) Copyright Oleksandr Kornijcuk 2026.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <cstdint>
#include <cmath>
#include <vector>
#include <random>
#include <algorithm>
#include <utility>
#include <boost/math/statistics/chatterjee_correlation.hpp>
#include <boost/math/tools/random_vector.hpp>
#include <boost/math/constants/constants.hpp>
#include "math_unit_test.hpp"

// The Chatterjee correlation is invariant under:
// - Shuffles. (X_i, Y_i) -> (X_sigma(i), Y_sigma(i)), where sigma is a permutation.
// - Strictly monotone transformations: (X_i, Y_i) -> (f(X_i), g(Y_i)) where f' > 0 and g' > 0.
//

using boost::math::statistics::chatterjee_correlation;
using boost::math::statistics::chatterjee_correlation_mnn;

template <typename Real>
void properties()
{
    std::size_t vector_size = 256;
    std::mt19937_64 mt(123521);
    std::uniform_real_distribution<Real> unif(-1, 1);
    std::vector<Real> X(vector_size);
    std::vector<Real> Y(vector_size);

    for (std::size_t i = 0; i < vector_size; ++i)
    {
        X[i] = unif(mt);
        Y[i] = unif(mt);
    }
    std::sort(X.begin(), X.end());
    Real coeff1 = chatterjee_correlation(X, Y);
    // The minimum possible value of En(X, Y) is -1/2 + O(1/n)
    CHECK_GE(coeff1, Real(-0.5));
    CHECK_LE(coeff1, Real(1));

    // Now apply a monotone function to the data
    for (std::size_t i = 0; i < vector_size; ++i)
    {
        X[i] = Real(2.3)*X[i] - Real(7.3);
        Y[i] = Real(7.6)*Y[i] - Real(8.6);
    }
    auto coeff3 = chatterjee_correlation(X, Y);
    CHECK_EQUAL(coeff1, coeff3);

    // If there are no ties among the Yis, the maximum possible value of Xi(X, Y) is (n - 2)/(n + 1), which is attained if Yi = Xi for all i
    auto coeff = chatterjee_correlation(X, X);
    // These floating point numbers are computed by two different methods, so we can expect some floating point error:
    const auto n = X.size();
    CHECK_ULP_CLOSE(coeff, Real(n-2)/Real(n+1), 1);
    std::sort(Y.begin(), Y.end());
    coeff = chatterjee_correlation(Y, Y);
    CHECK_ULP_CLOSE(coeff, Real(n-2)/Real(n+1), 1);
}

template <typename Real>
void test_spots()
{
    // Rank Order: Result will be 1 - 3*3 / (4^2 - 1) = 1 - 9/15 = 0.6
    std::vector<Real> x = {1, 2, 3, 4};
    std::vector<Real> y = {1, 2, 3, 4};
    CHECK_ULP_CLOSE(chatterjee_correlation(x, y), 1 - Real(9)/15, 1);

    // Reverse rank order should be the same as above
    y = {4, 3, 2, 1};
    CHECK_ULP_CLOSE(chatterjee_correlation(x, y), 1 - Real(9)/15, 1);

    // Alternating order: 1 - 3*5 / (4^2 - 1) = 1 - 15/15 = 0
    y = {1, 3, 2, 4};
    CHECK_ULP_CLOSE(chatterjee_correlation(x, y), Real(0), 1);

    // All ties will yield quiet NaN
    y = {1, 1, 1, 1};
    CHECK_NAN(chatterjee_correlation(x, y));
}

// Closed forms for the M-NN statistic xi_{n,M} of Lin and Han (2021), Remark 2.5.
// These are exact (no external reference needed) and are used to validate the implementation.

// Y = f(X) with f strictly increasing: xi_{n,M} = 1 - (3(M+1)/4) / (n + (M+1)/4).
template <typename Real>
Real xi_mnn_increasing(std::size_t n, std::size_t M)
{
    const Real nr = static_cast<Real>(n);
    const Real Mr = static_cast<Real>(M);
    return Real(1) - (Real(3)*(Mr + 1)/4) / (nr + (Mr + 1)/4);
}

// Y = f(X) with f strictly decreasing:
// xi_{n,M} = 1 - 3(M+1)[(5n+1)/4 - (2M+1)/3] / ((n+1)(n + (M+1)/4)).
template <typename Real>
Real xi_mnn_decreasing(std::size_t n, std::size_t M)
{
    const Real nr = static_cast<Real>(n);
    const Real Mr = static_cast<Real>(M);
    return Real(1) - (Real(3)*(Mr + 1)*((Real(5)*nr + 1)/4 - (Real(2)*Mr + 1)/3))
                     / ((nr + 1)*(nr + (Mr + 1)/4));
}

template <typename Real>
void test_mnn_spots()
{
    // Exact small spot values, computed independently as rationals.
    std::vector<Real> x = {1, 2, 3, 4};
    std::vector<Real> y = {1, 2, 3, 4};

    // M = 2: -2 + 6*20 / (5 * 9.5) = 10/19
    CHECK_ULP_CLOSE(chatterjee_correlation_mnn(x, y, 2), Real(10)/19, 4);
    // M = 3: 2/5
    CHECK_ULP_CLOSE(chatterjee_correlation_mnn(x, y, 3), Real(2)/5, 4);

    // Reverse order, M = 2: -34/95. Note: unlike the M = 1 statistic, the M-NN statistic is NOT
    // invariant under reversal of Y for M > 1 (min{.,.} is not reversal-symmetric).
    y = {4, 3, 2, 1};
    CHECK_ULP_CLOSE(chatterjee_correlation_mnn(x, y, 2), Real(-34)/95, 4);

    // Alternating order, M = 2: 2/5
    y = {1, 3, 2, 4};
    CHECK_ULP_CLOSE(chatterjee_correlation_mnn(x, y, 2), Real(2)/5, 4);

    // n = 5 increasing, M = 2: 14/23
    std::vector<Real> x5 = {1, 2, 3, 4, 5};
    std::vector<Real> y5 = {1, 2, 3, 4, 5};
    CHECK_ULP_CLOSE(chatterjee_correlation_mnn(x5, y5, 2), Real(14)/23, 4);

    // Constant Y yields quiet NaN.
    y = {1, 1, 1, 1};
    CHECK_NAN(chatterjee_correlation_mnn(x, y, 2));
}

template <typename Real>
void test_mnn_extremal()
{
    // Build an X sorted ascending with no ties.
    const std::size_t n = 200;
    std::vector<Real> x(n);
    for (std::size_t i = 0; i < n; ++i)
    {
        x[i] = static_cast<Real>(i);
    }

    // Strictly increasing dependence Y = X.
    for (std::size_t M : {std::size_t(1), std::size_t(2), std::size_t(5), std::size_t(20), std::size_t(50)})
    {
        const Real got = chatterjee_correlation_mnn(x, x, M);
        CHECK_ULP_CLOSE(got, xi_mnn_increasing<Real>(n, M), 4);
    }

    // Strictly decreasing dependence Y = -X.
    std::vector<Real> y(n);
    for (std::size_t i = 0; i < n; ++i)
    {
        y[i] = -x[i];
    }
    for (std::size_t M : {std::size_t(1), std::size_t(2), std::size_t(5), std::size_t(20), std::size_t(50)})
    {
        const Real got = chatterjee_correlation_mnn(x, y, M);
        CHECK_ULP_CLOSE(got, xi_mnn_decreasing<Real>(n, M), 4);
    }
}

template <typename Real>
void test_mnn_properties()
{
    const std::size_t n = 256;
    std::mt19937_64 mt(987654);
    std::uniform_real_distribution<Real> unif(-1, 1);
    std::vector<Real> X(n);
    std::vector<Real> Y(n);
    for (std::size_t i = 0; i < n; ++i)
    {
        X[i] = unif(mt);
        Y[i] = unif(mt);
    }
    std::sort(X.begin(), X.end());

    const std::size_t M = 16;
    const Real coeff1 = chatterjee_correlation_mnn(X, Y, M);

    // The finite-sample range is [-1/2, 1] up to a bias of order M/n; use generous bounds.
    CHECK_GE(coeff1, Real(-1));
    CHECK_LE(coeff1, Real(1));

    // Invariance under strictly increasing transforms of X and Y.
    for (std::size_t i = 0; i < n; ++i)
    {
        X[i] = Real(2.3)*X[i] - Real(7.3);
        Y[i] = Real(7.6)*Y[i] - Real(8.6);
    }
    const Real coeff2 = chatterjee_correlation_mnn(X, Y, M);
    CHECK_EQUAL(coeff1, coeff2);
}

#ifdef BOOST_MATH_EXEC_COMPATIBLE

template <typename Real, typename ExecutionPolicy>
void test_threaded(ExecutionPolicy&& exec)
{
    std::vector<Real> x = boost::math::generate_random_vector<Real>(1024, 2);
    std::vector<Real> y = boost::math::generate_random_vector<Real>(1024, 1);

    std::sort(std::forward<ExecutionPolicy>(exec), x.begin(), x.end());

    auto seq_ans = chatterjee_correlation(x, y);
    auto par_ans = chatterjee_correlation(exec, x, y);

    CHECK_ULP_CLOSE(seq_ans, par_ans, 1);
};

template <typename Real, typename ExecutionPolicy>
void test_mnn_threaded(ExecutionPolicy&& exec)
{
    std::vector<Real> x = boost::math::generate_random_vector<Real>(1024, 2);
    std::vector<Real> y = boost::math::generate_random_vector<Real>(1024, 1);

    std::sort(std::forward<ExecutionPolicy>(exec), x.begin(), x.end());

    for (std::size_t M : {std::size_t(1), std::size_t(8), std::size_t(32)})
    {
        auto seq_ans = chatterjee_correlation_mnn(x, y, M);
        auto par_ans = chatterjee_correlation_mnn(exec, x, y, M);
        CHECK_ULP_CLOSE(seq_ans, par_ans, 1);
    }
};

#endif // BOOST_MATH_EXEC_COMPATIBLE

template <typename Real>
void test_paper()
{
    constexpr Real two_pi = boost::math::constants::two_pi<Real>();

    // Page 9 figure (a) y = x
    size_t seed = 3;
    std::vector<Real> x = boost::math::generate_random_uniform_vector<Real>(100, seed, -two_pi, two_pi);
    std::sort(x.begin(), x.end());
    auto result = chatterjee_correlation(x, x);
    CHECK_MOLLIFIED_CLOSE(result, Real(0.970), 0.005);

    // Page 9 figure (d) y = x^2
    std::vector<Real> y = x;
    for (auto& i : y)
    {
        i *= i;
    }

    result = chatterjee_correlation(x, y);
    CHECK_MOLLIFIED_CLOSE(result, Real(0.941), 0.005);

    // Page 9 figure (g) y = sin(x)
    for (std::size_t i {}; i < x.size(); ++i)
    {
        y[i] = std::sin(x[i]);
    }

    result = chatterjee_correlation(x, y);
    CHECK_MOLLIFIED_CLOSE(result, Real(0.885), 0.012);
}

int main(void)
{
    properties<float>();
    properties<double>();
    properties<long double>();

    test_spots<float>();
    test_spots<double>();
    test_spots<long double>();

    test_mnn_spots<float>();
    test_mnn_spots<double>();
    test_mnn_spots<long double>();

    test_mnn_extremal<float>();
    test_mnn_extremal<double>();
    test_mnn_extremal<long double>();

    test_mnn_properties<float>();
    test_mnn_properties<double>();
    test_mnn_properties<long double>();

    #ifdef BOOST_MATH_EXEC_COMPATIBLE

    test_threaded<float>(std::execution::par);
    test_threaded<double>(std::execution::par);
    test_threaded<long double>(std::execution::par);
    test_threaded<float>(std::execution::par_unseq);
    test_threaded<double>(std::execution::par_unseq);
    test_threaded<long double>(std::execution::par_unseq);

    test_mnn_threaded<float>(std::execution::par);
    test_mnn_threaded<double>(std::execution::par);
    test_mnn_threaded<long double>(std::execution::par);
    test_mnn_threaded<float>(std::execution::par_unseq);
    test_mnn_threaded<double>(std::execution::par_unseq);
    test_mnn_threaded<long double>(std::execution::par_unseq);

    #endif // BOOST_MATH_EXEC_COMPATIBLE

    test_paper<float>();
    test_paper<double>();
    test_paper<long double>();

    return boost::math::test::report_errors();
}
