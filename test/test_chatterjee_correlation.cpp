//  (C) Copyright Matt Borland 2022.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <cstdint>
#include <vector>
#include <random>
#include <algorithm>
#include <boost/math/statistics/chatterjee_correlation.hpp>
#include "math_unit_test.hpp"

// The Chatterjee correlation is invariant under:
// - Shuffles. (X_i, Y_i) -> (X_sigma(i), Y_sigma(i)), where sigma is a permutation.
// - Strictly monotone transformations: (X_i, Y_i) -> (f(X_i), g(Y_i)) where f' > 0 and g' > 0.
//

using boost::math::statistics::chatterjee_correlation;

template <typename Real>
void properties()
{
    std::mt19937_64 mt(123521);
    std::uniform_real_distribution<Real> unif(-1, 1);
    std::vector<Real> X(256);
    std::vector<Real> Y(256);

    for (std::size_t i = 0; i < X.size(); ++i) 
    {
        X[i] = unif(mt);
	    Y[i] = unif(mt);
    }
    std::sort(X.begin(), X.end());
    Real coeff1 = chatterjee_correlation(X, Y);
    // "it is not very hard to prove that the minimum possible value of ξn(X, Y) is −1/2 + O(1/n)"
    CHECK_GE(coeff1, Real(-0.5));
    CHECK_LE(coeff1, Real(1));

    // Now apply a monotone function to the data
    for (std::size_t i = 0; i < X.size(); ++i) 
    {
        X[i] = 2.3*X[i] - 7.3;
	    Y[i] = 7.6*Y[i] - 8.6;
    }
    auto coeff3 = chatterjee_correlation(X, Y);
    CHECK_EQUAL(coeff1, coeff3);

    // If there are no ties among the Yi’s, the maximum possible value of Xi(X, Y) is (n − 2)/(n + 1), which is attained if Yi = Xi for all i
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

int main(void)
{
    properties<float>();
    properties<double>();
    properties<long double>();

    test_spots<float>();
    test_spots<double>();
    test_spots<long double>();

    return boost::math::test::report_errors();
}
