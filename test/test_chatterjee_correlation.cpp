//  (C) Copyright Matt Borland 2022.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <cstdint>
#include <vector>
#include <random>
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
    for (size_t i = 0; i < X.size(); ++i) {
        X[i] = unif(mt);
	Y[i] = unif(mt);
    }
    Real coeff1 = chatterjee_correlation(X, Y);
    // "it is not very hard to prove that the minimum possible value of ξn(X, Y ) is −1/2 + O(1/n)"
    CHECK_GE(coeff1, -0.5);
    CHECK_LE(coeff1, 1);
    // This is not a tautology, because X Y are modified by the call.
    Real coeff2 = chatterjee_correlation(X, Y);
    // Note that since the chatterjee_correlation is computed with ranks,
    // we should get *exact* equality.
    CHECK_EQUAL(coeff1, coeff2);
    // Now apply a monotone function to the data
    for (size_t i = 0; i < X.size(); ++i) {
        X[i] = 2.3*X[i] - 7.3;
	Y[i] = 7.6*Y[i] - 8.6;
    }
    auto coeff3 = chatterjee_correlation(X, Y);
    CHECK_EQUAL(coeff1, coeff2);

    // Shuffle invariance:
    mt.seed(18214);
    std::shuffle(X.begin(), X.end(), mt);
    // Of course it has to be the *same* shuffle:
    mt.seed(18214);
    std::shuffle(Y.begin(), Y.end(), mt);
    auto coeff4 = chatterjee_correlation(X, Y);
    CHECK_EQUAL(coeff1, coeff4);
    // If there are no ties among the Yi’s, the maximum possible value of Xi(X, Y) is (n − 2)/(n + 1), which is attained if Yi = Xi for all i
    auto coeff = chatterjee_correlation(X, X);
    // These floating point numbers are computed by two different methods, so we can expect some floating point error:
    CHECK_ULP_CLOSE(coeff, Real(n-2)/Real(n+1), 1);
    coeff = chatterjee_correlation(Y, Y);
    CHECK_ULP_CLOSE(coeff, Real(n-2)/Real(n+1), 1);
}


int main(void)
{
    properties<float>();
    properties<double>();
    properties<long double>();

    return boost::math::test::report_errors();
}
