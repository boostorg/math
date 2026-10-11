//  (C) Copyright Nick Thompson, John Maddock 2023.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// Module build check for boost::math::tools::evaluate_polynomial_estrin. It
// consumes the exported Estrin polynomial evaluator through `import boost.math;`
// and confirms all three overloads (general container, std::array, explicit
// scratch space) resolve from a module consumer. Estrin evaluates
// p(z) = c[0] + c[1] z + c[2] z^2 + ... with coefficients in ascending order.
// The assertions use only mathematically-certain values: integer coefficients
// and integer arguments keep every intermediate exactly representable, so the
// results are exact.

#ifndef BOOST_MATH_BUILD_MODULE
#include <boost/math/tools/estrin.hpp>
#else
import boost.math;
#endif

#include "math_unit_test.hpp"

#ifndef BOOST_MATH_BUILD_MODULE
#include <vector>
#include <array>
#endif

using boost::math::tools::evaluate_polynomial_estrin;

// General RandomAccessContainer overload (allocates its own scratch vector).
template <class Real>
void test_vector_overload()
{
    // p(z) = 2 + 3 z + 4 z^2.
    std::vector<Real> coeffs {Real(2), Real(3), Real(4)};

    CHECK_ULP_CLOSE(Real(2), evaluate_polynomial_estrin(coeffs, Real(0)), 0);   // p(0) = c0
    CHECK_ULP_CLOSE(Real(9), evaluate_polynomial_estrin(coeffs, Real(1)), 0);   // p(1) = sum(coeffs)
    CHECK_ULP_CLOSE(Real(24), evaluate_polynomial_estrin(coeffs, Real(2)), 0);  // 2 + 6 + 16
    CHECK_ULP_CLOSE(Real(3), evaluate_polynomial_estrin(coeffs, Real(-1)), 0);  // 2 - 3 + 4

    // A single coefficient is the constant polynomial p(z) = c0 for all z.
    std::vector<Real> constant {Real(7)};
    CHECK_ULP_CLOSE(Real(7), evaluate_polynomial_estrin(constant, Real(0)), 0);
    CHECK_ULP_CLOSE(Real(7), evaluate_polynomial_estrin(constant, Real(5)), 0);

    // No coefficients evaluates to zero for any argument.
    std::vector<Real> empty;
    CHECK_EQUAL(Real(0), evaluate_polynomial_estrin(empty, Real(3)));
}

// std::array overload: same value, no heap allocation for the scratch space.
template <class Real>
void test_array_overload()
{
    // p(z) = 1 + 2 z + 3 z^2 + 4 z^3.
    std::array<Real, 4> coeffs {Real(1), Real(2), Real(3), Real(4)};

    CHECK_ULP_CLOSE(Real(1), evaluate_polynomial_estrin(coeffs, Real(0)), 0);   // p(0) = c0
    CHECK_ULP_CLOSE(Real(10), evaluate_polynomial_estrin(coeffs, Real(1)), 0);  // 1 + 2 + 3 + 4
    CHECK_ULP_CLOSE(Real(49), evaluate_polynomial_estrin(coeffs, Real(2)), 0);  // 1 + 4 + 12 + 32
}

// Explicit scratch-space overload: the caller supplies working storage of at
// least (n + 1) / 2 elements whose value type matches the argument.
template <class Real>
void test_scratch_overload()
{
    // p(z) = 5 + z^3.
    std::vector<Real> coeffs {Real(5), Real(0), Real(0), Real(1)};
    std::vector<Real> scratch(2);

    CHECK_ULP_CLOSE(Real(5), evaluate_polynomial_estrin(coeffs, scratch, Real(0)), 0);   // p(0) = 5
    CHECK_ULP_CLOSE(Real(6), evaluate_polynomial_estrin(coeffs, scratch, Real(1)), 0);   // 5 + 1
    CHECK_ULP_CLOSE(Real(13), evaluate_polynomial_estrin(coeffs, scratch, Real(2)), 0);  // 5 + 8
}

int main()
{
    test_vector_overload<float>();
    test_vector_overload<double>();

    test_array_overload<float>();
    test_array_overload<double>();

    test_scratch_overload<float>();
    test_scratch_overload<double>();

    return boost::math::test::report_errors();
}
