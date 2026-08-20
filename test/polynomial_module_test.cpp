//  (C) Copyright Jeremy Murphy 2015.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// Module build check for boost/math/tools/polynomial.hpp. Consuming the class
// template and its free operators through `import boost.math;` confirms that
// polynomial construction, coefficient access, Horner evaluation, the
// arithmetic operators, quotient_remainder, pow, prime and integrate all
// resolve from a module consumer (along with the std::vector / std::pair /
// std::initializer_list machinery they rely on). Coefficients are stored in
// ascending order, so index 0 is the constant term. The known-good expected
// values are taken from test_polynomial.cpp (its double cases); the calculus
// checks use mathematically-certain derivative and integral coefficients.
// Only double and float are exercised.

#ifndef BOOST_MATH_BUILD_MODULE
#include <boost/math/tools/polynomial.hpp>
#else
import boost.math;
#endif

#include "math_unit_test.hpp"

#ifndef BOOST_MATH_BUILD_MODULE
#include <cstddef>
#endif

using boost::math::tools::polynomial;
using boost::math::tools::quotient_remainder;

// Construction, degree, indexing and Horner evaluation of
// p(x) = 10 - 6x - 4x^2 + 3x^3.
template <class Real>
void test_evaluation(const char* type_name)
{
    std::cout << "Testing polynomial construction and evaluation on type " << type_name << std::endl;

    const polynomial<Real> a {Real(10), Real(-6), Real(-4), Real(3)};

    CHECK_EQUAL(a.degree(), static_cast<std::size_t>(3));
    CHECK_EQUAL(Real(10), a[0]);
    CHECK_EQUAL(Real(-6), a[1]);
    CHECK_EQUAL(Real(-4), a[2]);
    CHECK_EQUAL(Real(3), a[3]);

    // Evaluate at exactly representable points.
    CHECK_ULP_CLOSE(Real(10), a(Real(0)), 3);          // constant term
    CHECK_ULP_CLOSE(Real(3), a(Real(1)), 3);           // 10 - 6 - 4 + 3
    CHECK_ULP_CLOSE(Real(6), a.evaluate(Real(2)), 3);  // 10 - 12 - 16 + 24
}

// (10 - 6x - 4x^2 + 3x^3) + (-2 + x) = 8 - 5x - 4x^2 + 3x^3.
template <class Real>
void test_addition(const char* type_name)
{
    std::cout << "Testing polynomial addition on type " << type_name << std::endl;

    const polynomial<Real> a {Real(10), Real(-6), Real(-4), Real(3)};
    const polynomial<Real> b {Real(-2), Real(1)};

    const auto sum = a + b;
    CHECK_EQUAL(sum.degree(), static_cast<std::size_t>(3));
    CHECK_EQUAL(Real(8), sum[0]);
    CHECK_EQUAL(Real(-5), sum[1]);
    CHECK_EQUAL(Real(-4), sum[2]);
    CHECK_EQUAL(Real(3), sum[3]);
}

// a * a = 100 - 120x - 44x^2 + 108x^3 - 20x^4 - 24x^5 + 9x^6, and pow agrees.
template <class Real>
void test_multiplication(const char* type_name)
{
    std::cout << "Testing polynomial multiplication and pow on type " << type_name << std::endl;

    const polynomial<Real> a {Real(10), Real(-6), Real(-4), Real(3)};
    const Real expected[7] {Real(100), Real(-120), Real(-44), Real(108), Real(-20), Real(-24), Real(9)};

    const auto asq = a * a;
    CHECK_EQUAL(asq.degree(), static_cast<std::size_t>(6));
    for (std::size_t i {0}; i < asq.size(); ++i)
    {
        CHECK_EQUAL(expected[i], asq[i]);
    }

    // pow(a, 2) reproduces the product; pow(a, 0) is the multiplicative identity.
    const auto apow2 = boost::math::tools::pow(a, 2);
    for (std::size_t i {0}; i < apow2.size(); ++i)
    {
        CHECK_EQUAL(asq[i], apow2[i]);
    }
    const auto one = boost::math::tools::pow(a, 0);
    CHECK_EQUAL(one.degree(), static_cast<std::size_t>(0));
    CHECK_EQUAL(Real(1), one[0]);
}

// (3x^3 - 4x^2 - 6x + 10) / (x - 2) = 3x^2 + 2x - 2, remainder 6.
template <class Real>
void test_division(const char* type_name)
{
    std::cout << "Testing polynomial division on type " << type_name << std::endl;

    const polynomial<Real> a {Real(10), Real(-6), Real(-4), Real(3)};
    const polynomial<Real> b {Real(-2), Real(1)};

    const auto result = quotient_remainder(a, b);
    const polynomial<Real>& q = result.first;
    const polynomial<Real>& r = result.second;

    CHECK_EQUAL(q.degree(), static_cast<std::size_t>(2));
    CHECK_EQUAL(Real(-2), q[0]);
    CHECK_EQUAL(Real(2), q[1]);
    CHECK_EQUAL(Real(3), q[2]);

    CHECK_EQUAL(r.degree(), static_cast<std::size_t>(0));
    CHECK_EQUAL(Real(6), r[0]);

    // q * b + r reconstructs the dividend exactly.
    const auto reconstructed = q * b + r;
    for (std::size_t i {0}; i < a.size(); ++i)
    {
        CHECK_EQUAL(a[i], reconstructed[i]);
    }

    // Division by a scalar divides each coefficient; the remainder is zero.
    const auto a_over_3 = a / Real(3);
    CHECK_ULP_CLOSE(Real(10) / Real(3), a_over_3[0], 3);
    CHECK_ULP_CLOSE(Real(-2), a_over_3[1], 3);
    CHECK_ULP_CLOSE(Real(-4) / Real(3), a_over_3[2], 3);
    CHECK_ULP_CLOSE(Real(1), a_over_3[3], 3);

    const auto a_mod_3 = a % Real(3);
    CHECK_EQUAL(a_mod_3.size(), static_cast<std::size_t>(0));

    // A polynomial divided by itself is the multiplicative identity.
    const auto a_over_a = a / a;
    CHECK_EQUAL(a_over_a.degree(), static_cast<std::size_t>(0));
    CHECK_EQUAL(Real(1), a_over_a[0]);
}

// p(x) = 1 + x + x^2 + x^3 + x^4 has derivative 1 + 2x + 3x^2 + 4x^3 and, with
// the constant chosen so P(0) = 0, integral coefficients 0, 1, 1/2, 1/3, 1/4, 1/5.
template <class Real>
void test_calculus(const char* type_name)
{
    std::cout << "Testing polynomial prime and integrate on type " << type_name << std::endl;

    const polynomial<Real> p {Real(1), Real(1), Real(1), Real(1), Real(1)};

    const auto dp = p.prime();
    CHECK_EQUAL(dp.degree(), static_cast<std::size_t>(3));
    for (std::size_t i {0}; i < dp.size(); ++i)
    {
        CHECK_EQUAL(Real(i + 1), dp[i]);
    }

    const auto P = p.integrate();
    CHECK_EQUAL(Real(0), P[0]);
    for (std::size_t i {1}; i < P.size(); ++i)
    {
        CHECK_ULP_CLOSE(Real(1) / Real(i), P[i], 3);
    }
}

int main()
{
    test_evaluation<float>("float");
    test_evaluation<double>("double");

    test_addition<float>("float");
    test_addition<double>("double");

    test_multiplication<float>("float");
    test_multiplication<double>("double");

    test_division<float>("float");
    test_division<double>("double");

    test_calculus<float>("float");
    test_calculus<double>("double");

    return boost::math::test::report_errors();
}
