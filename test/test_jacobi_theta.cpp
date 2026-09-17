/*
 * Copyright Evan Miller, 2020
 * Use, modification and distribution are subject to the
 * Boost Software License, Version 1.0. (See accompanying file
 * LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
 */

#include <pch_light.hpp>
#include <boost/math/concepts/real_concept.hpp>
#include "test_jacobi_theta.hpp"

// Test file for the Jacobi Theta functions, a.k.a the four horsemen of the
// Jacobi elliptic integrals. At the moment only Wolfrma Alpha spot checks are
// used. We should generate extra-precise numbers with NTL::RR or some such.

void expected_results()
{
   //
   // Define the max and mean errors expected for
   // various compilers and platforms.
   //
   //
   add_expected_result(
      ".*",                          // compiler
      ".*",                          // stdlib
      ".*",                          // platform
      ".*",                  // test type(s)
      ".*Small Tau.*",      // test data group
      ".*", 1000, 200);  // test function

   add_expected_result(
      ".*",                          // compiler
      ".*",                          // stdlib
      ".*",                          // platform
      ".*",                  // test type(s)
      ".*WolframAlpha.*",      // test data group
      ".*", 60, 15);  // test function

   // Catch all cases come last:
   //
   add_expected_result(
      ".*",                          // compiler
      ".*",                          // stdlib
      ".*",                          // platform
      ".*",                  // test type(s)
      ".*",      // test data group
      ".*", 20, 5);  // test function
   //
   // Finish off by printing out the compiler/stdlib/platform names,
   // we do this to make it easier to mark up expected error rates.
   //
   std::cout << "Tests run with " << BOOST_COMPILER << ", "
      << BOOST_STDLIB << ", " << BOOST_PLATFORM << std::endl;
}

BOOST_AUTO_TEST_CASE( test_main )
{
    expected_results();
    BOOST_MATH_CONTROL_FP;
    BOOST_MATH_STD_USING

    using namespace boost::math;

    BOOST_CHECK_THROW(jacobi_theta1(0.0, 0.0), std::domain_error);
    BOOST_CHECK_THROW(jacobi_theta1(0.0, 1.0), std::domain_error);

    BOOST_CHECK_THROW(jacobi_theta2(0.0, 0.0), std::domain_error);
    BOOST_CHECK_THROW(jacobi_theta2(0.0, 1.0), std::domain_error);

    BOOST_CHECK_THROW(jacobi_theta3(0.0, 0.0), std::domain_error);
    BOOST_CHECK_THROW(jacobi_theta3(0.0, 1.0), std::domain_error);

    BOOST_CHECK_THROW(jacobi_theta4(0.0, 0.0), std::domain_error);
    BOOST_CHECK_THROW(jacobi_theta4(0.0, 1.0), std::domain_error);

    BOOST_CHECK_THROW(jacobi_theta1tau(0.0, 0.0), std::domain_error);
    BOOST_CHECK_THROW(jacobi_theta1tau(0.0, -1.0), std::domain_error);

    BOOST_CHECK_THROW(jacobi_theta2tau(0.0, 0.0), std::domain_error);
    BOOST_CHECK_THROW(jacobi_theta2tau(0.0, -1.0), std::domain_error);

    BOOST_CHECK_THROW(jacobi_theta3tau(0.0, 0.0), std::domain_error);
    BOOST_CHECK_THROW(jacobi_theta3tau(0.0, -1.0), std::domain_error);

    BOOST_CHECK_THROW(jacobi_theta4tau(0.0, 0.0), std::domain_error);
    BOOST_CHECK_THROW(jacobi_theta4tau(0.0, -1.0), std::domain_error);

    // Non-finite arguments must raise a domain error. Previously a NaN tau or q
    // (or a NaN or infinite z when tau < 1) never satisfied the convergence
    // test and the series loops ran forever.
    double nan = std::numeric_limits<double>::quiet_NaN();
    double inf = std::numeric_limits<double>::infinity();
    BOOST_CHECK_THROW(jacobi_theta1(1.0, nan), std::domain_error);
    BOOST_CHECK_THROW(jacobi_theta2(1.0, nan), std::domain_error);
    BOOST_CHECK_THROW(jacobi_theta3(1.0, nan), std::domain_error);
    BOOST_CHECK_THROW(jacobi_theta4(1.0, nan), std::domain_error);
    BOOST_CHECK_THROW(jacobi_theta3m1(1.0, nan), std::domain_error);
    BOOST_CHECK_THROW(jacobi_theta4m1(1.0, nan), std::domain_error);
    BOOST_CHECK_THROW(jacobi_theta1tau(1.0, nan), std::domain_error);
    BOOST_CHECK_THROW(jacobi_theta2tau(1.0, nan), std::domain_error);
    BOOST_CHECK_THROW(jacobi_theta3tau(1.0, nan), std::domain_error);
    BOOST_CHECK_THROW(jacobi_theta4tau(1.0, nan), std::domain_error);
    BOOST_CHECK_THROW(jacobi_theta3m1tau(1.0, nan), std::domain_error);
    BOOST_CHECK_THROW(jacobi_theta4m1tau(1.0, nan), std::domain_error);
    for (double tau : { 0.5, 1.5 }) {
        BOOST_CHECK_THROW(jacobi_theta1tau(nan, tau), std::domain_error);
        BOOST_CHECK_THROW(jacobi_theta2tau(nan, tau), std::domain_error);
        BOOST_CHECK_THROW(jacobi_theta3tau(nan, tau), std::domain_error);
        BOOST_CHECK_THROW(jacobi_theta4tau(nan, tau), std::domain_error);
        BOOST_CHECK_THROW(jacobi_theta3m1tau(nan, tau), std::domain_error);
        BOOST_CHECK_THROW(jacobi_theta4m1tau(nan, tau), std::domain_error);
        BOOST_CHECK_THROW(jacobi_theta1tau(inf, tau), std::domain_error);
        BOOST_CHECK_THROW(jacobi_theta2tau(-inf, tau), std::domain_error);
        BOOST_CHECK_THROW(jacobi_theta3tau(inf, tau), std::domain_error);
        BOOST_CHECK_THROW(jacobi_theta4tau(-inf, tau), std::domain_error);
        BOOST_CHECK_THROW(jacobi_theta3m1tau(inf, tau), std::domain_error);
        BOOST_CHECK_THROW(jacobi_theta4m1tau(-inf, tau), std::domain_error);
    }
    for (double q : { 0.5, 0.01 }) {
        BOOST_CHECK_THROW(jacobi_theta1(nan, q), std::domain_error);
        BOOST_CHECK_THROW(jacobi_theta2(inf, q), std::domain_error);
        BOOST_CHECK_THROW(jacobi_theta3(nan, q), std::domain_error);
        BOOST_CHECK_THROW(jacobi_theta4(-inf, q), std::domain_error);
        BOOST_CHECK_THROW(jacobi_theta3m1(nan, q), std::domain_error);
        BOOST_CHECK_THROW(jacobi_theta4m1(inf, q), std::domain_error);
    }

    double eps = std::numeric_limits<double>::epsilon();

    // theta1 with tau < 1 is evaluated via the imaginary transformation, where
    // it was previously computed as a difference of two nearly equal Gaussians
    // and lost precision proportional to 1/z (about 1e9 ulps at z = 1e-10).
    // Reference values: direct Fourier series in 40-digit arithmetic (mpmath).
    struct theta1_case { double z, tau, expected; };
    const theta1_case theta1_small_z[] = {
        { 1e-4, 0.5, 0.0001175932162452335795142692 },
        { 1e-8, 0.5, 1.175932162099660862653499e-8 },
        { -1e-8, 0.5, -1.175932162099660862653499e-8 },
        { 1e-10, 0.9, 9.76024299466749337952347e-11 },
        { 1e-6, 0.1, 2.455212638799910177735142e-8 },
        // Exercise the fold from [-pi, pi] into [-pi/2, pi/2] and the sign of the odd reflection
        { 3.0, 0.5, 0.1665975725928680619270942 },
        { -2.5, 0.5, -0.7534671707347227485840592 },
        { 1.5, 0.9, 0.987286730404979118338758 },
    };
    for (const theta1_case& c : theta1_small_z) {
        BOOST_CHECK_CLOSE_FRACTION(jacobi_theta1tau(c.z, c.tau), c.expected, 20 * eps);
    }

    // The q-parameterized functions previously went through tau = -log(q)/pi
    // and back through exp(), which amplified rounding by |log q| (about 100
    // ulps at q = 1e-200). They now evaluate q^n directly with pow() when q is
    // small enough for the direct series to be used.
    struct q_case { double q, theta1, theta2, theta3m1, theta4m1; };
    const q_case small_q[] = {
        { 1e-8, 0.01288435374475381934703355, 0.01529684374568976751542499, 1.650671229819356594481906e-8, -1.650671229819356594481904e-8 },
        { 1e-100, 1.288435374475382107345229e-25, 1.52968437456897685251172e-25, 1.650671229819356594481905e-100, -1.650671229819356594481905e-100 },
        { 1e-200, 1.288435374475382107345229e-50, 1.52968437456897685251172e-50, 1.650671229819356594481905e-200, -1.650671229819356594481905e-200 },
    };
    for (const q_case& c : small_q) {
        BOOST_CHECK_CLOSE_FRACTION(jacobi_theta1(0.7, c.q), c.theta1, 10 * eps);
        BOOST_CHECK_CLOSE_FRACTION(jacobi_theta2(0.7, c.q), c.theta2, 10 * eps);
        BOOST_CHECK_CLOSE_FRACTION(jacobi_theta3m1(0.3, c.q), c.theta3m1, 10 * eps);
        BOOST_CHECK_CLOSE_FRACTION(jacobi_theta4m1(0.3, c.q), c.theta4m1, 10 * eps);
        BOOST_CHECK_CLOSE_FRACTION(jacobi_theta3(0.3, c.q), 1 + c.theta3m1, 10 * eps);
        BOOST_CHECK_CLOSE_FRACTION(jacobi_theta4(0.3, c.q), 1 + c.theta4m1, 10 * eps);
    }

    for (double q=0.0078125; q<1.0; q += 0.0078125) { // = 1/128
        // The periodicity test shifts z by the rounded constant two_pi, which
        // differs from the true period by about eps. For large q the theta
        // functions are steep enough (their logarithmic derivative is of
        // order 1/tau = -pi/ln q) that this shift changes them by more than
        // the rounding of the evaluation itself, so allow for it.
        double periodicity_tol = 100 * eps + 4 * constants::pi<double>() * constants::pi<double>() * eps / -log(q);
        for (double z=-8.0; z<=8.0; z += 0.125) {
            test_periodicity(z, q, periodicity_tol);
            test_argument_translation(z, q, 100 * eps);
            test_sums_of_squares(z, q, 100 * eps);
            // The addition formula is complicated, cut it some extra slack
            test_addition_formulas(z, constants::ln_two<double>(), q, sqrt(sqrt(eps)));
            test_duplication_formula(z, q, 100 * eps);
            test_transformations_of_nome(z, q, 100 * eps);
            // Watson's identities subtract products of order 1 to leave a much
            // smaller remainder, so single-ulp differences in the individual
            // theta values are amplified by up to a couple of orders of magnitude.
            test_watsons_identities(z, 0.5, q, 200 * eps);
            test_landen_transformations(z, -log(q)/constants::pi<double>(), sqrt(eps));
            test_elliptic_functions(z, q, 5 * sqrt(eps));
        }
        test_elliptic_integrals(q, 10 * eps);
    }

    test_special_values(eps);

    for (double s=0.125; s<3.0; s+=0.125) {
        // The integrals sum thousands of theta values, so allow a few ulps
        test_mellin_transforms(2.0 + s, eps, 6 * eps);
        test_laplace_transforms(s, eps, 4 * eps);
    }

    // Laplace transforms at fixed z, for all four functions. The z values
    // include the small-z regime where theta1 with tau < 1 used to lose
    // precision, and the ranges cover both series branches in tau.
    for (double a : { 0.5, 1.0, 2.0, 5.0 }) {
        for (double z : { 1e-6, 1e-3, 0.1, 0.5, 1.0, 1.5, -0.7, 2.5, 3.0 }) {
            test_laplace_transforms_in_z(static_cast<float>(a), static_cast<float>(z),
                std::numeric_limits<float>::epsilon(), 25 * std::numeric_limits<float>::epsilon());
            test_laplace_transforms_in_z(a, z, eps, 25 * eps);
        }
    }

    test_spots(0.0F, "float");
    test_spots(0.0, "double");
#ifndef BOOST_MATH_NO_LONG_DOUBLE_MATH_FUNCTIONS
    test_spots(0.0L, "long double");
#ifndef BOOST_MATH_NO_REAL_CONCEPT_TESTS
    test_spots(concepts::real_concept(0), "real_concept");
#endif
#else
   std::cout << "<note>The long double tests have been disabled on this platform "
      "either because the long double overloads of the usual math functions are "
      "not available at all, or because they are too inaccurate for these tests "
      "to pass.</note>" << std::endl;
#endif
}
