//  Copyright Matt Borland, 2026
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
//  See: https://github.com/boostorg/math/issues/1668
//
//  digamma and polygamma reduced the trigonometric part of their reflection from 1 - x after it had been rounded.
//  Whenever 1 - x was inexact that gave a pole error at non-integer x, or lost digits at ordinary negative x.
//  The two terms of the reflection also cancel near each zero of digamma at negative x, where nothing was left of them.

#include <boost/math/special_functions/digamma.hpp>
#include <boost/math/special_functions/polygamma.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/math/concepts/real_concept.hpp>
#include <boost/math/policies/policy.hpp>
#include <boost/math/tools/precision.hpp>
#include <cmath>
#include <limits>
#include <stdexcept>

#include "math_unit_test.hpp"

// Each type has to be evaluated at its own precision: 1 - x is exact once x has been promoted.
using no_promotion = boost::math::policies::policy<boost::math::policies::promote_float<false>, boost::math::policies::promote_double<false>>;

template <class T>
void test()
{
    using boost::math::digamma;
    using boost::math::polygamma;
    using std::ldexp;

    const no_promotion pol {};
    const auto digits {boost::math::tools::digits<T>()};
    const auto eps {boost::math::tools::epsilon<T>()};
    const T tol {10 * eps};
    const auto pi {boost::math::constants::pi<T>()};
    const auto euler {boost::math::constants::euler<T>()};

    // The largest half integer: 1 - x rounds to an integer, where the cotangent and its even derivatives are zero.
    const T half {-(ldexp(T(1), digits - 1) - T(0.5))};
    const T pi_sq {pi * pi};
    CHECK_RELATIVE_ERROR(T((digits - 1) * boost::math::constants::ln_two<T>()), digamma(half, pol), tol);
    CHECK_RELATIVE_ERROR(T(-ldexp(T(1), -2 * (digits - 1))), polygamma(2, half, pol), tol);
    CHECK_RELATIVE_ERROR(T(2 * pi_sq * pi_sq), polygamma(3, half, pol), tol);
    CHECK_RELATIVE_ERROR(T(16 * pi_sq * pi_sq * pi_sq), polygamma(5, half, pol), tol);

    // The nearest values either side of a pole: 1 - x rounds to the integer, which used to be a pole error.
    const T below {-1 - eps};
    const T above {-2 + eps};
    CHECK_RELATIVE_ERROR(T(1 / eps + 1 - euler), digamma(below, pol), tol);
    CHECK_RELATIVE_ERROR(T(-1 / eps + T(1.5) - euler), digamma(above, pol), tol);

    // 1 - x rounds to 8.5, but the cotangent is not zero at x. To first order in eps, with digamma(8.5) and trigamma(8.5):
    const T near_half {-(T(7.5) + 4 * eps)};
    const T digamma_at_half {static_cast<T>(2.08009081757942012140262451060208803L)};
    const T trigamma_at_half {static_cast<T>(0.124838118918926021992296403327307295L)};
    CHECK_RELATIVE_ERROR(T(digamma_at_half + 4 * eps * (trigamma_at_half - pi_sq)), digamma(near_half, pol), tol);

    // The negative integers are still poles:
    CHECK_THROW(digamma(T(-1), pol), std::domain_error);
    CHECK_THROW(digamma(T(-2), pol), std::domain_error);
    CHECK_THROW(digamma(T(half - T(0.5)), pol), std::domain_error);

    if (digits >= 53)
    {
        // As reported: an ordinary double, with references from cpp_bin_float at 200 digits.
        // The odd powers of the cosine at n = 2, 4 and 24 see the sign of the reduction, and n = 24 is past the tabulated orders.
        const T x {static_cast<T>(-3.5520310993343904)};
        CHECK_RELATIVE_ERROR(static_cast<T>(0.883580268308463629616854778380608097L), digamma(x, pol), tol);
        CHECK_RELATIVE_ERROR(static_cast<T>(-10.5660939413826720964429875017171912L), polygamma(2, x, pol), tol);
        CHECK_RELATIVE_ERROR(static_cast<T>(216.419719412729058854230663349408100L), polygamma(3, x, pol), tol);
        CHECK_RELATIVE_ERROR(static_cast<T>(-863.396380879822761916680300227263654L), polygamma(4, x, pol), tol);
        CHECK_RELATIVE_ERROR(static_cast<T>(-322967669936316919447938532705641.366L), polygamma(24, x, pol), T(500 * eps));
    }
}

// float and double evaluate digamma in double-word arithmetic where its terms cancel.
template <class T, class Policy>
void test_near_zeros(const Policy& pol)
{
    using boost::math::digamma;

    const T tol {10 * boost::math::tools::epsilon<T>()};

    // The floats nearest some of the zeros, with references from cpp_bin_float at 200 digits.
    CHECK_RELATIVE_ERROR(static_cast<T>(-2.60254955029935875940781389197928384e-07L), digamma(static_cast<T>(-0.504083037f), pol), tol);
    CHECK_RELATIVE_ERROR(static_cast<T>(-1.42264490245067774754155064414771635e-07L), digamma(static_cast<T>(-1.57349849f), pol), tol);
    CHECK_RELATIVE_ERROR(static_cast<T>(-4.81110764929540668960511107411511341e-08L), digamma(static_cast<T>(-2.61072087f), pol), tol);
    CHECK_RELATIVE_ERROR(static_cast<T>(-2.28795078852302616779304775385416664e-06L), digamma(static_cast<T>(-7.68778849f), pol), tol);
    CHECK_RELATIVE_ERROR(static_cast<T>(-7.05501600621932775870983753486253531e-06L), digamma(static_cast<T>(-19.743082f), pol), tol);
    CHECK_RELATIVE_ERROR(static_cast<T>(-3.93884894889086536608283212640219964e-05L), digamma(static_cast<T>(-55.7892532f), pol), tol);
    CHECK_RELATIVE_ERROR(static_cast<T>(0.00133123500459704095919847177449172123L), digamma(static_cast<T>(-1000.86414f), pol), tol);

    // Further from a zero, where only a digit is lost:
    CHECK_RELATIVE_ERROR(static_cast<T>(-0.157662052808143238767682222467571768L), digamma(static_cast<T>(-2.625), pol), tol);
}

int main()
{
    test<float>();
    test<double>();
    // By default a float is evaluated as a double:
    test_near_zeros<float>(boost::math::policies::policy<>());
    test_near_zeros<float>(no_promotion());
    test_near_zeros<double>(no_promotion());
#ifndef BOOST_MATH_NO_LONG_DOUBLE_MATH_FUNCTIONS
    test<long double>();
    // And so does a long double that is a double:
    if (std::numeric_limits<long double>::digits == std::numeric_limits<double>::digits)
    {
        test_near_zeros<long double>(no_promotion());
    }
#endif
#ifndef BOOST_MATH_NO_REAL_CONCEPT_TESTS
    test<boost::math::concepts::real_concept>();
#endif

    return boost::math::test::report_errors();
}
