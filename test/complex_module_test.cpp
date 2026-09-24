//  (C) Copyright John Maddock 2005.
//  (C) Copyright Matt Borland 2026.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// Module build check for the complex component <boost/math/complex.hpp>.
// Consuming boost.math as a named module, it exercises the exported entry
// points acos, asin, acosh, asinh, atan, atanh and fabs. Every expected value
// is mathematically certain: exact Pythagorean hypotenuses for fabs, the C99
// special cases at the origin, and standard exact arc-function values. These
// mirror the double/float spot checks in complex_test.cpp.

#ifndef BOOST_MATH_BUILD_MODULE
#include <boost/math/complex.hpp>
#else
import boost.math;
#endif

#ifndef BOOST_MATH_BUILD_MODULE
#include <complex>
#include <iostream>
#endif
#include "math_unit_test.hpp"

// The boost::math complex overloads are deprecated in favour of the C++11
// std:: versions but remain exported; silence the deprecation notice.
#if defined(_MSC_VER)
#  pragma warning(disable : 4996)
#elif defined(__clang__)
#  pragma clang diagnostic ignored "-Wdeprecated-declarations"
#elif defined(__GNUC__)
#  pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#endif

template <class Real>
void test_spots(const char* type_name)
{
    std::cout << "Testing complex functions for type " << type_name << std::endl;

    // Calls are fully qualified: an unqualified name plus a std::complex
    // argument would be ambiguous with the C++11 std:: overload via ADL.

    // Reference angles held in double precision (the module and, textually,
    // <boost/math/complex.hpp> both make boost::math::constants available).
    const double half_pi {boost::math::constants::half_pi<double>()};
    const double third_pi {boost::math::constants::third_pi<double>()};
    const double sixth_pi {boost::math::constants::sixth_pi<double>()};
    const double quarter_pi {boost::math::constants::quarter_pi<double>()};

    using ct = std::complex<Real>;
    const Real zero_tol {static_cast<Real>(1e-6)};

    // fabs(z) = hypot(re, im): exact Pythagorean triples.
    CHECK_ULP_CLOSE(5.0, boost::math::fabs(ct(3, 4)), 0);
    CHECK_ULP_CLOSE(13.0, boost::math::fabs(ct(5, 12)), 2);
    CHECK_EQUAL(static_cast<Real>(0), boost::math::fabs(ct(0, 0)));

    // acos(0) = pi/2 (C99 spot); acos(1/2) = pi/3. A real input stays real.
    ct r {boost::math::acos(ct(0, 0))};
    CHECK_ULP_CLOSE(half_pi, r.real(), 2);
    CHECK_ABSOLUTE_ERROR(0.0, r.imag(), zero_tol);
    r = boost::math::acos(ct(0.5, 0));
    CHECK_ULP_CLOSE(third_pi, r.real(), 4);
    CHECK_ABSOLUTE_ERROR(0.0, r.imag(), zero_tol);

    // asin(1/2) = pi/6; asin(0) = 0.
    r = boost::math::asin(ct(0.5, 0));
    CHECK_ULP_CLOSE(sixth_pi, r.real(), 4);
    CHECK_ABSOLUTE_ERROR(0.0, r.imag(), zero_tol);
    r = boost::math::asin(ct(0, 0));
    CHECK_ABSOLUTE_ERROR(0.0, r.real(), zero_tol);
    CHECK_ABSOLUTE_ERROR(0.0, r.imag(), zero_tol);

    // acosh(0) = i*pi/2 (C99 spot); asinh(0) = 0.
    r = boost::math::acosh(ct(0, 0));
    CHECK_ABSOLUTE_ERROR(0.0, r.real(), zero_tol);
    CHECK_ULP_CLOSE(half_pi, r.imag(), 2);
    r = boost::math::asinh(ct(0, 0));
    CHECK_ABSOLUTE_ERROR(0.0, r.real(), zero_tol);
    CHECK_ABSOLUTE_ERROR(0.0, r.imag(), zero_tol);

    // atan(1) = pi/4; atan(0) = 0; atanh(0) = 0.
    r = boost::math::atan(ct(1, 0));
    CHECK_ULP_CLOSE(quarter_pi, r.real(), 2);
    CHECK_ABSOLUTE_ERROR(0.0, r.imag(), zero_tol);
    r = boost::math::atan(ct(0, 0));
    CHECK_ABSOLUTE_ERROR(0.0, r.real(), zero_tol);
    CHECK_ABSOLUTE_ERROR(0.0, r.imag(), zero_tol);
    r = boost::math::atanh(ct(0, 0));
    CHECK_ABSOLUTE_ERROR(0.0, r.real(), zero_tol);
    CHECK_ABSOLUTE_ERROR(0.0, r.imag(), zero_tol);

    // Principal-branch identity: acos(z) + asin(z) = pi/2 for any z.
    const ct z {static_cast<Real>(0.5), static_cast<Real>(0.25)};
    const ct s {boost::math::acos(z) + boost::math::asin(z)};
    CHECK_ULP_CLOSE(half_pi, s.real(), 4);
    CHECK_ABSOLUTE_ERROR(0.0, s.imag(), zero_tol);
}

int main()
{
    test_spots<double>("double");
    test_spots<float>("float");

    return boost::math::test::report_errors();
}
