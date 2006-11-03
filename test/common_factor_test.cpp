//  Boost GCD & LCM common_factor.hpp test program  --------------------------//

//  (C) Copyright Daryle Walker 2001, 2006.
//  Distributed under the Boost Software License, Version 1.0. (See
//  accompanying file LICENSE_1_0.txt or copy at
//  http://www.boost.org/LICENSE_1_0.txt)

//  See http://www.boost.org for most recent version including documentation.

//  Revision History
//  02 Nov 2006  Change to Boost.Test's unit test system
//  07 Nov 2001  Initial version (Daryle Walker)

#define BOOST_TEST_MAIN  "Boost.Math GCD & LCM unit tests"

#include <boost/config.hpp>              // for BOOST_MSVC
#include <boost/math/common_factor.hpp>  // for boost::math::gcd, etc.
#include <boost/mpl/list.hpp>            // for boost::mpl::list
#include <boost/test/unit_test.hpp>
#include <boost/test/test_case_template.hpp>


namespace {

// TODO: add new signed and unsigned types here, with and without numeric_limits
// specialization; add polynominal/non-real type; especially after any switch to
// the binary-GCD algorithm for built-in types

// Various types to test with each GCD/LCM
typedef ::boost::mpl::list<short, int, long>  builtin_signed_test_types;
typedef ::boost::mpl::list<unsigned short, unsigned, unsigned long>
  builtin_unsigned_test_types;

}  // namespace


// GCD tests
BOOST_AUTO_TEST_SUITE( gcd_test_suite )

// GCD on built-in signed integer types
BOOST_AUTO_TEST_CASE_TEMPLATE( gcd_int_test, T, builtin_signed_test_types )
{
#ifndef BOOST_MSVC
    using boost::math::gcd;
#else
    using namespace boost::math;
#endif

    // Originally from Boost.Rational tests
    BOOST_CHECK_EQUAL( gcd<T>(  1,  -1), static_cast<T>( 1) );
    BOOST_CHECK_EQUAL( gcd<T>( -1,   1), static_cast<T>( 1) );
    BOOST_CHECK_EQUAL( gcd<T>(  1,   1), static_cast<T>( 1) );
    BOOST_CHECK_EQUAL( gcd<T>( -1,  -1), static_cast<T>( 1) );
    BOOST_CHECK_EQUAL( gcd<T>(  0,   0), static_cast<T>( 0) );
    BOOST_CHECK_EQUAL( gcd<T>(  7,   0), static_cast<T>( 7) );
    BOOST_CHECK_EQUAL( gcd<T>(  0,   9), static_cast<T>( 9) );
    BOOST_CHECK_EQUAL( gcd<T>( -7,   0), static_cast<T>( 7) );
    BOOST_CHECK_EQUAL( gcd<T>(  0,  -9), static_cast<T>( 9) );
    BOOST_CHECK_EQUAL( gcd<T>( 42,  30), static_cast<T>( 6) );
    BOOST_CHECK_EQUAL( gcd<T>(  6,  -9), static_cast<T>( 3) );
    BOOST_CHECK_EQUAL( gcd<T>(-10, -10), static_cast<T>(10) );
    BOOST_CHECK_EQUAL( gcd<T>(-25, -10), static_cast<T>( 5) );
    BOOST_CHECK_EQUAL( gcd<T>(  3,   7), static_cast<T>( 1) );
    BOOST_CHECK_EQUAL( gcd<T>(  8,   9), static_cast<T>( 1) );
    BOOST_CHECK_EQUAL( gcd<T>(  7,  49), static_cast<T>( 7) );
}

// GCD on built-in unsigned integer types
BOOST_AUTO_TEST_CASE_TEMPLATE(gcd_unsigned_test, T, builtin_unsigned_test_types)
{
#ifndef BOOST_MSVC
    using boost::math::gcd;
#else
    using namespace boost::math;
#endif

    BOOST_CHECK_EQUAL( gcd<T>(  1u,   1u), static_cast<T>( 1u) );
    BOOST_CHECK_EQUAL( gcd<T>(  0u,   0u), static_cast<T>( 0u) );
    BOOST_CHECK_EQUAL( gcd<T>(  7u,   0u), static_cast<T>( 7u) );
    BOOST_CHECK_EQUAL( gcd<T>(  0u,   9u), static_cast<T>( 9u) );
    BOOST_CHECK_EQUAL( gcd<T>( 42u,  30u), static_cast<T>( 6u) );
    BOOST_CHECK_EQUAL( gcd<T>(  3u,   7u), static_cast<T>( 1u) );
    BOOST_CHECK_EQUAL( gcd<T>(  8u,   9u), static_cast<T>( 1u) );
    BOOST_CHECK_EQUAL( gcd<T>(  7u,  49u), static_cast<T>( 7u) );
}

// GCD at compile-time
BOOST_AUTO_TEST_CASE( gcd_static_test )
{
#ifndef BOOST_MSVC
    using boost::math::static_gcd;
#else
    using namespace boost::math;
#endif

    // Can't use "BOOST_CHECK_EQUAL", otherwise the "value" member will be
    // disqualified as compile-time-only constant, needing explicit definition
    BOOST_CHECK( (static_gcd< 1,  1>::value) == 1 );
    BOOST_CHECK( (static_gcd< 0,  0>::value) == 0 );
    BOOST_CHECK( (static_gcd< 7,  0>::value) == 7 );
    BOOST_CHECK( (static_gcd< 0,  9>::value) == 9 );
    BOOST_CHECK( (static_gcd<42, 30>::value) == 6 );
    BOOST_CHECK( (static_gcd< 3,  7>::value) == 1 );
    BOOST_CHECK( (static_gcd< 8,  9>::value) == 1 );
    BOOST_CHECK( (static_gcd< 7, 49>::value) == 7 );
}

// TODO: non-built-in signed and unsigned integer tests, with and without
// numeric_limits specialization; polynominal tests; note any changes if
// built-ins switch to binary-GCD algorithm

BOOST_AUTO_TEST_SUITE_END()


// LCM tests
BOOST_AUTO_TEST_SUITE( lcm_test_suite )

// LCM on built-in signed integer types
BOOST_AUTO_TEST_CASE_TEMPLATE( lcm_int_test, T, builtin_signed_test_types )
{
#ifndef BOOST_MSVC
    using boost::math::lcm;
#else
    using namespace boost::math;
#endif

    // Originally from Boost.Rational tests
    BOOST_CHECK_EQUAL( lcm<T>(  1,  -1), static_cast<T>( 1) );
    BOOST_CHECK_EQUAL( lcm<T>( -1,   1), static_cast<T>( 1) );
    BOOST_CHECK_EQUAL( lcm<T>(  1,   1), static_cast<T>( 1) );
    BOOST_CHECK_EQUAL( lcm<T>( -1,  -1), static_cast<T>( 1) );
    BOOST_CHECK_EQUAL( lcm<T>(  0,   0), static_cast<T>( 0) );
    BOOST_CHECK_EQUAL( lcm<T>(  6,   0), static_cast<T>( 0) );
    BOOST_CHECK_EQUAL( lcm<T>(  0,   7), static_cast<T>( 0) );
    BOOST_CHECK_EQUAL( lcm<T>( -5,   0), static_cast<T>( 0) );
    BOOST_CHECK_EQUAL( lcm<T>(  0,  -4), static_cast<T>( 0) );
    BOOST_CHECK_EQUAL( lcm<T>( 18,  30), static_cast<T>(90) );
    BOOST_CHECK_EQUAL( lcm<T>( -6,   9), static_cast<T>(18) );
    BOOST_CHECK_EQUAL( lcm<T>(-10, -10), static_cast<T>(10) );
    BOOST_CHECK_EQUAL( lcm<T>( 25, -10), static_cast<T>(50) );
    BOOST_CHECK_EQUAL( lcm<T>(  3,   7), static_cast<T>(21) );
    BOOST_CHECK_EQUAL( lcm<T>(  8,   9), static_cast<T>(72) );
    BOOST_CHECK_EQUAL( lcm<T>(  7,  49), static_cast<T>(49) );
}

// LCM on built-in unsigned integer types
BOOST_AUTO_TEST_CASE_TEMPLATE(lcm_unsigned_test, T, builtin_unsigned_test_types)
{
#ifndef BOOST_MSVC
    using boost::math::lcm;
#else
    using namespace boost::math;
#endif

    BOOST_CHECK_EQUAL( lcm<T>(  1u,   1u), static_cast<T>( 1u) );
    BOOST_CHECK_EQUAL( lcm<T>(  0u,   0u), static_cast<T>( 0u) );
    BOOST_CHECK_EQUAL( lcm<T>(  6u,   0u), static_cast<T>( 0u) );
    BOOST_CHECK_EQUAL( lcm<T>(  0u,   7u), static_cast<T>( 0u) );
    BOOST_CHECK_EQUAL( lcm<T>( 18u,  30u), static_cast<T>(90u) );
    BOOST_CHECK_EQUAL( lcm<T>(  3u,   7u), static_cast<T>(21u) );
    BOOST_CHECK_EQUAL( lcm<T>(  8u,   9u), static_cast<T>(72u) );
    BOOST_CHECK_EQUAL( lcm<T>(  7u,  49u), static_cast<T>(49u) );
}

// LCM at compile-time
BOOST_AUTO_TEST_CASE( lcm_static_test )
{
#ifndef BOOST_MSVC
    using boost::math::static_lcm;
#else
    using namespace boost::math;
#endif

    // Can't use "BOOST_CHECK_EQUAL", otherwise the "value" member will be
    // disqualified as compile-time-only constant, needing explicit definition
    BOOST_CHECK( (static_lcm< 1,  1>::value) ==  1 );
    BOOST_CHECK( (static_lcm< 0,  0>::value) ==  0 );
    BOOST_CHECK( (static_lcm< 6,  0>::value) ==  0 );
    BOOST_CHECK( (static_lcm< 0,  7>::value) ==  0 );
    BOOST_CHECK( (static_lcm<18, 30>::value) == 90 );
    BOOST_CHECK( (static_lcm< 3,  7>::value) == 21 );
    BOOST_CHECK( (static_lcm< 8,  9>::value) == 72 );
    BOOST_CHECK( (static_lcm< 7, 49>::value) == 49 );
}

// TODO: see GCD to-do

BOOST_AUTO_TEST_SUITE_END()
