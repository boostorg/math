//  Boost cayley_element_test.cpp test file  ---------------------------------//

//  Copyright 2005 Daryle Walker.  Use, modification, and distribution are
//  subject to the Boost Software License, Version 1.0.  (See accompanying file
//  LICENSE_1_0.txt or a copy at <http://www.boost.org/LICENSE_1_0.txt>.)

//  See <http://www.boost.org/libs/math/> for the library's home page.

//  Revision History
//   23 Jul 2005  Initial version (Daryle Walker)

#define BOOST_AUTO_TEST_MAIN  "cayley element test"

#include <boost/math/cayley_element.hpp>  // for boost:...:cayley_element, etc.

#include <boost/test/auto_unit_test.hpp>             // for BOOST_CHECK_EQUAL
#include <boost/test/floating_point_comparison.hpp>  // for BOOST_CHECK_CLOSE

#include <boost/rational.hpp>  // for boost::rational

#include <cmath>    // for std::asin
#include <cstddef>  // for std::size_t


// Types
using boost::math::cayley_element;
using boost::math::negated_cayley_element;
using boost::math::scaled_cayley_element;

typedef boost::rational<int>  ri_type;

typedef         cayley_element   ce;
typedef negated_cayley_element  nce;

typedef scaled_cayley_element<int>       ci_type;
typedef scaled_cayley_element<unsigned>  cu_type;
typedef scaled_cayley_element<double>    cd_type;
typedef scaled_cayley_element<ri_type>   cr_type;


// Testing preparations
BOOST_TEST_DONT_PRINT_LOG_VALUE( ce );
BOOST_TEST_DONT_PRINT_LOG_VALUE( nce );
BOOST_TEST_DONT_PRINT_LOG_VALUE( ci_type );
BOOST_TEST_DONT_PRINT_LOG_VALUE( cu_type );
BOOST_TEST_DONT_PRINT_LOG_VALUE( cd_type );
BOOST_TEST_DONT_PRINT_LOG_VALUE( cr_type );


// Element constructors, accessors, mutators, and Boolean conversion
BOOST_AUTO_UNIT_TEST( cayley_constructor_test )
{
    {
        ce  t1( 3 );
        BOOST_CHECK_EQUAL( 3u, t1.basis() );
        BOOST_CHECK( t1 );

        t1.basis( 5 );
        BOOST_CHECK_EQUAL( 5u, t1.basis() );
        BOOST_CHECK( t1 );
    }

    {
        nce  t2( 7 );
        BOOST_CHECK_EQUAL( 7u, t2.basis() );
        BOOST_CHECK( not t2.negative() );
        BOOST_CHECK( t2 );

        t2.basis( 11 );
        t2.negative( true );
        BOOST_CHECK_EQUAL( 11u, t2.basis() );
        BOOST_CHECK( t2.negative() );
        BOOST_CHECK( t2 );

        nce const  t3( 13, true );
        BOOST_CHECK_EQUAL( 13u, t3.basis() );
        BOOST_CHECK( t3.negative() );
        BOOST_CHECK( t3 );

        ce const   t4( 17 );
        nce const  t5( t4 );
        BOOST_CHECK_EQUAL( t4.basis(), t5.basis() );
        BOOST_CHECK( not t5.negative() );
        BOOST_CHECK( t4 && t5 );
    }

    {
        ci_type  t6( 19, 2 );
        BOOST_CHECK_EQUAL( 19u, t6.basis() );
        BOOST_CHECK_EQUAL( 2, t6.scale() );
        BOOST_CHECK( t6 );

        t6.basis( 23 );
        t6.scale( 0 );
        BOOST_CHECK_EQUAL( 23u, t6.basis() );
        BOOST_CHECK_EQUAL( 0, t6.scale() );
        BOOST_CHECK( not t6 );

        nce const      t7( 29, true );
        ci_type const  t8( t7 );
        BOOST_CHECK_EQUAL( t7.basis(), t8.basis() );
        BOOST_CHECK_EQUAL( -1, t8.scale() );
        BOOST_CHECK( t7 && t8 );

        ce const       t9( 31 );
        ci_type const  t10( t9 );
        BOOST_CHECK_EQUAL( t9.basis(), t10.basis() );
        BOOST_CHECK_EQUAL( 1, t10.scale() );
        BOOST_CHECK( t9 && t10 );

        ci_type const  t11( 3 );
        BOOST_CHECK_EQUAL( 0u, t11.basis() );
        BOOST_CHECK_EQUAL( 3, t11.scale() );
        BOOST_CHECK( t11 );

        ci_type const  t12;
        BOOST_CHECK_EQUAL( 0u, t12.basis() );
        BOOST_CHECK_EQUAL( 0, t12.scale() );
        BOOST_CHECK( not t12 );

        ci_type const  t13( nce(ce( 37 )) );
        BOOST_CHECK_EQUAL( 37u, t13.basis() );
        BOOST_CHECK_EQUAL( 1, t13.scale() );
        BOOST_CHECK( t13 );

        cd_type const  t14( ci_type(41u, -4) );
        BOOST_CHECK_EQUAL( 41u, t14.basis() );
        BOOST_CHECK_CLOSE( -4.0, t14.scale(), 0.1 );
        BOOST_CHECK( t14 );
    }
}

// Rung & index counting
BOOST_AUTO_UNIT_TEST( cayley_rung_index_test )
{
    BOOST_CHECK_EQUAL( 0u, ce::maximum_index_for_rung(0) );
    BOOST_CHECK_EQUAL( 0u, ce::minimum_rung_for_index(0) );

    BOOST_CHECK_EQUAL( 1u, ce::maximum_index_for_rung(1) );
    BOOST_CHECK_EQUAL( 1u, ce::minimum_rung_for_index(1) );

    BOOST_CHECK_EQUAL( 3u, ce::maximum_index_for_rung(2) );
    BOOST_CHECK_EQUAL( 2u, ce::minimum_rung_for_index(2) );
    BOOST_CHECK_EQUAL( 2u, ce::minimum_rung_for_index(3) );

    BOOST_CHECK_EQUAL( 7u, ce::maximum_index_for_rung(3) );
    BOOST_CHECK_EQUAL( 3u, ce::minimum_rung_for_index(4) );
    BOOST_CHECK_EQUAL( 3u, ce::minimum_rung_for_index(5) );
    BOOST_CHECK_EQUAL( 3u, ce::minimum_rung_for_index(6) );
    BOOST_CHECK_EQUAL( 3u, ce::minimum_rung_for_index(7) );

    BOOST_CHECK_EQUAL( 15u, ce::maximum_index_for_rung(4) );
    BOOST_CHECK_EQUAL( 4u, ce::minimum_rung_for_index(8) );
    BOOST_CHECK_EQUAL( 4u, ce::minimum_rung_for_index(9) );
    BOOST_CHECK_EQUAL( 4u, ce::minimum_rung_for_index(10) );
    BOOST_CHECK_EQUAL( 4u, ce::minimum_rung_for_index(11) );
    BOOST_CHECK_EQUAL( 4u, ce::minimum_rung_for_index(12) );
    BOOST_CHECK_EQUAL( 4u, ce::minimum_rung_for_index(13) );
    BOOST_CHECK_EQUAL( 4u, ce::minimum_rung_for_index(14) );
    BOOST_CHECK_EQUAL( 4u, ce::minimum_rung_for_index(15) );

    BOOST_CHECK_EQUAL( 31u, ce::maximum_index_for_rung(5) );
    BOOST_CHECK_EQUAL( 5u, ce(16).minimum_rung() );
    BOOST_CHECK_EQUAL( 5u, ce(17).minimum_rung() );
    BOOST_CHECK_EQUAL( 5u, ce(18).minimum_rung() );
    BOOST_CHECK_EQUAL( 5u, ce(19).minimum_rung() );
    BOOST_CHECK_EQUAL( 5u, ce(20).minimum_rung() );
    BOOST_CHECK_EQUAL( 5u, nce(21).minimum_rung() );
    BOOST_CHECK_EQUAL( 5u, nce(22).minimum_rung() );
    BOOST_CHECK_EQUAL( 5u, nce(23).minimum_rung() );
    BOOST_CHECK_EQUAL( 5u, nce(24).minimum_rung() );
    BOOST_CHECK_EQUAL( 5u, nce(25).minimum_rung() );
    BOOST_CHECK_EQUAL( 5u, cu_type(26, 1u).minimum_rung() );
    BOOST_CHECK_EQUAL( 5u, cu_type(27, 2u).minimum_rung() );
    BOOST_CHECK_EQUAL( 5u, cu_type(28, 3u).minimum_rung() );
    BOOST_CHECK_EQUAL( 5u, cu_type(29, 5u).minimum_rung() );
    BOOST_CHECK_EQUAL( 5u, cu_type(30, 7u).minimum_rung() );
    BOOST_CHECK_EQUAL( 5u, cu_type(31, 11u).minimum_rung() );
}

// Equality operators
BOOST_AUTO_UNIT_TEST( cayley_equal_op_test )
{
    ce const  t1( 0 ), t2( 1 ), t3( 0 );

    BOOST_CHECK_EQUAL( t1, t1 );
    BOOST_CHECK_EQUAL( t2, t2 );
    BOOST_CHECK_EQUAL( t3, t3 );

    BOOST_CHECK_EQUAL( t1, t3 );
    BOOST_CHECK_EQUAL( t3, t1 );

    BOOST_CHECK( !(t1 == t2) );
    BOOST_CHECK( !(t2 == t1) );
    BOOST_CHECK( !(t2 == t3) );
    BOOST_CHECK( !(t3 == t2) );

    BOOST_CHECK( !(t1 != t1) );
    BOOST_CHECK( !(t2 != t2) );
    BOOST_CHECK( !(t3 != t3) );
    BOOST_CHECK( !(t1 != t3) );
    BOOST_CHECK( !(t3 != t1) );
    BOOST_CHECK( t1 != t2 );
    BOOST_CHECK( t2 != t1 );
    BOOST_CHECK( t2 != t3 );
    BOOST_CHECK( t3 != t2 );

    nce const  t4( 2 ), t5( 3 ), t6( 2 ), t7( 2, true ), t8( 5, true );

    BOOST_CHECK_EQUAL( t4, t4 );  BOOST_CHECK_EQUAL( t5, t5 );
    BOOST_CHECK_EQUAL( t6, t6 );  BOOST_CHECK_EQUAL( t7, t7 );
    BOOST_CHECK_EQUAL( t8, t8 );

    BOOST_CHECK_EQUAL( t4, t6 );  BOOST_CHECK_EQUAL( t6, t4 );

    BOOST_CHECK( !(t4 == t5) );  BOOST_CHECK( !(t4 == t7) );
    BOOST_CHECK( !(t4 == t8) );
    BOOST_CHECK( !(t5 == t4) );  BOOST_CHECK( !(t5 == t6) );
    BOOST_CHECK( !(t5 == t7) );  BOOST_CHECK( !(t5 == t8) );
    BOOST_CHECK( !(t6 == t5) );  BOOST_CHECK( !(t6 == t7) );
    BOOST_CHECK( !(t6 == t8) );
    BOOST_CHECK( !(t7 == t4) );  BOOST_CHECK( !(t7 == t5) );
    BOOST_CHECK( !(t7 == t6) );  BOOST_CHECK( !(t7 == t8) );
    BOOST_CHECK( !(t8 == t4) );  BOOST_CHECK( !(t8 == t5) );
    BOOST_CHECK( !(t8 == t6) );  BOOST_CHECK( !(t8 == t7) );

    BOOST_CHECK( !(t4 != t4) );  BOOST_CHECK( !(t5 != t5) );
    BOOST_CHECK( !(t6 != t6) );  BOOST_CHECK( !(t7 != t7) );
    BOOST_CHECK( !(t8 != t8) );

    BOOST_CHECK( !(t4 != t6) );  BOOST_CHECK( !(t6 != t4) );

    BOOST_CHECK( t4 != t5 );  BOOST_CHECK( t4 != t7 );
    BOOST_CHECK( t4 != t8 );
    BOOST_CHECK( t5 != t4 );  BOOST_CHECK( t5 != t6 );
    BOOST_CHECK( t5 != t7 );  BOOST_CHECK( t5 != t8 );
    BOOST_CHECK( t6 != t5 );  BOOST_CHECK( t6 != t7 );
    BOOST_CHECK( t6 != t8 );
    BOOST_CHECK( t7 != t4 );  BOOST_CHECK( t7 != t5 );
    BOOST_CHECK( t7 != t6 );  BOOST_CHECK( t7 != t8 );
    BOOST_CHECK( t8 != t4 );  BOOST_CHECK( t8 != t5 );
    BOOST_CHECK( t8 != t6 );  BOOST_CHECK( t8 != t7 );

    cu_type const  t9, t10( 7, 0 ), t11( 11, 0 );
    cu_type const  t12( 1 ), t13( 11, 1 ), t14( 11, 1 );

    BOOST_CHECK_EQUAL( t9, t9 );    BOOST_CHECK_EQUAL( t10, t10 );
    BOOST_CHECK_EQUAL( t11, t11 );  BOOST_CHECK_EQUAL( t12, t12 );
    BOOST_CHECK_EQUAL( t13, t13 );  BOOST_CHECK_EQUAL( t14, t14 );

    BOOST_CHECK_EQUAL( t9, t10 );   BOOST_CHECK_EQUAL( t9, t11 );
    BOOST_CHECK_EQUAL( t10, t9 );   BOOST_CHECK_EQUAL( t11, t9 );
    BOOST_CHECK_EQUAL( t10, t11 );  BOOST_CHECK_EQUAL( t11, t10 );
    BOOST_CHECK_EQUAL( t13, t14 );  BOOST_CHECK_EQUAL( t14, t13 );

    BOOST_CHECK( !(t9 == t12) );   BOOST_CHECK( !(t9 == t13) );
    BOOST_CHECK( !(t9 == t14) );
    BOOST_CHECK( !(t10 == t12) );  BOOST_CHECK( !(t10 == t13) );
    BOOST_CHECK( !(t10 == t14) );
    BOOST_CHECK( !(t11 == t12) );  BOOST_CHECK( !(t11 == t13) );
    BOOST_CHECK( !(t11 == t14) );
    BOOST_CHECK( !(t12 == t9) );   BOOST_CHECK( !(t12 == t10) );
    BOOST_CHECK( !(t12 == t11) );  BOOST_CHECK( !(t12 == t13) );
    BOOST_CHECK( !(t12 == t14) );
    BOOST_CHECK( !(t13 == t9) );   BOOST_CHECK( !(t13 == t10) );
    BOOST_CHECK( !(t13 == t11) );  BOOST_CHECK( !(t13 == t12) );
    BOOST_CHECK( !(t14 == t9) );   BOOST_CHECK( !(t14 == t10) );
    BOOST_CHECK( !(t14 == t11) );  BOOST_CHECK( !(t14 == t12) );

    BOOST_CHECK( !(t9 != t9) );    BOOST_CHECK( !(t10 != t10) );
    BOOST_CHECK( !(t11 != t11) );  BOOST_CHECK( !(t12 != t12) );
    BOOST_CHECK( !(t13 != t13) );  BOOST_CHECK( !(t14 != t14) );

    BOOST_CHECK( !(t9 != t10) );   BOOST_CHECK( !(t9 != t11) );
    BOOST_CHECK( !(t10 != t9) );   BOOST_CHECK( !(t11 != t9) );
    BOOST_CHECK( !(t10 != t11) );  BOOST_CHECK( !(t11 != t10) );
    BOOST_CHECK( !(t13 != t14) );  BOOST_CHECK( !(t14 != t13) );

    BOOST_CHECK( t9 != t12 );   BOOST_CHECK( t9 != t13 );
    BOOST_CHECK( t9 != t14 );
    BOOST_CHECK( t10 != t12 );  BOOST_CHECK( t10 != t13 );
    BOOST_CHECK( t10 != t14 );
    BOOST_CHECK( t11 != t12 );  BOOST_CHECK( t11 != t13 );
    BOOST_CHECK( t11 != t14 );
    BOOST_CHECK( t12 != t9 );   BOOST_CHECK( t12 != t10 );
    BOOST_CHECK( t12 != t11 );  BOOST_CHECK( t12 != t13 );
    BOOST_CHECK( t12 != t14 );
    BOOST_CHECK( t13 != t9 );   BOOST_CHECK( t13 != t10 );
    BOOST_CHECK( t13 != t11 );  BOOST_CHECK( t13 != t12 );
    BOOST_CHECK( t14 != t9 );   BOOST_CHECK( t14 != t10 );
    BOOST_CHECK( t14 != t11 );  BOOST_CHECK( t14 != t12 );
}

// Unary operators
BOOST_AUTO_UNIT_TEST( cayley_unary_op_test )
{
    ce const  t1( 0 ), t2( 1 );

    BOOST_CHECK( !(!t1) );
    BOOST_CHECK( !(!t2) );

    nce const  t2a( t2 ), t3( 2 ), t4( 3, true );

    BOOST_CHECK( !(!t2a) );
    BOOST_CHECK( !(!t3) );
    BOOST_CHECK( !(!t4) );

    ci_type const  t2b( t2 ), t4a( t4 );
    cd_type const  t5( 5, 1.0 ), t6( 7, 0.0 );
    ci_type const  t7( 11, -2 );

    BOOST_CHECK( !(!t2b) );
    BOOST_CHECK( !(!t4a) );
    BOOST_CHECK( !(!t5) );
    BOOST_CHECK( !t6 );
    BOOST_CHECK( !(!t7) );

    BOOST_CHECK_EQUAL( +t1, t1 );
    BOOST_CHECK_EQUAL( +t2, t2 );

    BOOST_CHECK_EQUAL( +t2a, t2a );
    BOOST_CHECK_EQUAL( +t3, t3 );
    BOOST_CHECK_EQUAL( +t4, t4 );

    BOOST_CHECK_EQUAL( +t2b, t2b );
    BOOST_CHECK_EQUAL( +t4a, t4a );
    BOOST_CHECK_EQUAL( +t5, t5 );
    BOOST_CHECK_EQUAL( +t6, t6 );
    BOOST_CHECK_EQUAL( +t7, t7 );

    BOOST_CHECK_EQUAL( -t1, nce(0, true) );
    BOOST_CHECK_EQUAL( -t2, nce(1, true) );

    BOOST_CHECK_EQUAL( -t2a, nce(1, true) );
    BOOST_CHECK_EQUAL( -t3, nce(2, true) );
    BOOST_CHECK_EQUAL( -t4, nce(3) );

    BOOST_CHECK_EQUAL( -t2b, ci_type(1, -1) );
    BOOST_CHECK_EQUAL( -t4a, ci_type(3, 1) );
    BOOST_CHECK_EQUAL( -t5, cd_type(5, -1.0) );
    BOOST_CHECK_EQUAL( -t6, t6 );
    BOOST_CHECK_EQUAL( -t7, ci_type(11, 2) );
}

// Unary-operator-like member functions
BOOST_AUTO_UNIT_TEST( cayley_unary_self_test )
{
    ce const  t1( 0 );
    ce        t2 = t1;

    t2.same_self();
    BOOST_CHECK_EQUAL( +t1, t2 );

    nce const  t3( 1 );
    nce        t4 = t3;

    t4.same_self();
    BOOST_CHECK_EQUAL( +t3, t4 );
    t4.negate_self();
    BOOST_CHECK_EQUAL( -t3, t4 );

    ci_type const  t5( 2, 1 ), t6( 2, 0 );
    ci_type        t7 = t5, t8 = t6;

    t7.same_self();
    BOOST_CHECK_EQUAL( +t5, t7 );
    t7.negate_self();
    BOOST_CHECK_EQUAL( -t5, t7 );
    t7.not_self();
    BOOST_CHECK_EQUAL( t6, t7 );
    t8.same_self();
    BOOST_CHECK_EQUAL( +t6, t8 );
    t8.negate_self();
    BOOST_CHECK_EQUAL( -t6, t8 );
    t8.not_self();
    BOOST_CHECK_EQUAL( t5, t8 );
}

// Condition functions
BOOST_AUTO_UNIT_TEST( cayley_condition_function_test )
{
    using boost::math::dot_product;
    using boost::math::conj;
    using boost::math::abs;
    using boost::math::arg;
    using boost::math::norm;
    using boost::math::sgn;
    using boost::math::reciprocal;
    using boost::bad_rational;

    // dot product
    BOOST_CHECK_EQUAL( 1, dot_product(ce( 5 ), ce( 5 )) );
    BOOST_CHECK_EQUAL( 0, dot_product(ce( 2 ), ce( 3 )) );

    BOOST_CHECK_EQUAL(  1, dot_product(nce( 5 ), nce( 5 )) );
    BOOST_CHECK_EQUAL(  0, dot_product(nce( 2 ), nce( 3 )) );
    BOOST_CHECK_EQUAL( -1, dot_product(nce( 5, true ), nce( 5 )) );
    BOOST_CHECK_EQUAL( -1, dot_product(nce( 5 ), nce( 5, true )) );
    BOOST_CHECK_EQUAL(  1, dot_product(nce( 5, true ), nce( 5, true )) );
    BOOST_CHECK_EQUAL(  0, dot_product(nce( 2, true ), nce( 3 )) );
    BOOST_CHECK_EQUAL(  0, dot_product(nce( 2 ), nce( 3, true )) );
    BOOST_CHECK_EQUAL(  0, dot_product(nce( 2, true ), nce( 3, true )) );

    BOOST_CHECK_EQUAL(  1, dot_product(ci_type( 5, 1 ), ci_type( 5, 1 )) );
    BOOST_CHECK_EQUAL(  0, dot_product(ci_type( 2, 1 ), ci_type( 3, 1 )) );
    BOOST_CHECK_EQUAL( -1, dot_product(ci_type( 5, -1 ), ci_type( 5, 1 )) );
    BOOST_CHECK_EQUAL( -1, dot_product(ci_type( 5, 1 ), ci_type( 5, -1 )) );
    BOOST_CHECK_EQUAL(  1, dot_product(ci_type( 5, -1 ), ci_type( 5, -1 )) );
    BOOST_CHECK_EQUAL(  0, dot_product(ci_type( 2, -1 ), ci_type( 3, 1 )) );
    BOOST_CHECK_EQUAL(  0, dot_product(ci_type( 2, 1 ), ci_type( 3, -1 )) );
    BOOST_CHECK_EQUAL(  0, dot_product(ci_type( 2, -1 ), ci_type( 3, -1 )) );
    BOOST_CHECK_EQUAL(  0, dot_product(ci_type( 7, 0 ), ci_type( 7, 4 )) );
    BOOST_CHECK_EQUAL( -8, dot_product(ci_type( 7, -2 ), ci_type( 7, 4 )) );
    BOOST_CHECK_EQUAL(  6, dot_product(ci_type( 11, 3 ), ci_type( 11, 2 )) );
    BOOST_CHECK_EQUAL( 10, dot_product(ci_type( 13, -2 ), ci_type( 13, -5 )) );

    // conjugate
    ce const  t1( 0 ), t2( 1 );

    BOOST_CHECK_EQUAL( conj(t1), t1 );
    BOOST_CHECK_EQUAL( conj(t2), -t2 );

    nce const  t1a( t1 ), t2a( t2 ), t3( 2, true ), t4( 0, true );

    BOOST_CHECK_EQUAL( conj(t1a), t1a );
    BOOST_CHECK_EQUAL( conj(t2a), -t2a );
    BOOST_CHECK_EQUAL( conj(t3), -t3 );
    BOOST_CHECK_EQUAL( conj(t4), t4 );

    ci_type const  t1b( t1 ), t2b( t2 ), t3a( t3 ), t4a( t4 );
    ci_type const  t5( 0, 6 ), t6( 0, -8 ), t7( 3, 9 ), t8( 5, -12 );
    ci_type const  t9( 0, 0 ), t10( 7, 0 );

    BOOST_CHECK_EQUAL( conj(t1b), t1b );
    BOOST_CHECK_EQUAL( conj(t2b), -t2b );
    BOOST_CHECK_EQUAL( conj(t3a), -t3a );
    BOOST_CHECK_EQUAL( conj(t4a), t4a );
    BOOST_CHECK_EQUAL( conj(t5), t5 );
    BOOST_CHECK_EQUAL( conj(t6), t6 );
    BOOST_CHECK_EQUAL( conj(t7), -t7 );
    BOOST_CHECK_EQUAL( conj(t8), -t8 );
    BOOST_CHECK_EQUAL( conj(t9), t9 );
    BOOST_CHECK_EQUAL( conj(t10), -t10 );

    // absolute value (2-norm)
    BOOST_CHECK_EQUAL( abs(t1), 1 );
    BOOST_CHECK_EQUAL( abs(t2), 1 );

    BOOST_CHECK_EQUAL( abs(t1a), 1 );
    BOOST_CHECK_EQUAL( abs(t2a), 1 );
    BOOST_CHECK_EQUAL( abs(t3), 1 );
    BOOST_CHECK_EQUAL( abs(t4), 1 );

    BOOST_CHECK_EQUAL( abs(t1b), 1 );
    BOOST_CHECK_EQUAL( abs(t2b), 1 );
    BOOST_CHECK_EQUAL( abs(t3a), 1 );
    BOOST_CHECK_EQUAL( abs(t4a), 1 );
    BOOST_CHECK_EQUAL( abs(t5), 6 );
    BOOST_CHECK_EQUAL( abs(t6), 8 );
    BOOST_CHECK_EQUAL( abs(t7), 9 );
    BOOST_CHECK_EQUAL( abs(t8), 12 );
    BOOST_CHECK_EQUAL( abs(t9), 0 );
    BOOST_CHECK_EQUAL( abs(t10), 0 );

    // argument angle
    double const  pi_2 = std::asin( 1.0 ), pi = 2.0 * pi_2;

    BOOST_CHECK_CLOSE( arg(t1), 0.0, 0.1 );
    BOOST_CHECK_CLOSE( arg(t2), pi_2, 0.1 );

    BOOST_CHECK_CLOSE( arg(t1a), 0.0, 0.1 );
    BOOST_CHECK_CLOSE( arg(t2a), pi_2, 0.1 );
    BOOST_CHECK_CLOSE( arg(t3), pi_2, 0.1 );
    BOOST_CHECK_CLOSE( arg(t4), pi, 0.1 );

    BOOST_CHECK_CLOSE( arg(cd_type( t1b )), 0.0, 0.1 );
    BOOST_CHECK_CLOSE( arg(cd_type( t2b )), pi_2, 0.1 );
    BOOST_CHECK_CLOSE( arg(cd_type( t3a )), pi_2, 0.1 );
    BOOST_CHECK_CLOSE( arg(cd_type( t4a )), pi, 0.1 );
    BOOST_CHECK_CLOSE( arg(cd_type( t5 )), 0.0, 0.1 );
    BOOST_CHECK_CLOSE( arg(cd_type( t6 )), pi, 0.1 );
    BOOST_CHECK_CLOSE( arg(cd_type( t7 )), pi_2, 0.1 );
    BOOST_CHECK_CLOSE( arg(cd_type( t8 )), pi_2, 0.1 );
    BOOST_CHECK_CLOSE( arg(cd_type( t9 )), 0.0, 0.1 );
    BOOST_CHECK_CLOSE( arg(cd_type( t10 )), 0.0, 0.1 );

    // norm (quadratic, Cayley)
    BOOST_CHECK_EQUAL( norm(t1), 1 );
    BOOST_CHECK_EQUAL( norm(t2), 1 );

    BOOST_CHECK_EQUAL( norm(t1a), 1 );
    BOOST_CHECK_EQUAL( norm(t2a), 1 );
    BOOST_CHECK_EQUAL( norm(t3), 1 );
    BOOST_CHECK_EQUAL( norm(t4), 1 );

    BOOST_CHECK_EQUAL( norm(t1b), 1 );
    BOOST_CHECK_EQUAL( norm(t2b), 1 );
    BOOST_CHECK_EQUAL( norm(t3a), 1 );
    BOOST_CHECK_EQUAL( norm(t4a), 1 );
    BOOST_CHECK_EQUAL( norm(t5), 36 );
    BOOST_CHECK_EQUAL( norm(t6), 64 );
    BOOST_CHECK_EQUAL( norm(t7), 81 );
    BOOST_CHECK_EQUAL( norm(t8), 144 );
    BOOST_CHECK_EQUAL( norm(t9), 0 );
    BOOST_CHECK_EQUAL( norm(t10), 0 );

    // sign
    BOOST_CHECK_EQUAL( sgn(t1), t1 );
    BOOST_CHECK_EQUAL( sgn(t2), t2 );

    BOOST_CHECK_EQUAL( sgn(t1a), t1a );
    BOOST_CHECK_EQUAL( sgn(t2a), t2a );
    BOOST_CHECK_EQUAL( sgn(t3), t3 );
    BOOST_CHECK_EQUAL( sgn(t4), t4 );

    BOOST_CHECK_EQUAL( sgn(t1b), t1b );
    BOOST_CHECK_EQUAL( sgn(t2b), t2b );
    BOOST_CHECK_EQUAL( sgn(t3a), t3a );
    BOOST_CHECK_EQUAL( sgn(t4a), t4a );
    BOOST_CHECK_EQUAL( sgn(t5), ci_type(0, 1) );
    BOOST_CHECK_EQUAL( sgn(t6), ci_type(0, -1) );
    BOOST_CHECK_EQUAL( sgn(t7), ci_type(3, 1) );
    BOOST_CHECK_EQUAL( sgn(t8), ci_type(5, -1) );
    BOOST_CHECK_EQUAL( sgn(t9), t9 );
    BOOST_CHECK_EQUAL( sgn(t10), t10 );

    // reciprocal
    BOOST_CHECK_EQUAL( reciprocal(t1), t1 );
    BOOST_CHECK_EQUAL( reciprocal(t2), -t2 );

    BOOST_CHECK_EQUAL( reciprocal(t1a), t1a );
    BOOST_CHECK_EQUAL( reciprocal(t2a), -t2a );
    BOOST_CHECK_EQUAL( reciprocal(t3), -t3 );
    BOOST_CHECK_EQUAL( reciprocal(t4), t4 );

    BOOST_CHECK_EQUAL( reciprocal(cr_type( t1b )), cr_type(t1b) );
    BOOST_CHECK_EQUAL( reciprocal(cr_type( t2b )), cr_type(-t2b) );
    BOOST_CHECK_EQUAL( reciprocal(cr_type( t3a )), cr_type(-t3a) );
    BOOST_CHECK_EQUAL( reciprocal(cr_type( t4a )), cr_type(t4a) );
    BOOST_CHECK_EQUAL( reciprocal(cr_type( t5 )), cr_type(0, ri_type( 1, 6 )) );
    BOOST_CHECK_EQUAL( reciprocal(cr_type( t6 )), cr_type(0, ri_type( -1,
     8 )) );
    BOOST_CHECK_EQUAL( reciprocal(cr_type( t7 )), cr_type(3, ri_type( -1,
     9 )) );
    BOOST_CHECK_EQUAL( reciprocal(cr_type( t8 )), cr_type(5, ri_type( 1,
     12 )) );
    BOOST_CHECK_THROW( reciprocal(cr_type( t9 )), bad_rational );
    BOOST_CHECK_THROW( reciprocal(cr_type( t10 )), bad_rational );

    // other norms (supremum [inf-norm] and 1-norm)
    BOOST_CHECK_EQUAL( sup(t1), 1 );  BOOST_CHECK_EQUAL( ll(t1), 1 );
    BOOST_CHECK_EQUAL( sup(t2), 1 );  BOOST_CHECK_EQUAL( ll(t2), 1 );

    BOOST_CHECK_EQUAL( sup(t1a), 1 );  BOOST_CHECK_EQUAL( ll(t1a), 1 );
    BOOST_CHECK_EQUAL( sup(t2a), 1 );  BOOST_CHECK_EQUAL( ll(t2a), 1 );
    BOOST_CHECK_EQUAL( sup(t3), 1 );   BOOST_CHECK_EQUAL( ll(t3), 1 );
    BOOST_CHECK_EQUAL( sup(t4), 1 );   BOOST_CHECK_EQUAL( ll(t4), 1 );

    BOOST_CHECK_EQUAL( sup(t1b), 1 );  BOOST_CHECK_EQUAL( ll(t1b), 1 );
    BOOST_CHECK_EQUAL( sup(t2b), 1 );  BOOST_CHECK_EQUAL( ll(t2b), 1 );
    BOOST_CHECK_EQUAL( sup(t3a), 1 );  BOOST_CHECK_EQUAL( ll(t3a), 1 );
    BOOST_CHECK_EQUAL( sup(t4a), 1 );  BOOST_CHECK_EQUAL( ll(t4a), 1 );
    BOOST_CHECK_EQUAL( sup(t5), 6 );   BOOST_CHECK_EQUAL( ll(t5), 6 );
    BOOST_CHECK_EQUAL( sup(t6), 8 );   BOOST_CHECK_EQUAL( ll(t6), 8 );
    BOOST_CHECK_EQUAL( sup(t7), 9 );   BOOST_CHECK_EQUAL( ll(t7), 9 );
    BOOST_CHECK_EQUAL( sup(t8), 12 );  BOOST_CHECK_EQUAL( ll(t8), 12 );
    BOOST_CHECK_EQUAL( sup(t9), 0 );   BOOST_CHECK_EQUAL( ll(t9), 0 );
    BOOST_CHECK_EQUAL( sup(t10), 0 );  BOOST_CHECK_EQUAL( ll(t10), 0 );
}

// Self-conditioning member functions
BOOST_AUTO_UNIT_TEST( cayley_condition_self_test )
{
    using boost::math::conj;
    using boost::math::sgn;
    using boost::math::reciprocal;
    using boost::bad_rational;

    ce const  c1( 0 ), c2( 1 );
    ce        t1 = c1, t2 = c2;

    t1.sign_self();   BOOST_CHECK_EQUAL( sgn(c1), t1 );
    t2.sign_self();   BOOST_CHECK_EQUAL( sgn(c2), t2 );

    nce const  c3( c1 ), c4( 0, true ), c5( c2 ), c6( 1, true );
    nce        t3 = c3, t4 = c4, t5 = c5, t6 = c6;

    t3.sign_self();   BOOST_CHECK_EQUAL( sgn(c3), t3 );
    t4.sign_self();   BOOST_CHECK_EQUAL( sgn(c4), t4 );
    t5.sign_self();   BOOST_CHECK_EQUAL( sgn(c5), t5 );
    t6.sign_self();   BOOST_CHECK_EQUAL( sgn(c6), t6 );

    t3 = c3;  t3.conjugate_self();   BOOST_CHECK_EQUAL( conj(c3), t3 );
    t4 = c4;  t4.conjugate_self();   BOOST_CHECK_EQUAL( conj(c4), t4 );
    t5 = c5;  t5.conjugate_self();   BOOST_CHECK_EQUAL( conj(c5), t5 );
    t6 = c6;  t6.conjugate_self();   BOOST_CHECK_EQUAL( conj(c6), t6 );

    t3 = c3;  t3.reciprocate_self();   BOOST_CHECK_EQUAL( reciprocal(c3), t3 );
    t4 = c4;  t4.reciprocate_self();   BOOST_CHECK_EQUAL( reciprocal(c4), t4 );
    t5 = c5;  t5.reciprocate_self();   BOOST_CHECK_EQUAL( reciprocal(c5), t5 );
    t6 = c6;  t6.reciprocate_self();   BOOST_CHECK_EQUAL( reciprocal(c6), t6 );

    ci_type const  c7( c3 ), c8( c4 ), c9( c5 ), c10( c6 );
    ci_type        t7 = c7, t8 = c8, t9 = c9, t10 = c10;

    t7.reciprocate_self();  BOOST_CHECK_EQUAL( reciprocal(c7), t7 );
    t8.reciprocate_self();  BOOST_CHECK_EQUAL( reciprocal(c8), t8 );
    t9.reciprocate_self();  BOOST_CHECK_EQUAL( reciprocal(c9), t9 );
    t10.reciprocate_self();  BOOST_CHECK_EQUAL( reciprocal(c10), t10 );

    t7 = c7;  t7.sign_self();  BOOST_CHECK_EQUAL( sgn(c7), t7 );
    t8 = c8;  t8.sign_self();  BOOST_CHECK_EQUAL( sgn(c8), t8 );
    t9 = c9;  t9.sign_self();  BOOST_CHECK_EQUAL( sgn(c9), t9 );
    t10 = c10;  t10.sign_self();  BOOST_CHECK_EQUAL( sgn(c10), t10 );

    t7 = c7;  t7.conjugate_self();  BOOST_CHECK_EQUAL( conj(c7), t7 );
    t8 = c8;  t8.conjugate_self();  BOOST_CHECK_EQUAL( conj(c8), t8 );
    t9 = c9;  t9.conjugate_self();  BOOST_CHECK_EQUAL( conj(c9), t9 );
    t10 = c10;  t10.conjugate_self();  BOOST_CHECK_EQUAL( conj(c10), t10 );

    cr_type const  c11( 0, ri_type(4) ), c12( 0, ri_type(-3) ), c13( 0, 0 );
    cr_type const  c14( 2, ri_type(-5) ), c15( 2, ri_type(7) ), c16( 2, 0 );
    cr_type        t11 = c11, t12 = c12, t13 = c13;
    cr_type        t14 = c14, t15 = c15, t16 = c16;

    t11.reciprocate_self();  BOOST_CHECK_EQUAL( reciprocal(c11), t11 );
    t12.reciprocate_self();  BOOST_CHECK_EQUAL( reciprocal(c12), t12 );
    BOOST_CHECK_THROW( t13.reciprocate_self(), bad_rational );
    t14.reciprocate_self();  BOOST_CHECK_EQUAL( reciprocal(c14), t14 );
    t15.reciprocate_self();  BOOST_CHECK_EQUAL( reciprocal(c15), t15 );
    BOOST_CHECK_THROW( t16.reciprocate_self(), bad_rational );

    t11 = c11;  t11.conjugate_self();  BOOST_CHECK_EQUAL( conj(c11), t11 );
    t12 = c12;  t12.conjugate_self();  BOOST_CHECK_EQUAL( conj(c12), t12 );
    t13 = c13;  t13.conjugate_self();  BOOST_CHECK_EQUAL( conj(c13), t13 );
    t14 = c14;  t14.conjugate_self();  BOOST_CHECK_EQUAL( conj(c14), t14 );
    t15 = c15;  t15.conjugate_self();  BOOST_CHECK_EQUAL( conj(c15), t15 );
    t16 = c16;  t16.conjugate_self();  BOOST_CHECK_EQUAL( conj(c16), t16 );

    t11 = c11;  t11.sign_self();  BOOST_CHECK_EQUAL( sgn(c11), t11 );
    t12 = c12;  t12.sign_self();  BOOST_CHECK_EQUAL( sgn(c12), t12 );
    t13 = c13;  t13.sign_self();  BOOST_CHECK_EQUAL( sgn(c13), t13 );
    t14 = c14;  t14.sign_self();  BOOST_CHECK_EQUAL( sgn(c14), t14 );
    t15 = c15;  t15.sign_self();  BOOST_CHECK_EQUAL( sgn(c15), t15 );
    t16 = c16;  t16.sign_self();  BOOST_CHECK_EQUAL( sgn(c16), t16 );
}

// Scalar multiplication and division
BOOST_AUTO_UNIT_TEST( cayley_scalar_muldiv_op_test )
{
    using boost::bad_rational;

    BOOST_CHECK_EQUAL( 6 * ci_type(0, 7), ci_type(0, 42) );
    BOOST_CHECK_EQUAL( ci_type(1, 5) * -3, ci_type(1, -15) );

    BOOST_CHECK_EQUAL( ri_type(4) / cr_type(0, 3), cr_type(0, ri_type(4, 3)) );
    BOOST_CHECK_EQUAL( ri_type(12) / cr_type(2, 3), cr_type(2, ri_type(-4)) );
    BOOST_CHECK_THROW( ri_type(4) / cr_type(0, 0), bad_rational );
    BOOST_CHECK_THROW( ri_type(12) / cr_type(2, 0), bad_rational );
    BOOST_CHECK_EQUAL( cr_type(0, 5) / ri_type(6), cr_type(0, ri_type(5, 6)) );
    BOOST_CHECK_EQUAL( cr_type(3, 7) / ri_type(6), cr_type(3, ri_type(7, 6)) );
    BOOST_CHECK_THROW( cr_type(0, 5) / ri_type(0), bad_rational );
    BOOST_CHECK_THROW( cr_type(3, 7) / ri_type(0), bad_rational );

    BOOST_CHECK_EQUAL( 2 * nce(0), ci_type(0, 2) );
    BOOST_CHECK_EQUAL( nce(0) * 2, ci_type(0, 2) );
    BOOST_CHECK_EQUAL( 2 * nce(5), ci_type(5, 2) );
    BOOST_CHECK_EQUAL( nce(5) * 2, ci_type(5, 2) );
    BOOST_CHECK_EQUAL( 2 * nce(0, true), ci_type(0, -2) );
    BOOST_CHECK_EQUAL( nce(0, true) * 2, ci_type(0, -2) );
    BOOST_CHECK_EQUAL( 2 * nce(5, true), ci_type(5, -2) );
    BOOST_CHECK_EQUAL( nce(5, true) * 2, ci_type(5, -2) );

    BOOST_CHECK_EQUAL( ri_type(5) / nce(0), cr_type(0, ri_type(5)) );
    BOOST_CHECK_EQUAL( ri_type(5) / nce(0, true), cr_type(0, ri_type(-5)) );
    BOOST_CHECK_EQUAL( nce(0) / ri_type(5), cr_type(0, ri_type(1, 5)) );
    BOOST_CHECK_EQUAL( nce(0, true) / ri_type(5), cr_type(0, ri_type(-1, 5)) );
    BOOST_CHECK_THROW( nce(0) / ri_type(0), bad_rational );
    BOOST_CHECK_THROW( nce(0, true) / ri_type(0), bad_rational );
    BOOST_CHECK_EQUAL( ri_type(8) / nce(7), cr_type(7, ri_type(-8)) );
    BOOST_CHECK_EQUAL( ri_type(8) / nce(7, true), cr_type(7, ri_type(8)) );
    BOOST_CHECK_EQUAL( nce(7) / ri_type(8), cr_type(7, ri_type(1, 8)) );
    BOOST_CHECK_EQUAL( nce(7, true) / ri_type(8), cr_type(7, ri_type(-1, 8)) );
    BOOST_CHECK_THROW( nce(7) / ri_type(0), bad_rational );
    BOOST_CHECK_THROW( nce(7, true) / ri_type(0), bad_rational );

    BOOST_CHECK_EQUAL( ri_type(7) / ce(0), cr_type(0, ri_type(7)) );
    BOOST_CHECK_EQUAL( ce(0) / ri_type(7), cr_type(0, ri_type(1, 7)) );
    BOOST_CHECK_THROW( ce(0) / ri_type(0), bad_rational );
    BOOST_CHECK_EQUAL( ri_type(9) / ce(11), cr_type(11, ri_type(-9)) );
    BOOST_CHECK_EQUAL( ce(11) / ri_type(9), cr_type(11, ri_type(1, 9)) );
    BOOST_CHECK_THROW( ce(11) / ri_type(0), bad_rational );
}

// Element multiplication and division
BOOST_AUTO_UNIT_TEST( cayley_element_muldiv_op_test )
{
    using std::size_t;
    using boost::bad_rational;

    // basic elements
    size_t const  top_rung = 4u, top_index = ( 1u << top_rung );
    nce const     results[ top_index ][ top_index ] =
        {
            { +ce( 0), +ce( 1), +ce( 2), +ce( 3), +ce( 4), +ce( 5), +ce( 6),
               +ce( 7), +ce( 8), +ce( 9), +ce(10), +ce(11), +ce(12), +ce(13),
               +ce(14), +ce(15) },
            { +ce( 1), -ce( 0), +ce( 3), -ce( 2), +ce( 5), -ce( 4), -ce( 7),
               +ce( 6), +ce( 9), -ce( 8), -ce(11), +ce(10), -ce(13), +ce(12),
               +ce(15), -ce(14) },
            { +ce( 2), -ce( 3), -ce( 0), +ce( 1), +ce( 6), +ce( 7), -ce( 4),
               -ce( 5), +ce(10), +ce(11), -ce( 8), -ce( 9), -ce(14), -ce(15),
               +ce(12), +ce(13) },
            { +ce( 3), +ce( 2), -ce( 1), -ce( 0), +ce( 7), -ce( 6), +ce( 5),
               -ce( 4), +ce(11), -ce(10), +ce( 9), -ce( 8), -ce(15), +ce(14),
               -ce(13), +ce(12) },
            { +ce( 4), -ce( 5), -ce( 6), -ce( 7), -ce( 0), +ce( 1), +ce( 2),
               +ce( 3), +ce(12), +ce(13), +ce(14), +ce(15), -ce( 8), -ce( 9),
               -ce(10), -ce(11) },
            { +ce( 5), +ce( 4), -ce( 7), +ce( 6), -ce( 1), -ce( 0), -ce( 3),
               +ce( 2), +ce(13), -ce(12), +ce(15), -ce(14), +ce( 9), -ce( 8),
               +ce(11), -ce(10) },
            { +ce( 6), +ce( 7), +ce( 4), -ce( 5), -ce( 2), +ce( 3), -ce( 0),
               -ce( 1), +ce(14), -ce(15), -ce(12), +ce(13), +ce(10), -ce(11),
               -ce( 8), +ce( 9) },
            { +ce( 7), -ce( 6), +ce( 5), +ce( 4), -ce( 3), -ce( 2), +ce( 1),
               -ce( 0), +ce(15), +ce(14), -ce(13), -ce(12), +ce(11), +ce(10),
               -ce( 9), -ce( 8) },
            { +ce( 8), -ce( 9), -ce(10), -ce(11), -ce(12), -ce(13), -ce(14),
               -ce(15), -ce( 0), +ce( 1), +ce( 2), +ce( 3), +ce( 4), +ce( 5),
               +ce( 6), +ce( 7) },
            { +ce( 9), +ce( 8), -ce(11), +ce(10), -ce(13), +ce(12), +ce(15),
               -ce(14), -ce( 1), -ce( 0), -ce( 3), +ce( 2), -ce( 5), +ce( 4),
               +ce( 7), -ce( 6) },
            { +ce(10), +ce(11), +ce( 8), -ce( 9), -ce(14), -ce(15), +ce(12),
               +ce(13), -ce( 2), +ce( 3), -ce( 0), -ce( 1), -ce( 6), -ce( 7),
               +ce( 4), +ce( 5) },
            { +ce(11), -ce(10), +ce( 9), +ce( 8), -ce(15), +ce(14), -ce(13),
               +ce(12), -ce( 3), -ce( 2), +ce( 1), -ce( 0), -ce( 7), +ce( 6),
               -ce( 5), +ce( 4) },
            { +ce(12), +ce(13), +ce(14), +ce(15), +ce( 8), -ce( 9), -ce(10),
               -ce(11), -ce( 4), +ce( 5), +ce( 6), +ce( 7), -ce( 0), -ce( 1),
               -ce( 2), -ce( 3) },
            { +ce(13), -ce(12), +ce(15), -ce(14), +ce( 9), +ce( 8), +ce(11),
               -ce(10), -ce( 5), -ce( 4), +ce( 7), -ce( 6), +ce( 1), -ce( 0),
               +ce( 3), -ce( 2) },
            { +ce(14), -ce(15), -ce(12), +ce(13), +ce(10), -ce(11), +ce( 8),
               +ce( 9), -ce( 6), -ce( 7), -ce( 4), +ce( 5), +ce( 2), -ce( 3),
               -ce( 0), +ce( 1) },
            { +ce(15), +ce(14), -ce(13), -ce(12), +ce(11), +ce(10), -ce( 9),
               +ce( 8), -ce( 7), +ce( 6), -ce( 5), -ce( 4), +ce( 3), +ce( 2),
               -ce( 1), -ce( 0) }
        };

    for ( size_t i = 0 ; i < top_index ; ++i )
        for ( size_t j = 0 ; j < top_index ; ++j )
            BOOST_CHECK_MESSAGE( ce(i) * ce(j) == results[i][j], "e_" << i
             << " * e_" << j << " -> " << ((ce(i) * ce(j)).negative() ? '-'
             : '+') << "e_" << (ce(i) * ce(j)).basis() << " != "
             << (results[i][j].negative() ? '-' : '+') << "e_"
             << results[i][j].basis() << '.' );

    // signed elements
    BOOST_CHECK_EQUAL( +nce(0) * +nce(0), +nce(0) );
    BOOST_CHECK_EQUAL( -nce(0) * -nce(0), +nce(0) );
    BOOST_CHECK_EQUAL( -nce(0) * +nce(0), -nce(0) );
    BOOST_CHECK_EQUAL( +nce(0) * -nce(0), -nce(0) );

    BOOST_CHECK_EQUAL( +nce(1) * +nce(1), -nce(0) );
    BOOST_CHECK_EQUAL( -nce(1) * -nce(1), -nce(0) );
    BOOST_CHECK_EQUAL( -nce(1) * +nce(1), +nce(0) );
    BOOST_CHECK_EQUAL( +nce(1) * -nce(1), +nce(0) );

    BOOST_CHECK_EQUAL( +nce(0) * +nce(1), +nce(1) );
    BOOST_CHECK_EQUAL( -nce(0) * -nce(1), +nce(1) );
    BOOST_CHECK_EQUAL( -nce(0) * +nce(1), -nce(1) );
    BOOST_CHECK_EQUAL( +nce(0) * -nce(1), -nce(1) );
    BOOST_CHECK_EQUAL( +nce(1) * +nce(0), +nce(1) );
    BOOST_CHECK_EQUAL( -nce(1) * -nce(0), +nce(1) );
    BOOST_CHECK_EQUAL( -nce(1) * +nce(0), -nce(1) );
    BOOST_CHECK_EQUAL( +nce(1) * -nce(0), -nce(1) );

    BOOST_CHECK_EQUAL( +nce(2) * +nce(5), +nce(7) );
    BOOST_CHECK_EQUAL( -nce(2) * -nce(5), +nce(7) );
    BOOST_CHECK_EQUAL( -nce(2) * +nce(5), -nce(7) );
    BOOST_CHECK_EQUAL( +nce(2) * -nce(5), -nce(7) );
    BOOST_CHECK_EQUAL( +nce(5) * +nce(2), -nce(7) );
    BOOST_CHECK_EQUAL( -nce(5) * -nce(2), -nce(7) );
    BOOST_CHECK_EQUAL( -nce(5) * +nce(2), +nce(7) );
    BOOST_CHECK_EQUAL( +nce(5) * -nce(2), +nce(7) );

    // scaled elements
    BOOST_CHECK_EQUAL( ci_type(0,  2) * ci_type(0,  3),   6 * ce(0) );
    BOOST_CHECK_EQUAL( ci_type(0, -7) * ci_type(0, -5),  35 * ce(0) );
    BOOST_CHECK_EQUAL( ci_type(0, -5) * ci_type(0,  3), -15 * ce(0) );
    BOOST_CHECK_EQUAL( ci_type(0, 11) * ci_type(0, -2), -22 * ce(0) );

    BOOST_CHECK_EQUAL( ci_type(1,   3) * ci_type(1,  7), -21 * ce(0) );
    BOOST_CHECK_EQUAL( ci_type(1,  -2) * ci_type(1, -5), -10 * ce(0) );
    BOOST_CHECK_EQUAL( ci_type(1, -11) * ci_type(1, 13), 143 * ce(0) );
    BOOST_CHECK_EQUAL( ci_type(1,   7) * ci_type(1, -3),  21 * ce(0) );

    BOOST_CHECK_EQUAL( ci_type(0,   2) * ci_type(1,  19),  38 * ce(1) );
    BOOST_CHECK_EQUAL( ci_type(0,  -3) * ci_type(1, -17),  51 * ce(1) );
    BOOST_CHECK_EQUAL( ci_type(0,  -5) * ci_type(1,  13), -65 * ce(1) );
    BOOST_CHECK_EQUAL( ci_type(0,   7) * ci_type(1, -11), -77 * ce(1) );
    BOOST_CHECK_EQUAL( ci_type(1,  11) * ci_type(0,   7),  77 * ce(1) );
    BOOST_CHECK_EQUAL( ci_type(1, -13) * ci_type(0,  -5),  65 * ce(1) );
    BOOST_CHECK_EQUAL( ci_type(1, -17) * ci_type(0,   3), -51 * ce(1) );
    BOOST_CHECK_EQUAL( ci_type(1,  19) * ci_type(0,  -2), -38 * ce(1) );

    BOOST_CHECK_EQUAL( ci_type(2,  23) * ci_type(5,  53),  1219 * ce(7) );
    BOOST_CHECK_EQUAL( ci_type(2, -29) * ci_type(5, -47),  1363 * ce(7) );
    BOOST_CHECK_EQUAL( ci_type(2, -31) * ci_type(5,  43), -1333 * ce(7) );
    BOOST_CHECK_EQUAL( ci_type(2,  37) * ci_type(5, -41), -1517 * ce(7) );
    BOOST_CHECK_EQUAL( ci_type(5,  41) * ci_type(2,  37), -1517 * ce(7) );
    BOOST_CHECK_EQUAL( ci_type(5, -43) * ci_type(2, -31), -1333 * ce(7) );
    BOOST_CHECK_EQUAL( ci_type(5, -47) * ci_type(2,  29),  1363 * ce(7) );
    BOOST_CHECK_EQUAL( ci_type(5,  53) * ci_type(2, -23),  1219 * ce(7) );

    // division
    BOOST_CHECK_EQUAL( ce(0) / ce(0), ce(0) );
    BOOST_CHECK_EQUAL( ce(1) / ce(0), ce(1) );
    BOOST_CHECK_EQUAL( ce(0) / ce(1), -ce(1) );
    BOOST_CHECK_EQUAL( nce(2, true) / nce(3, false), nce(1, false) );
    BOOST_CHECK_EQUAL( cr_type(14, 7) / cr_type(9, -4), cr_type(7,
     ri_type( -7, 4 )) );

    BOOST_CHECK_THROW( cr_type(9, 7) / cr_type(1, 0), bad_rational );

    // mixed-type operations
    BOOST_CHECK_EQUAL( ce(0) * nce(0, true), nce(0, true) );
    BOOST_CHECK_EQUAL( nce(1, true) * ce(1), nce(0) );
    BOOST_CHECK_EQUAL( ce(2) * ci_type(2, 5), ci_type(0, -5) );
    BOOST_CHECK_EQUAL( ci_type(3, -7) * ce(4), ci_type(7, -7) );
    BOOST_CHECK_EQUAL( nce(5) * ci_type(7, 10), ci_type(2, 10) );
    BOOST_CHECK_EQUAL( ci_type(9, -11) * nce(14, true), ci_type(7, 11) );

    BOOST_CHECK_EQUAL( ce(0) / nce(0, true), nce(0, true) );
    BOOST_CHECK_EQUAL( nce(1, true) / ce(1), nce(0, true) );
    BOOST_CHECK_EQUAL( ce(2) / cr_type(2, 5), cr_type(0, ri_type(1, 5)) );
    BOOST_CHECK_EQUAL( cr_type(3, -7) / ce(4), cr_type(7, 7) );
    BOOST_CHECK_EQUAL( nce(5) / cr_type(7, 10), cr_type(2, ri_type(1, -10)) );
    BOOST_CHECK_EQUAL( cr_type(9, -11) / nce(14, true), cr_type(7, -11) );

    BOOST_CHECK_THROW( ce(7) / cr_type(3, 0), bad_rational );
    BOOST_CHECK_THROW( nce(11, true) / cr_type(0, 0), bad_rational );
}

// Scalar multiplication and division member functions
BOOST_AUTO_UNIT_TEST( cayley_scalar_muldiv_self_test )
{
    using boost::bad_rational;

    ci_type const  c1( 0, 0 ), c2( 0, 2 ), c3( 0, -3 );
    ci_type const  c4( 1, 0 ), c5( 2, 5 ), c6( 3, -7 );
    ci_type        t1 = c1, t2 = c2, t3 = c3, t4 = c4, t5 = c5, t6 = c6;

    t1 *= 11;  BOOST_CHECK_EQUAL( t1, c1 * 11 );
    t2 *= 11;  BOOST_CHECK_EQUAL( t2, c2 * 11 );
    t3 *= 11;  BOOST_CHECK_EQUAL( t3, c3 * 11 );
    t4 *= 11;  BOOST_CHECK_EQUAL( t4, c4 * 11 );
    t5 *= 11;  BOOST_CHECK_EQUAL( t5, c5 * 11 );
    t6 *= 11;  BOOST_CHECK_EQUAL( t6, c6 * 11 );

    t1 *= 0;  BOOST_CHECK_EQUAL( t1, ci_type(0) );
    t2 *= 0;  BOOST_CHECK_EQUAL( t2, ci_type(0) );
    t3 *= 0;  BOOST_CHECK_EQUAL( t3, ci_type(0) );
    t4 *= 0;  BOOST_CHECK_EQUAL( t4, ci_type(0) );
    t5 *= 0;  BOOST_CHECK_EQUAL( t5, ci_type(0) );
    t6 *= 0;  BOOST_CHECK_EQUAL( t6, ci_type(0) );

    cr_type const  c1a(c1), c2a(c2), c3a(c3), c4a(c4), c5a(c5), c6a(c6);
    cr_type        t1a = c1a, t2a = c2a, t3a = c3a;
    cr_type        t4a = c4a, t5a = c5a, t6a = c6a;

    t1a /= ri_type(-13);  BOOST_CHECK_EQUAL( t1a, c1a / ri_type(-13) );
    t2a /= ri_type(-13);  BOOST_CHECK_EQUAL( t2a, c2a / ri_type(-13) );
    t3a /= ri_type(-13);  BOOST_CHECK_EQUAL( t3a, c3a / ri_type(-13) );
    t4a /= ri_type(-13);  BOOST_CHECK_EQUAL( t4a, c4a / ri_type(-13) );
    t5a /= ri_type(-13);  BOOST_CHECK_EQUAL( t5a, c5a / ri_type(-13) );
    t6a /= ri_type(-13);  BOOST_CHECK_EQUAL( t6a, c6a / ri_type(-13) );

    BOOST_CHECK_THROW( t1a /= ri_type(0), bad_rational );
    BOOST_CHECK_THROW( t2a /= ri_type(0), bad_rational );
    BOOST_CHECK_THROW( t3a /= ri_type(0), bad_rational );
    BOOST_CHECK_THROW( t4a /= ri_type(0), bad_rational );
    BOOST_CHECK_THROW( t5a /= ri_type(0), bad_rational );
    BOOST_CHECK_THROW( t6a /= ri_type(0), bad_rational );
}

// Element multiplication and division member functions
BOOST_AUTO_UNIT_TEST( cayley_element_muldiv_self_test )
{
    // signed elements
    nce const  c1( 0 ), c2( 0, true ), c3( 1 ), c4( 2, true );
    nce        t1 = c1, t2 = c2, t3 = c3, t4 = c4;

    t1 *= c3;  BOOST_CHECK_EQUAL( t1, c1 * c3 );
    t2 *= c4;  BOOST_CHECK_EQUAL( t2, c2 * c4 );
    t3 *= c1;  BOOST_CHECK_EQUAL( t3, c3 * c1 );
    t4 *= c2;  BOOST_CHECK_EQUAL( t4, c4 * c2 );

    t1 = c1;  t1 /= c4;  BOOST_CHECK_EQUAL( t1, c1 / c4 );
    t2 = c2;  t2 /= c3;  BOOST_CHECK_EQUAL( t2, c2 / c3 );
    t3 = c3;  t3 /= c2;  BOOST_CHECK_EQUAL( t3, c3 / c2 );
    t4 = c4;  t4 /= c1;  BOOST_CHECK_EQUAL( t4, c4 / c1 );

    // scaled elements
    cr_type const  c5( 0, 15 ), c6( 1, -4 );
    cr_type const  c7( 2, ri_type(-11, 7) ), c8( 3, ri_type(2, 3) );
    cr_type        t5 = c5, t6 = c6, t7 = c7, t8 = c8;

    t5 *= c7;  BOOST_CHECK_EQUAL( t5, c5 * c7 );
    t6 *= c8;  BOOST_CHECK_EQUAL( t6, c6 * c8 );
    t7 *= c5;  BOOST_CHECK_EQUAL( t7, c7 * c5 );
    t8 *= c6;  BOOST_CHECK_EQUAL( t8, c8 * c6 );

    t5 = c5;  t5 /= c8;  BOOST_CHECK_EQUAL( t5, c5 / c8 );
    t6 = c6;  t6 /= c7;  BOOST_CHECK_EQUAL( t6, c6 / c7 );
    t7 = c7;  t7 /= c6;  BOOST_CHECK_EQUAL( t7, c7 / c6 );
    t8 = c8;  t8 /= c5;  BOOST_CHECK_EQUAL( t8, c8 / c5 );

    BOOST_CHECK_THROW( t7 /= cr_type(5, 0), boost::bad_rational );

    // mixed elements (via conversion)
    t2 = c2;  BOOST_CHECK_EQUAL( t2 /= ce(1), c2 / ce(1) );
    t3 = c3;  BOOST_CHECK_EQUAL( t3 *= ce(2), c3 * ce(2) );
    t6 = c6;  BOOST_CHECK_EQUAL( t6 /= ce(3), c6 / ce(3) );
    t7 = c7;  BOOST_CHECK_EQUAL( t7 *= ce(5), c7 * ce(5) );

    t5 = c5;  BOOST_CHECK_EQUAL( t5 /= nce(7), c5 / nce(7) );
    t5 = c5;  BOOST_CHECK_EQUAL( t5 *= nce(11, true), c5 * nce(11, true) );
    t8 = c8;  BOOST_CHECK_EQUAL( t8 /= ce(13), c8 / ce(13) );
    t8 = c8;  BOOST_CHECK_EQUAL( t8 *= ce(17), c8 * ce(17) );
}

// Shifting
BOOST_AUTO_UNIT_TEST( cayley_element_shift_op_test )
{
    // basic elements
    BOOST_CHECK_EQUAL( ce(0) << 2, ce(2) );
    BOOST_CHECK_EQUAL( ce(2) << 3, ce(5) );

    BOOST_CHECK_EQUAL( ce(7) >> 5, ce(2) );
    BOOST_CHECK_EQUAL( ce(0) >> 0, ce(0) );

    // signed elements
    BOOST_CHECK_EQUAL( nce(0)       <<  5, nce(5) );
    BOOST_CHECK_EQUAL( nce(3, true) << 11, nce(14, true) );

    BOOST_CHECK_EQUAL( nce(11, true) >> 7, nce(4, true) );
    BOOST_CHECK_EQUAL( nce( 5)       >> 2, nce(3) );

    // scaled elements
    BOOST_CHECK_EQUAL( ci_type(0,  2) << 5, ci_type( 5,  2) );
    BOOST_CHECK_EQUAL( ci_type(5, -1) << 7, ci_type(12, -1) );

    BOOST_CHECK_EQUAL( ci_type(17, -20) >> 11, ci_type(6, -20) );
    BOOST_CHECK_EQUAL( ci_type( 3, 101) >>  1, ci_type(2, 101) );
}

// Shifting member functions
BOOST_AUTO_UNIT_TEST( cayley_element_shift_self_test )
{
    // basic elements
    ce const  c1( 0 ), c2( 1 ), c3( 10 );
    ce        t1 = c1, t2 = c2, t3 = c3;

    t1 <<= 3;  BOOST_CHECK_EQUAL( t1, c1 << 3 );
    t2 <<= 5;  BOOST_CHECK_EQUAL( t2, c2 << 5 );
    t3 <<= 7;  BOOST_CHECK_EQUAL( t3, c3 << 7 );

    t1 = c1;  t1 >>= 0;  BOOST_CHECK_EQUAL( t1, c1 >> 0 );
    t2 = c2;  t2 >>= 1;  BOOST_CHECK_EQUAL( t2, c2 >> 1 );
    t3 = c3;  t3 >>= 7;  BOOST_CHECK_EQUAL( t3, c3 >> 7 );

    // signed elements
    nce const  c4( 0, true ), c5( 3 ), c6( 11, true );
    nce        t4 = c4, t5 = c5, t6 = c6;

    t4 <<= 3;  BOOST_CHECK_EQUAL( t4, c4 << 3 );
    t5 <<= 5;  BOOST_CHECK_EQUAL( t5, c5 << 5 );
    t6 <<= 7;  BOOST_CHECK_EQUAL( t6, c6 << 7 );

    t4 = c4;  t4 >>= 0;  BOOST_CHECK_EQUAL( t4, c4 >> 0 );
    t5 = c5;  t5 >>= 2;  BOOST_CHECK_EQUAL( t5, c5 >> 2 );
    t6 = c6;  t6 >>= 5;  BOOST_CHECK_EQUAL( t6, c6 >> 5 );

    // scaled elements
    cr_type const  c7( 0, ri_type(2, 5) ), c8( 5, -3 ), c9( 13, 20 );
    cr_type        t7 = c7, t8 = c8, t9 = c9;

    t7 <<= 3;  BOOST_CHECK_EQUAL( t7, c7 << 3 );
    t8 <<= 5;  BOOST_CHECK_EQUAL( t8, c8 << 5 );
    t9 <<= 7;  BOOST_CHECK_EQUAL( t9, c9 << 7 );

    t7 = c7;  t7 >>= 0;  BOOST_CHECK_EQUAL( t7, c7 >> 0 );
    t8 = c8;  t8 >>= 3;  BOOST_CHECK_EQUAL( t8, c8 >> 3 );
    t9 = c9;  t9 >>= 7;  BOOST_CHECK_EQUAL( t9, c9 >> 7 );
}

// Component functions
BOOST_AUTO_UNIT_TEST( cayley_element_component_test )
{
    ce const       c1( 0 ), c2( 1 ), c3( 2 );
    nce const      c4( 0 ), c5( 1, true ), c6( 3 );
    cr_type const  c7( 0, 4 ), c8( 1, -3 ), c9( 5, ri_type(2, 7) );

    // real component
    BOOST_CHECK_EQUAL( real(c1), 1 );
    BOOST_CHECK_EQUAL( real(c2), 0 );
    BOOST_CHECK_EQUAL( real(c3), 0 );
    BOOST_CHECK_EQUAL( real(c4), 1 );  BOOST_CHECK_EQUAL( real(-c4), -1 );
    BOOST_CHECK_EQUAL( real(c5), 0 );  BOOST_CHECK_EQUAL( real(-c5),  0 );
    BOOST_CHECK_EQUAL( real(c6), 0 );  BOOST_CHECK_EQUAL( real(-c6),  0 );
    BOOST_CHECK_EQUAL( real(c7), 4 );
    BOOST_CHECK_EQUAL( real(c8), 0 );
    BOOST_CHECK_EQUAL( real(c9), 0 );

    // (complex) imaginary component
    BOOST_CHECK_EQUAL( imag(c1),  0 );
    BOOST_CHECK_EQUAL( imag(c2),  1 );
    BOOST_CHECK_EQUAL( imag(c3),  0 );
    BOOST_CHECK_EQUAL( imag(c4),  0 );  BOOST_CHECK_EQUAL( imag(-c4), 0 );
    BOOST_CHECK_EQUAL( imag(c5), -1 );  BOOST_CHECK_EQUAL( imag(-c5), 1 );
    BOOST_CHECK_EQUAL( imag(c6),  0 );  BOOST_CHECK_EQUAL( imag(-c6), 0 );
    BOOST_CHECK_EQUAL( imag(c7),  0 );
    BOOST_CHECK_EQUAL( imag(c8), -3 );
    BOOST_CHECK_EQUAL( imag(c9),  0 );

    // unreal components
    BOOST_CHECK_EQUAL( unreal(c9), c9 );
    BOOST_CHECK_EQUAL( unreal(c8), c8 );
    BOOST_CHECK_EQUAL( unreal(c7), cr_type() );
    BOOST_CHECK_EQUAL( unreal(c6), ci_type(c6) );
    BOOST_CHECK_EQUAL( unreal(c5), ci_type(c5) );
    BOOST_CHECK_EQUAL( unreal(c4), ci_type() );
    BOOST_CHECK_EQUAL( unreal(c3), ci_type(c3) );
    BOOST_CHECK_EQUAL( unreal(c2), ci_type(c2) );
    BOOST_CHECK_EQUAL( unreal(c1), ci_type() );
}

// Element to integer scalar power
BOOST_AUTO_UNIT_TEST( cayley_element_integer_power_test )
{
    ce const   c1( 0 ), c2( 1 ), c3( 2 );
    nce const  c1a( c1 ), c2a( c2 ), c3a( c3 ), c4 = -c1, c5 = -c2, c6 = -c3;

    BOOST_CHECK_EQUAL( pow(c1,  0), ce(0) );
    BOOST_CHECK_EQUAL( pow(c1,  1), ce(0) );
    BOOST_CHECK_EQUAL( pow(c1,  2), ce(0) );
    BOOST_CHECK_EQUAL( pow(c1,  3), ce(0) );
    BOOST_CHECK_EQUAL( pow(c1,  4), ce(0) );
    BOOST_CHECK_EQUAL( pow(c1, -1), ce(0) );
    BOOST_CHECK_EQUAL( pow(c1, -2), ce(0) );
    BOOST_CHECK_EQUAL( pow(c1, -3), ce(0) );
    BOOST_CHECK_EQUAL( pow(c1, -4), ce(0) );

    BOOST_CHECK_EQUAL( pow(c2,  0), +ce(0) );
    BOOST_CHECK_EQUAL( pow(c2,  1), +ce(1) );
    BOOST_CHECK_EQUAL( pow(c2,  2), -ce(0) );
    BOOST_CHECK_EQUAL( pow(c2,  3), -ce(1) );
    BOOST_CHECK_EQUAL( pow(c2,  4), +ce(0) );
    BOOST_CHECK_EQUAL( pow(c2, -1), -ce(1) );
    BOOST_CHECK_EQUAL( pow(c2, -2), -ce(0) );
    BOOST_CHECK_EQUAL( pow(c2, -3), +ce(1) );
    BOOST_CHECK_EQUAL( pow(c2, -4), +ce(0) );

    BOOST_CHECK_EQUAL( pow(c3,  0), +ce(0) );
    BOOST_CHECK_EQUAL( pow(c3,  1), +ce(2) );
    BOOST_CHECK_EQUAL( pow(c3,  2), -ce(0) );
    BOOST_CHECK_EQUAL( pow(c3,  3), -ce(2) );
    BOOST_CHECK_EQUAL( pow(c3,  4), +ce(0) );
    BOOST_CHECK_EQUAL( pow(c3, -1), -ce(2) );
    BOOST_CHECK_EQUAL( pow(c3, -2), -ce(0) );
    BOOST_CHECK_EQUAL( pow(c3, -3), +ce(2) );
    BOOST_CHECK_EQUAL( pow(c3, -4), +ce(0) );

    BOOST_CHECK_EQUAL( pow(c1a,  0), ce(0) );
    BOOST_CHECK_EQUAL( pow(c1a,  1), ce(0) );
    BOOST_CHECK_EQUAL( pow(c1a,  2), ce(0) );
    BOOST_CHECK_EQUAL( pow(c1a,  3), ce(0) );
    BOOST_CHECK_EQUAL( pow(c1a,  4), ce(0) );
    BOOST_CHECK_EQUAL( pow(c1a, -1), ce(0) );
    BOOST_CHECK_EQUAL( pow(c1a, -2), ce(0) );
    BOOST_CHECK_EQUAL( pow(c1a, -3), ce(0) );
    BOOST_CHECK_EQUAL( pow(c1a, -4), ce(0) );

    BOOST_CHECK_EQUAL( pow(c2a,  0), +ce(0) );
    BOOST_CHECK_EQUAL( pow(c2a,  1), +ce(1) );
    BOOST_CHECK_EQUAL( pow(c2a,  2), -ce(0) );
    BOOST_CHECK_EQUAL( pow(c2a,  3), -ce(1) );
    BOOST_CHECK_EQUAL( pow(c2a,  4), +ce(0) );
    BOOST_CHECK_EQUAL( pow(c2a, -1), -ce(1) );
    BOOST_CHECK_EQUAL( pow(c2a, -2), -ce(0) );
    BOOST_CHECK_EQUAL( pow(c2a, -3), +ce(1) );
    BOOST_CHECK_EQUAL( pow(c2a, -4), +ce(0) );

    BOOST_CHECK_EQUAL( pow(c3a,  0), +ce(0) );
    BOOST_CHECK_EQUAL( pow(c3a,  1), +ce(2) );
    BOOST_CHECK_EQUAL( pow(c3a,  2), -ce(0) );
    BOOST_CHECK_EQUAL( pow(c3a,  3), -ce(2) );
    BOOST_CHECK_EQUAL( pow(c3a,  4), +ce(0) );
    BOOST_CHECK_EQUAL( pow(c3a, -1), -ce(2) );
    BOOST_CHECK_EQUAL( pow(c3a, -2), -ce(0) );
    BOOST_CHECK_EQUAL( pow(c3a, -3), +ce(2) );
    BOOST_CHECK_EQUAL( pow(c3a, -4), +ce(0) );

    BOOST_CHECK_EQUAL( pow(c4,  0), +nce(0) );
    BOOST_CHECK_EQUAL( pow(c4,  1), -nce(0) );
    BOOST_CHECK_EQUAL( pow(c4,  2), +nce(0) );
    BOOST_CHECK_EQUAL( pow(c4,  3), -nce(0) );
    BOOST_CHECK_EQUAL( pow(c4,  4), +nce(0) );
    BOOST_CHECK_EQUAL( pow(c4, -1), -nce(0) );
    BOOST_CHECK_EQUAL( pow(c4, -2), +nce(0) );
    BOOST_CHECK_EQUAL( pow(c4, -3), -nce(0) );
    BOOST_CHECK_EQUAL( pow(c4, -4), +nce(0) );

    BOOST_CHECK_EQUAL( pow(c5,  0), +nce(0) );
    BOOST_CHECK_EQUAL( pow(c5,  1), -nce(1) );
    BOOST_CHECK_EQUAL( pow(c5,  2), -nce(0) );
    BOOST_CHECK_EQUAL( pow(c5,  3), +nce(1) );
    BOOST_CHECK_EQUAL( pow(c5,  4), +nce(0) );
    BOOST_CHECK_EQUAL( pow(c5, -1), +nce(1) );
    BOOST_CHECK_EQUAL( pow(c5, -2), -nce(0) );
    BOOST_CHECK_EQUAL( pow(c5, -3), -nce(1) );
    BOOST_CHECK_EQUAL( pow(c5, -4), +nce(0) );

    BOOST_CHECK_EQUAL( pow(c6,  0), +nce(0) );
    BOOST_CHECK_EQUAL( pow(c6,  1), -nce(2) );
    BOOST_CHECK_EQUAL( pow(c6,  2), -nce(0) );
    BOOST_CHECK_EQUAL( pow(c6,  3), +nce(2) );
    BOOST_CHECK_EQUAL( pow(c6,  4), +nce(0) );
    BOOST_CHECK_EQUAL( pow(c6, -1), +nce(2) );
    BOOST_CHECK_EQUAL( pow(c6, -2), -nce(0) );
    BOOST_CHECK_EQUAL( pow(c6, -3), -nce(2) );
    BOOST_CHECK_EQUAL( pow(c6, -4), +nce(0) );

    cd_type const  c7( 0, 2.0 ), c8( 1, -3.0 ), c9( 2, 5.0 );
    cd_type        r;

    r = pow( c7, 0 );  BOOST_CHECK_EQUAL( r.basis(), 0u );
     BOOST_CHECK_CLOSE( r.scale(), 1.0, 0.1 );
    r = pow( c7, 1 );  BOOST_CHECK_EQUAL( r.basis(), 0u );
     BOOST_CHECK_CLOSE( r.scale(), 2.0, 0.1 );
    r = pow( c7, 2 );  BOOST_CHECK_EQUAL( r.basis(), 0u );
     BOOST_CHECK_CLOSE( r.scale(), 4.0, 0.1 );
    r = pow( c7, 3 );  BOOST_CHECK_EQUAL( r.basis(), 0u );
     BOOST_CHECK_CLOSE( r.scale(), 8.0, 0.1 );
    r = pow( c7, 4 );  BOOST_CHECK_EQUAL( r.basis(), 0u );
     BOOST_CHECK_CLOSE( r.scale(), 16.0, 0.1 );
    r = pow( c7, -1 );  BOOST_CHECK_EQUAL( r.basis(), 0u );
     BOOST_CHECK_CLOSE( r.scale(), 0.5, 0.1 );
    r = pow( c7, -2 );  BOOST_CHECK_EQUAL( r.basis(), 0u );
     BOOST_CHECK_CLOSE( r.scale(), 0.25, 0.1 );
    r = pow( c7, -3 );  BOOST_CHECK_EQUAL( r.basis(), 0u );
     BOOST_CHECK_CLOSE( r.scale(), 0.125, 0.1 );
    r = pow( c7, -4 );  BOOST_CHECK_EQUAL( r.basis(), 0u );
     BOOST_CHECK_CLOSE( r.scale(), 0.0625, 0.1 );

    r = pow( c8, 0 );  BOOST_CHECK_EQUAL( r.basis(), 0u );
     BOOST_CHECK_CLOSE( r.scale(), 1.0, 0.1 );
    r = pow( c8, 1 );  BOOST_CHECK_EQUAL( r.basis(), 1u );
     BOOST_CHECK_CLOSE( r.scale(), -3.0, 0.1 );
    r = pow( c8, 2 );  BOOST_CHECK_EQUAL( r.basis(), 0u );
     BOOST_CHECK_CLOSE( r.scale(), -9.0, 0.1 );
    r = pow( c8, 3 );  BOOST_CHECK_EQUAL( r.basis(), 1u );
     BOOST_CHECK_CLOSE( r.scale(), 27.0, 0.1 );
    r = pow( c8, 4 );  BOOST_CHECK_EQUAL( r.basis(), 0u );
     BOOST_CHECK_CLOSE( r.scale(), 81.0, 0.1 );
    r = pow( c8, -1 );  BOOST_CHECK_EQUAL( r.basis(), 1u );
     BOOST_CHECK_CLOSE( r.scale(), 1.0 / 3.0, 0.1 );
    r = pow( c8, -2 );  BOOST_CHECK_EQUAL( r.basis(), 0u );
     BOOST_CHECK_CLOSE( r.scale(), 1.0 / -9.0, 0.1 );
    r = pow( c8, -3 );  BOOST_CHECK_EQUAL( r.basis(), 1u );
     BOOST_CHECK_CLOSE( r.scale(), 1.0 / -27.0, 0.1 );
    r = pow( c8, -4 );  BOOST_CHECK_EQUAL( r.basis(), 0u );
     BOOST_CHECK_CLOSE( r.scale(), 1.0 / 81.0, 0.1 );

    r = pow( c9, 0 );  BOOST_CHECK_EQUAL( r.basis(), 0u );
     BOOST_CHECK_CLOSE( r.scale(), 1.0, 0.1 );
    r = pow( c9, 1 );  BOOST_CHECK_EQUAL( r.basis(), 2u );
     BOOST_CHECK_CLOSE( r.scale(), 5.0, 0.1 );
    r = pow( c9, 2 );  BOOST_CHECK_EQUAL( r.basis(), 0u );
     BOOST_CHECK_CLOSE( r.scale(), -25.0, 0.1 );
    r = pow( c9, 3 );  BOOST_CHECK_EQUAL( r.basis(), 2u );
     BOOST_CHECK_CLOSE( r.scale(), -125.0, 0.1 );
    r = pow( c9, 4 );  BOOST_CHECK_EQUAL( r.basis(), 0u );
     BOOST_CHECK_CLOSE( r.scale(), 625.0, 0.1 );
    r = pow( c9, -1 );  BOOST_CHECK_EQUAL( r.basis(), 2u );
     BOOST_CHECK_CLOSE( r.scale(), -0.2, 0.1 );
    r = pow( c9, -2 );  BOOST_CHECK_EQUAL( r.basis(), 0u );
     BOOST_CHECK_CLOSE( r.scale(), -0.04, 0.1 );
    r = pow( c9, -3 );  BOOST_CHECK_EQUAL( r.basis(), 2u );
     BOOST_CHECK_CLOSE( r.scale(), 0.008, 0.1 );
    r = pow( c9, -4 );  BOOST_CHECK_EQUAL( r.basis(), 0u );
     BOOST_CHECK_CLOSE( r.scale(), 0.0016, 0.1 );
}
