//  Boost big_whole_test.cpp test file  --------------------------------------//

//  Copyright 2004 Daryle Walker.  Use, modification, and distribution are
//  subject to the Boost Software License, Version 1.0.  (See accompanying file
//  LICENSE_1_0.txt or a copy at <http://www.boost.org/LICENSE_1_0.txt>.)

//  See <http://www.boost.org/libs/math/> for the library's home page.

//  Revision History
//   16 Feb 2004  Initial version (Daryle Walker)

#include <boost/math/big_whole.hpp>  // for boost::math::big_whole, etc.
#include <boost/test/unit_test.hpp>  // for main, BOOST_CHECK_EQUAL, etc.

#include <algorithm>  // for std::copy
#include <cstddef>    // for std::size_t
#include <limits>     // for std::numeric_limits
#include <set>        // for std::set
#include <stdexcept>  // for std::range_error, std::length_error, etc.
#include <valarray>   // for std::valarray


// Control to allow extremely large memory allocations
#ifdef CONTROL_USE_BIG_MEMORY
#define PRIVATE_USE_BIG_MEMORY  1
#else
#define PRIVATE_USE_BIG_MEMORY  0
#endif


// Use internal knowledge of big_whole (i.e. cheat) to force situations
// where multiple-word representations have to be used
typedef unsigned int                       word_type;
typedef std::numeric_limits<word_type>  wlimits_type;


// Helper function to compare valarrays
template < typename T >
bool
equal_valarrays
(
    std::valarray<T> const &  lhs,
    std::valarray<T> const &  rhs
)
{
    if ( lhs.size() == rhs.size() )
    {
        std::valarray<std::size_t>  s( lhs.size() );

        s[ lhs == rhs ] = 1;
        return s.sum() == s.size();
    }
    else
    {
        return false;
    }
}

// Helper function to insert values into sets
template < class SetType, typename ValueType >
void
insert_value_range
(
    SetType &  s,
    ValueType  start,
    ValueType  finish
)
{
    for ( ValueType i = start ; i <= finish ; ++i )
    {
        s.insert( i );
    }
}

// Helper function to remove values from sets
template < class SetType, typename ValueType >
void
erase_value_range
(
    SetType &  s,
    ValueType  start,
    ValueType  finish
)
{
    for ( ValueType i = start ; i <= finish ; ++i )
    {
        s.erase( i );
    }
}

// Helper function to convert sets to valarrays
template < typename T >
std::valarray<T>
set_to_valarray
(
    std::set<T> const &  s
)
{
    std::valarray<T>  temp( s.size() );

    std::copy( s.begin(), s.end(), &temp[0] );
    return temp;
}


// Unit test for the basics
void
basic_bigwhole_unit_test
(
)
{
    using boost::math::big_whole;
    using std::valarray;
    using std::size_t;

    typedef valarray<bool>    va_bool_t;
    typedef valarray<size_t>  va_size_t;

    // Default construction
    big_whole  x1;
    BOOST_CHECK_EQUAL( 0u, x1.to_uintmax() );

    // Converting assignment
    x1.assign( 5u );
    BOOST_CHECK_EQUAL( 5u, x1.to_uintmax() );

    // Converting construction
    big_whole  x2 = 17;
    BOOST_CHECK_EQUAL( 17u, x2.to_uintmax() );

    // Copy construction
    big_whole  x3( x1 );
    BOOST_CHECK_EQUAL( 5u, x3.to_uintmax() );

    // Assignment operator
    x1 = x2;
    BOOST_CHECK_EQUAL( 17u, x1.to_uintmax() );

    // Swapping
    swap( x1, x3 );
    BOOST_CHECK_EQUAL( 5u, x1.to_uintmax() );
    BOOST_CHECK_EQUAL( 17u, x3.to_uintmax() );

    // Copying assignment
    x2.assign( big_whole() );
    BOOST_CHECK_EQUAL( 0u, x2.to_uintmax() );

    // Bit-vector conversion
    va_bool_t const  x1_b = x1.to_bit_vector();
    bool const       x1_b_check[] = { true, false, true };
    size_t const     x1_b_size = sizeof( x1_b_check ) / sizeof( x1_b_check[0] );
    BOOST_CHECK( equal_valarrays(va_bool_t( x1_b_check, x1_b_size ), x1_b) );

    BOOST_CHECK_EQUAL( 0u, x2.to_bit_vector().size() );

    va_bool_t const  x3_b = x3.to_bit_vector();
    bool const       x3_b_check[] = { true, false, false, false, true };
    size_t const     x3_b_size = sizeof( x3_b_check ) / sizeof( x3_b_check[0] );
    BOOST_CHECK( equal_valarrays(va_bool_t( x3_b_check, x3_b_size ), x3_b) );

    // Bit-index conversion
    va_size_t const  x1_i = x1.to_bit_indices();
    size_t const     x1_i_check[] = { 0, 2 };
    size_t const     x1_i_size = sizeof( x1_i_check ) / sizeof( x1_i_check[0] );
    BOOST_CHECK( equal_valarrays(va_size_t( x1_i_check, x1_i_size ), x1_i) );
    BOOST_CHECK_EQUAL( x1_b_size, 1u + x1_i.max() );

    BOOST_CHECK_EQUAL( 0u, x2.to_bit_indices().size() );

    va_size_t const  x3_i = x3.to_bit_indices();
    size_t const     x3_i_check[] = { 0, 4 };
    size_t const     x3_i_size = sizeof( x3_i_check ) / sizeof( x3_i_check[0] );
    BOOST_CHECK( equal_valarrays(va_size_t( x3_i_check, x3_i_size ), x3_i) );
    BOOST_CHECK_EQUAL( x3_b_size, 1u + x3_i.max() );

    // Bit-vector construction and assignment
    big_whole  x4( x1_b );
    BOOST_CHECK_EQUAL( 5u, x4.to_uintmax() );

    x4.reconfigure( x3_b );
    BOOST_CHECK_EQUAL( 17u, x4.to_uintmax() );

    x4.reconfigure( va_bool_t() );
    BOOST_CHECK_EQUAL( 0u, x4.to_uintmax() );

    // Bit-index construction and assignment
    big_whole  x5( x3_i );
    BOOST_CHECK_EQUAL( 17u, x5.to_uintmax() );

    x5.reconfigure( x1_i );
    BOOST_CHECK_EQUAL( 5u, x5.to_uintmax() );

    x5.reconfigure( va_size_t() );
    BOOST_CHECK_EQUAL( 0u, x5.to_uintmax() );

    // Minimum-required bit length
    BOOST_CHECK_EQUAL( x1_b_size, x1.length() );
    BOOST_CHECK_EQUAL( 0u, x2.length() );
    BOOST_CHECK_EQUAL( x3_b_size, x3.length() );

    // Bit count
    BOOST_CHECK_EQUAL( x1_i_size, x1.count() );
    BOOST_CHECK_EQUAL( 0u, x2.count() );
    BOOST_CHECK_EQUAL( x3_i_size, x3.count() );

    BOOST_CHECK( x1.any() );
    BOOST_CHECK( !x2.any() );
    BOOST_CHECK( x3.any() );

    BOOST_CHECK( !x1.none() );
    BOOST_CHECK( x2.none() );
    BOOST_CHECK( !x3.none() );

    // Bit testing
    BOOST_CHECK( x1.test(0) && !x1.test(1) && x1.test(2) && !x1.test(3)
     && !x1.test(4) && !x1.test(5) && !x1.test(wlimits_type::digits) );
    BOOST_CHECK( !x2.test(0) && !x2.test(1) && !x2.test(2) && !x2.test(3)
     && !x2.test(4) && !x2.test(5) && !x2.test(wlimits_type::digits) );
    BOOST_CHECK( x3.test(0) && !x3.test(1) && !x3.test(2) && !x3.test(3)
     && x3.test(4) && !x3.test(5) && !x3.test(wlimits_type::digits) );

    // Boolean test
    BOOST_CHECK( x1 );
    BOOST_CHECK( !x2 );
    BOOST_CHECK( x3 );
}

// Unit test for "tests"
void
bigwhole_multi_bit_check_unit_test
(
)
{
    using boost::math::big_whole;
    using std::size_t;

    size_t const  size_max = std::numeric_limits<size_t>::max();

    // Non-zero tests
    big_whole const  x1( 74u );
    size_t const     l1 = x1.length();

    BOOST_CHECK_EQUAL( 7u, l1 );

    BOOST_CHECK_EQUAL( 74u, x1.tests(0, size_max).to_uintmax() );
    BOOST_CHECK_EQUAL( 74u, x1.tests(0, l1).to_uintmax() );
    BOOST_CHECK_EQUAL( 74u, x1.tests(0, l1 - 1).to_uintmax() );
    BOOST_CHECK_EQUAL( 10u, x1.tests(0, 5).to_uintmax() );
    BOOST_CHECK_EQUAL( 10u, x1.tests(0, 4).to_uintmax() );
    BOOST_CHECK_EQUAL( 10u, x1.tests(0, 3).to_uintmax() );
    BOOST_CHECK_EQUAL( 2u, x1.tests(0, 2).to_uintmax() );
    BOOST_CHECK_EQUAL( 2u, x1.tests(0, 1).to_uintmax() );
    BOOST_CHECK_EQUAL( 0u, x1.tests(0, 0).to_uintmax() );

    BOOST_CHECK_EQUAL( 37u, x1.tests(1, l1).to_uintmax() );
    BOOST_CHECK_EQUAL( 18u, x1.tests(2, l1).to_uintmax() );
    BOOST_CHECK_EQUAL( 9u, x1.tests(3, l1).to_uintmax() );
    BOOST_CHECK_EQUAL( 4u, x1.tests(4, l1).to_uintmax() );
    BOOST_CHECK_EQUAL( 2u, x1.tests(5, l1).to_uintmax() );
    BOOST_CHECK_EQUAL( 1u, x1.tests(l1 - 1, l1).to_uintmax() );
    BOOST_CHECK_EQUAL( 0u, x1.tests(l1, l1).to_uintmax() );

    BOOST_CHECK( x1.tests(1, 1) );
    BOOST_CHECK( !x1.tests(2, 2) );
    BOOST_CHECK( x1.tests(3, 3) );
    BOOST_CHECK( !x1.tests(4, 4) );
    BOOST_CHECK( !x1.tests(5, 5) );
    BOOST_CHECK( x1.tests(l1 - 1, l1 - 1) );

    // Zero tests
    big_whole const  x2;

    BOOST_CHECK( !x2.tests(0, 4) );
    BOOST_CHECK( !x2.tests(2, size_max) );
    BOOST_CHECK( !x2.tests(0, size_max) );
    BOOST_CHECK( !x2.tests(3, 3) );
}

// Unit test for reversing
void
bigwhole_reverse_unit_test
(
)
{
    using boost::math::big_whole;
    using std::size_t;

    typedef std::valarray<size_t>  va_size_t;

    // Non-zero tests
    big_whole const  x1( 1 );

    BOOST_CHECK_EQUAL( 1u, x1.reverse().to_uintmax() );
    BOOST_CHECK_EQUAL( 1u, x1.length() );

    BOOST_CHECK_EQUAL( 1u, x1.reverse(0).to_uintmax() );
    BOOST_CHECK_EQUAL( 2u, x1.reverse(1).to_uintmax() );
    BOOST_CHECK_EQUAL( 4u, x1.reverse(2).to_uintmax() );
    BOOST_CHECK_EQUAL( 128u, x1.reverse(7).to_uintmax() );

    big_whole const  x2( 5 );

    BOOST_CHECK_EQUAL( 5u, x2.reverse().to_uintmax() );
    BOOST_CHECK_EQUAL( 3u, x2.length() );

    BOOST_CHECK_EQUAL( 1u, x2.reverse(0).to_uintmax() );
    BOOST_CHECK_EQUAL( 2u, x2.reverse(1).to_uintmax() );
    BOOST_CHECK_EQUAL( 5u, x2.reverse(2).to_uintmax() );
    BOOST_CHECK_EQUAL( 10u, x2.reverse(3).to_uintmax() );
    BOOST_CHECK_EQUAL( 20u, x2.reverse(4).to_uintmax() );
    BOOST_CHECK_EQUAL( 160u, x2.reverse(7).to_uintmax() );

    big_whole const  x3( 74 );

    BOOST_CHECK_EQUAL( 41u, x3.reverse().to_uintmax() );
    BOOST_CHECK_EQUAL( 7u, x3.length() );

    BOOST_CHECK_EQUAL( 0u, x3.reverse(0).to_uintmax() );
    BOOST_CHECK_EQUAL( 1u, x3.reverse(1).to_uintmax() );
    BOOST_CHECK_EQUAL( 2u, x3.reverse(2).to_uintmax() );
    BOOST_CHECK_EQUAL( 5u, x3.reverse(3).to_uintmax() );
    BOOST_CHECK_EQUAL( 10u, x3.reverse(4).to_uintmax() );
    BOOST_CHECK_EQUAL( 20u, x3.reverse(5).to_uintmax() );
    BOOST_CHECK_EQUAL( 41u, x3.reverse(6).to_uintmax() );
    BOOST_CHECK_EQUAL( 82u, x3.reverse(7).to_uintmax() );
    BOOST_CHECK_EQUAL( 164u, x3.reverse(8).to_uintmax() );
    BOOST_CHECK_EQUAL( 656u, x3.reverse(10).to_uintmax() );

    // Zero tests
    big_whole const  x4;

    BOOST_CHECK( !x4.length() );
    BOOST_CHECK( !x4.reverse() );
    BOOST_CHECK( !x4.reverse(0) );
    BOOST_CHECK( !x4.reverse(2 * wlimits_type::digits) );

    // Multi-word tests
    size_t const     x5_i[] = { wlimits_type::digits - 1, wlimits_type::digits + 1 };
    size_t const     x5_s = sizeof( x5_i ) / sizeof( x5_i[0] );
    va_size_t const  x5_v( x5_i, x5_s );
    big_whole        x5( x5_v );

    BOOST_CHECK_EQUAL( 5u, x5.reverse().to_uintmax() );

    BOOST_CHECK_EQUAL( 5u, x5.reverse(wlimits_type::digits + 1).to_uintmax() );
    BOOST_CHECK_EQUAL( 2u, x5.reverse(wlimits_type::digits).to_uintmax() );
    BOOST_CHECK_EQUAL( 1u, x5.reverse(wlimits_type::digits - 1).to_uintmax() );
    BOOST_CHECK_EQUAL( 0u, x5.reverse(wlimits_type::digits - 2).to_uintmax() );
    BOOST_CHECK_EQUAL( 0u, x5.reverse(wlimits_type::digits - 3).to_uintmax() );
    BOOST_CHECK_EQUAL( 10u, x5.reverse(wlimits_type::digits + 2).to_uintmax() );
    BOOST_CHECK_EQUAL( 20u, x5.reverse(wlimits_type::digits + 3).to_uintmax() );
}

// Unit test for resetting every bit
void
bigwhole_all_bit_reset_unit_test
(
)
{
    using boost::math::big_whole;
    using std::size_t;

    typedef std::valarray<size_t>  va_size_t;

    // zero
    big_whole  a1;

    BOOST_CHECK( !a1 );

    a1.reset();
    BOOST_CHECK( !a1 );

    // one bit set
    big_whole  a2( 8 );

    BOOST_CHECK( a2 );

    a2.reset();
    BOOST_CHECK( !a2 );

    // multiple bits set
    big_whole  a3( 25 );

    BOOST_CHECK( a3 );

    a3.reset();
    BOOST_CHECK( !a3 );

    // two words
    size_t const  a4_i[] = { 0, wlimits_type::digits + 1 };
    size_t const  a4_s = sizeof( a4_i ) / sizeof( a4_i[0] );
    big_whole     a4( va_size_t(a4_i, a4_s) );

    BOOST_CHECK( a4 );

    a4.reset();
    BOOST_CHECK( !a4 );

    // more-than-two words
    size_t const  a5_i[] = { 1, wlimits_type::digits + 3, wlimits_type::digits + 4, 2 * wlimits_type::digits + 5 };
    size_t const  a5_s = sizeof( a5_i ) / sizeof( a5_i[0] );
    big_whole     a5( va_size_t(a5_i, a5_s) );

    BOOST_CHECK( a5 );

    a5.reset();
    BOOST_CHECK( !a5 );
}

// Unit test for resetting single bits
void
bigwhole_single_bit_reset_unit_test
(
)
{
    using boost::math::big_whole;
    using std::size_t;

    typedef std::valarray<size_t>  va_size_t;

    // zero
    big_whole  a1;

    BOOST_CHECK_EQUAL( 0u, a1.to_uintmax() );

    a1.reset( 3 );
    BOOST_CHECK_EQUAL( 0u, a1.to_uintmax() );

    // one bit set
    big_whole  a2( 8 );

    BOOST_CHECK_EQUAL( 8u, a2.to_uintmax() );

    a2.reset( 2 );
    BOOST_CHECK_EQUAL( 8u, a2.to_uintmax() );

    a2.reset( 3 );
    BOOST_CHECK_EQUAL( 0u, a2.to_uintmax() );

    // multiple bits set
    big_whole  a3( 25 );

    BOOST_CHECK_EQUAL( 25u, a3.to_uintmax() );

    a3.reset( 0 );
    BOOST_CHECK_EQUAL( 24u, a3.to_uintmax() );

    a3.reset( 1 );
    BOOST_CHECK_EQUAL( 24u, a3.to_uintmax() );

    a3.reset( 4 );
    BOOST_CHECK_EQUAL( 8u, a3.to_uintmax() );

    // two words
    size_t const     a4_old_i[] = { 0, wlimits_type::digits + 1 };
    size_t const     a4_old_s = sizeof( a4_old_i ) / sizeof( a4_old_i[0] );
    va_size_t const  a4_old( a4_old_i, a4_old_s );
    big_whole        a4( a4_old );

    BOOST_CHECK( equal_valarrays(a4_old, a4.to_bit_indices()) );

    a4.reset( 5 );
    BOOST_CHECK( equal_valarrays(a4_old, a4.to_bit_indices()) );

    a4.reset( 0 );
    BOOST_CHECK( equal_valarrays(va_size_t( a4_old_i + 1, a4_old_s - 1 ), a4.to_bit_indices()) );

    // more-than-two words
    size_t const     a5_old_i[] = { 1, wlimits_type::digits + 3, wlimits_type::digits + 4, 2 * wlimits_type::digits + 5 };
    size_t const     a5_old_s = sizeof( a5_old_i ) / sizeof( a5_old_i[0] );
    va_size_t const  a5_old( a5_old_i, a5_old_s );
    big_whole        a5( a5_old );

    BOOST_CHECK( equal_valarrays(a5_old, a5.to_bit_indices()) );

    a5.reset( 4 );
    BOOST_CHECK( equal_valarrays(a5_old, a5.to_bit_indices()) );

    a5.reset( 2 * wlimits_type::digits + 5 );
    BOOST_CHECK( equal_valarrays(va_size_t( a5_old_i, a5_old_s - 1 ), a5.to_bit_indices()) );
}

// Unit test for resetting a group of bits
void
bigwhole_group_bit_reset_unit_test
(
)
{
    using boost::math::big_whole;
    using std::size_t;

    typedef std::valarray<size_t>  va_size_t;

    // zero
    big_whole  a1;

    BOOST_CHECK_EQUAL( 0u, a1.to_uintmax() );

    a1.reset( 3, 7 );
    BOOST_CHECK_EQUAL( 0u, a1.to_uintmax() );

    // one bit set
    big_whole  a2( 8 );

    BOOST_CHECK_EQUAL( 8u, a2.to_uintmax() );

    a2.reset( 6, 9 );
    BOOST_CHECK_EQUAL( 8u, a2.to_uintmax() );

    a2.reset( 2, 3 );
    BOOST_CHECK_EQUAL( 0u, a2.to_uintmax() );

    // multiple bits set
    big_whole  a3( 25 );

    BOOST_CHECK_EQUAL( 25u, a3.to_uintmax() );

    a3.reset( 3, 3 );
    BOOST_CHECK_EQUAL( 17u, a3.to_uintmax() );

    a3.reset( 1, 2 );
    BOOST_CHECK_EQUAL( 17u, a3.to_uintmax() );

    a3.reset( 2, 6 );
    BOOST_CHECK_EQUAL( 1u, a3.to_uintmax() );

    // two words
    size_t const     a4_old_i[] = { 0, wlimits_type::digits + 1 };
    size_t const     a4_old_s = sizeof( a4_old_i ) / sizeof( a4_old_i[0] );
    va_size_t const  a4_old( a4_old_i, a4_old_s );
    big_whole        a4( a4_old );

    BOOST_CHECK( equal_valarrays(a4_old, a4.to_bit_indices()) );

    a4.reset( 5, 12 );
    BOOST_CHECK( equal_valarrays(a4_old, a4.to_bit_indices()) );

    a4.reset( 9, 2 * wlimits_type::digits );
    BOOST_CHECK( equal_valarrays(va_size_t( a4_old_i, a4_old_s - 1 ), a4.to_bit_indices()) );

    // more-than-two words
    size_t const     a5_old_i[] = { 1, wlimits_type::digits + 3, wlimits_type::digits + 4, 2 * wlimits_type::digits + 5 };
    size_t const     a5_old_s = sizeof( a5_old_i ) / sizeof( a5_old_i[0] );
    va_size_t const  a5_old( a5_old_i, a5_old_s );
    big_whole        a5( a5_old );

    BOOST_CHECK( equal_valarrays(a5_old, a5.to_bit_indices()) );

    a5.reset( 3, 12 );
    BOOST_CHECK( equal_valarrays(a5_old, a5.to_bit_indices()) );

    size_t const  a5_new_i[] = { 1, 2 * wlimits_type::digits + 5 };
    size_t const  a5_new_s = sizeof( a5_new_i ) / sizeof( a5_new_i[0] );

    a5.reset( wlimits_type::digits - 1, 2 * wlimits_type::digits + 1 );
    BOOST_CHECK( equal_valarrays(va_size_t( a5_new_i, a5_new_s ), a5.to_bit_indices()) );

    a5.reset( 2 * wlimits_type::digits + 1, 3 * wlimits_type::digits );
    BOOST_CHECK( equal_valarrays(va_size_t( 1u, 1 ), a5.to_bit_indices()) );

    a5.reset( wlimits_type::digits + 1, 3 * wlimits_type::digits );
    BOOST_CHECK( equal_valarrays(va_size_t( 1u, 1 ), a5.to_bit_indices()) );
}

// Unit test for setting single bits
void
bigwhole_single_bit_set_unit_test
(
)
{
    using boost::math::big_whole;
    using std::size_t;

    typedef std::valarray<size_t>  va_size_t;

    // zero
    big_whole  a1;

    BOOST_CHECK_EQUAL( 0u, a1.to_uintmax() );

    a1.set( 3 );
    BOOST_CHECK_EQUAL( 8u, a1.to_uintmax() );

    // one bit set
    big_whole  a2( 8 );

    BOOST_CHECK_EQUAL( 8u, a2.to_uintmax() );

    a2.set( 2 );
    BOOST_CHECK_EQUAL( 12u, a2.to_uintmax() );

    a2.set( 3 );
    BOOST_CHECK_EQUAL( 12u, a2.to_uintmax() );

    // multiple bits set
    big_whole  a3( 25 );

    BOOST_CHECK_EQUAL( 25u, a3.to_uintmax() );

    a3.set( 0 );
    BOOST_CHECK_EQUAL( 25u, a3.to_uintmax() );

    a3.set( 1 );
    BOOST_CHECK_EQUAL( 27u, a3.to_uintmax() );

    a3.set( 4 );
    BOOST_CHECK_EQUAL( 27u, a3.to_uintmax() );

    // two words
    size_t const     a4_old_i[] = { 0, wlimits_type::digits + 1 };
    size_t const     a4_old_s = sizeof( a4_old_i ) / sizeof( a4_old_i[0] );
    va_size_t const  a4_old( a4_old_i, a4_old_s );
    big_whole        a4( a4_old );

    BOOST_CHECK( equal_valarrays(a4_old, a4.to_bit_indices()) );

    size_t const     a4_new_i[] = { 0, 5, wlimits_type::digits + 1 };
    size_t const     a4_new_s = sizeof( a4_new_i ) / sizeof( a4_new_i[0] );
    va_size_t const  a4_new( a4_new_i, a4_new_s );

    a4.set( 5 );
    BOOST_CHECK( equal_valarrays(a4_new, a4.to_bit_indices()) );

    a4.set( 0 );
    BOOST_CHECK( equal_valarrays(a4_new, a4.to_bit_indices()) );

    // more-than-two words
    size_t const     a5_old_i[] = { 1, wlimits_type::digits + 3, wlimits_type::digits + 4, 2 * wlimits_type::digits + 5 };
    size_t const     a5_old_s = sizeof( a5_old_i ) / sizeof( a5_old_i[0] );
    va_size_t const  a5_old( a5_old_i, a5_old_s );
    big_whole        a5( a5_old );

    BOOST_CHECK( equal_valarrays(a5_old, a5.to_bit_indices()) );

    size_t const     a5_new_i[] = { 1, 4, wlimits_type::digits + 3, wlimits_type::digits + 4, 2 * wlimits_type::digits + 5 };
    size_t const     a5_new_s = sizeof( a5_new_i ) / sizeof( a5_new_i[0] );
    va_size_t const  a5_new( a5_new_i, a5_new_s );

    a5.set( 4 );
    BOOST_CHECK( equal_valarrays(a5_new, a5.to_bit_indices()) );

    a5.set( 2 * wlimits_type::digits + 5 );
    BOOST_CHECK( equal_valarrays(a5_new, a5.to_bit_indices()) );
}

// Unit test for setting a group of bits
void
bigwhole_group_bit_set_unit_test
(
)
{
    using boost::math::big_whole;
    using std::size_t;

    typedef std::valarray<size_t>  va_size_t;
    typedef std::set<size_t>       st_size_t;

    // zero
    big_whole  a1;

    BOOST_CHECK_EQUAL( 0u, a1.to_uintmax() );

    a1.set( 3, 7 );
    BOOST_CHECK_EQUAL( 248u, a1.to_uintmax() );

    // one bit set
    big_whole  a2( 8 );

    BOOST_CHECK_EQUAL( 8u, a2.to_uintmax() );

    a2.set( 6, 9 );
    BOOST_CHECK_EQUAL( 968u, a2.to_uintmax() );

    a2.set( 2, 3 );
    BOOST_CHECK_EQUAL( 972u, a2.to_uintmax() );

    // multiple bits set
    big_whole  a3( 25 );

    BOOST_CHECK_EQUAL( 25u, a3.to_uintmax() );

    a3.set( 3, 3 );
    BOOST_CHECK_EQUAL( 25u, a3.to_uintmax() );

    a3.set( 1, 2 );
    BOOST_CHECK_EQUAL( 31u, a3.to_uintmax() );

    a3.set( 2, 6 );
    BOOST_CHECK_EQUAL( 127u, a3.to_uintmax() );

    // two words
    size_t const     a4_old_i[] = { 0, wlimits_type::digits + 1 };
    size_t const     a4_old_s = sizeof( a4_old_i ) / sizeof( a4_old_i[0] );
    va_size_t const  a4_old( a4_old_i, a4_old_s );
    big_whole        a4( a4_old );
    st_size_t        a4_new( a4_old_i, a4_old_i + a4_old_s );

    BOOST_CHECK( equal_valarrays(a4_old, a4.to_bit_indices()) );

    a4.set( 5, 12 );
    insert_value_range( a4_new, 5, 12 );
    BOOST_CHECK( equal_valarrays(set_to_valarray( a4_new ), a4.to_bit_indices()) );

    a4.set( 9, 2 * wlimits_type::digits );
    insert_value_range( a4_new, 9, 2 * wlimits_type::digits );
    BOOST_CHECK( equal_valarrays(set_to_valarray( a4_new ), a4.to_bit_indices()) );

    // more-than-two words
    size_t const     a5_old_i[] = { 1, wlimits_type::digits + 3, wlimits_type::digits + 4, 2 * wlimits_type::digits + 5 };
    size_t const     a5_old_s = sizeof( a5_old_i ) / sizeof( a5_old_i[0] );
    va_size_t const  a5_old( a5_old_i, a5_old_s );
    big_whole        a5( a5_old );
    st_size_t        a5_new( a5_old_i, a5_old_i + a5_old_s );

    BOOST_CHECK( equal_valarrays(a5_old, a5.to_bit_indices()) );

    a5.set( 3, 12 );
    insert_value_range( a5_new, 3, 12 );
    BOOST_CHECK( equal_valarrays(set_to_valarray( a5_new ), a5.to_bit_indices()) );

    a5.set( wlimits_type::digits - 1, 2 * wlimits_type::digits + 1 );
    insert_value_range( a5_new, wlimits_type::digits - 1, 2 * wlimits_type::digits + 1 );
    BOOST_CHECK( equal_valarrays(set_to_valarray( a5_new ), a5.to_bit_indices()) );

    a5.set( 2 * wlimits_type::digits + 1, 3 * wlimits_type::digits );
    insert_value_range( a5_new, 2 * wlimits_type::digits + 1, 3 * wlimits_type::digits );
    BOOST_CHECK( equal_valarrays(set_to_valarray( a5_new ), a5.to_bit_indices()) );

    a5.set( wlimits_type::digits + 1, 3 * wlimits_type::digits );
    insert_value_range( a5_new, wlimits_type::digits + 1, 3 * wlimits_type::digits );
    BOOST_CHECK( equal_valarrays(set_to_valarray( a5_new ), a5.to_bit_indices()) );
}

// Unit test for flipping single bits
void
bigwhole_single_bit_flip_unit_test
(
)
{
    using boost::math::big_whole;
    using std::size_t;

    typedef std::valarray<size_t>  va_size_t;

    // zero
    big_whole  a1;

    BOOST_CHECK_EQUAL( 0u, a1.to_uintmax() );

    a1.flip( 3 );
    BOOST_CHECK_EQUAL( 8u, a1.to_uintmax() );

    // one bit set
    big_whole  a2( 8 );

    BOOST_CHECK_EQUAL( 8u, a2.to_uintmax() );

    a2.flip( 2 );
    BOOST_CHECK_EQUAL( 12u, a2.to_uintmax() );

    a2.flip( 3 );
    BOOST_CHECK_EQUAL( 4u, a2.to_uintmax() );

    // multiple bits set
    big_whole  a3( 25 );

    BOOST_CHECK_EQUAL( 25u, a3.to_uintmax() );

    a3.flip( 0 );
    BOOST_CHECK_EQUAL( 24u, a3.to_uintmax() );

    a3.flip( 1 );
    BOOST_CHECK_EQUAL( 26u, a3.to_uintmax() );

    a3.flip( 4 );
    BOOST_CHECK_EQUAL( 10u, a3.to_uintmax() );

    // two words
    size_t const     a4_old_i[] = { 0, wlimits_type::digits + 1 };
    size_t const     a4_old_s = sizeof( a4_old_i ) / sizeof( a4_old_i[0] );
    va_size_t const  a4_old( a4_old_i, a4_old_s );
    big_whole        a4( a4_old );

    BOOST_CHECK( equal_valarrays(a4_old, a4.to_bit_indices()) );

    size_t const     a4_new_i[] = { 0, 5, wlimits_type::digits + 1 };
    size_t const     a4_new_s = sizeof( a4_new_i ) / sizeof( a4_new_i[0] );
    va_size_t const  a4_new1( a4_new_i, a4_new_s );
    va_size_t const  a4_new2( a4_new_i + 1, a4_new_s - 1 );

    a4.flip( 5 );
    BOOST_CHECK( equal_valarrays(a4_new1, a4.to_bit_indices()) );

    a4.flip( 0 );
    BOOST_CHECK( equal_valarrays(a4_new2, a4.to_bit_indices()) );

    // more-than-two words
    size_t const     a5_old_i[] = { 1, wlimits_type::digits + 3, wlimits_type::digits + 4, 2 * wlimits_type::digits + 5 };
    size_t const     a5_old_s = sizeof( a5_old_i ) / sizeof( a5_old_i[0] );
    va_size_t const  a5_old( a5_old_i, a5_old_s );
    big_whole        a5( a5_old );

    BOOST_CHECK( equal_valarrays(a5_old, a5.to_bit_indices()) );

    size_t const     a5_new_i[] = { 1, 4, wlimits_type::digits + 3, wlimits_type::digits + 4, 2 * wlimits_type::digits + 5 };
    size_t const     a5_new_s = sizeof( a5_new_i ) / sizeof( a5_new_i[0] );
    va_size_t const  a5_new1( a5_new_i, a5_new_s );
    va_size_t const  a5_new2( a5_new_i, a5_new_s - 1 );

    a5.flip( 4 );
    BOOST_CHECK( equal_valarrays(a5_new1, a5.to_bit_indices()) );

    a5.flip( 2 * wlimits_type::digits + 5 );
    BOOST_CHECK( equal_valarrays(a5_new2, a5.to_bit_indices()) );
}

// Unit test for flipping a group of bits
void
bigwhole_group_bit_flip_unit_test
(
)
{
    using boost::math::big_whole;
    using std::size_t;

    typedef std::valarray<size_t>  va_size_t;
    typedef std::set<size_t>       st_size_t;

    // zero
    big_whole  a1;

    BOOST_CHECK_EQUAL( 0u, a1.to_uintmax() );

    a1.flip( 3, 7 );
    BOOST_CHECK_EQUAL( 248u, a1.to_uintmax() );

    a1.flip( 2, 4 );
    BOOST_CHECK_EQUAL( 228u, a1.to_uintmax() );

    // one bit set
    big_whole  a2( 8 );

    BOOST_CHECK_EQUAL( 8u, a2.to_uintmax() );

    a2.flip( 6, 9 );
    BOOST_CHECK_EQUAL( 968u, a2.to_uintmax() );

    a2.flip( 2, 3 );
    BOOST_CHECK_EQUAL( 964u, a2.to_uintmax() );

    // multiple bits set
    big_whole  a3( 25 );

    BOOST_CHECK_EQUAL( 25u, a3.to_uintmax() );

    a3.flip( 3, 3 );
    BOOST_CHECK_EQUAL( 17u, a3.to_uintmax() );

    a3.flip( 1, 2 );
    BOOST_CHECK_EQUAL( 23u, a3.to_uintmax() );

    a3.flip( 2, 6 );
    BOOST_CHECK_EQUAL( 107u, a3.to_uintmax() );

    // two words
    size_t const     a4_old_i[] = { 0, wlimits_type::digits + 1 };
    size_t const     a4_old_s = sizeof( a4_old_i ) / sizeof( a4_old_i[0] );
    va_size_t const  a4_old( a4_old_i, a4_old_s );
    big_whole        a4( a4_old );
    st_size_t        a4_new( a4_old_i, a4_old_i + a4_old_s );

    BOOST_CHECK( equal_valarrays(a4_old, a4.to_bit_indices()) );

    a4.flip( 5, 12 );
    insert_value_range( a4_new, 5, 12 );
    BOOST_CHECK( equal_valarrays(set_to_valarray( a4_new ), a4.to_bit_indices()) );

    a4.flip( 9, 2 * wlimits_type::digits );
    insert_value_range( a4_new, 9, 2 * wlimits_type::digits );
    erase_value_range( a4_new, 9, 12 );
    a4_new.erase( wlimits_type::digits + 1 );
    BOOST_CHECK( equal_valarrays(set_to_valarray( a4_new ), a4.to_bit_indices()) );

    // more-than-two words
    size_t const     a5_old_i[] = { 1, wlimits_type::digits + 3, wlimits_type::digits + 4, 2 * wlimits_type::digits + 5 };
    size_t const     a5_old_s = sizeof( a5_old_i ) / sizeof( a5_old_i[0] );
    va_size_t const  a5_old( a5_old_i, a5_old_s );
    big_whole        a5( a5_old );
    st_size_t        a5_new( a5_old_i, a5_old_i + a5_old_s );

    BOOST_CHECK( equal_valarrays(a5_old, a5.to_bit_indices()) );

    a5.flip( 3, 12 );
    insert_value_range( a5_new, 3, 12 );
    BOOST_CHECK( equal_valarrays(set_to_valarray( a5_new ), a5.to_bit_indices()) );

    a5.flip( wlimits_type::digits - 1, 2 * wlimits_type::digits + 1 );
    insert_value_range( a5_new, wlimits_type::digits - 1, 2 * wlimits_type::digits + 1 );
    a5_new.erase( wlimits_type::digits + 3 );
    a5_new.erase( wlimits_type::digits + 4 );
    BOOST_CHECK( equal_valarrays(set_to_valarray( a5_new ), a5.to_bit_indices()) );

    a5.flip( 2 * wlimits_type::digits + 1, 3 * wlimits_type::digits );
    insert_value_range( a5_new, 2 * wlimits_type::digits + 1, 3 * wlimits_type::digits );
    a5_new.erase( 2 * wlimits_type::digits + 1 );
    a5_new.erase( 2 * wlimits_type::digits + 5 );
    BOOST_CHECK( equal_valarrays(set_to_valarray( a5_new ), a5.to_bit_indices()) );

    a5.flip( wlimits_type::digits + 1, 3 * wlimits_type::digits );
    erase_value_range( a5_new, wlimits_type::digits + 1, 3 * wlimits_type::digits );
    a5_new.insert( wlimits_type::digits + 3 );
    a5_new.insert( wlimits_type::digits + 4 );
    a5_new.insert( 2 * wlimits_type::digits + 1 );
    a5_new.insert( 2 * wlimits_type::digits + 5 );
    BOOST_CHECK( equal_valarrays(set_to_valarray( a5_new ), a5.to_bit_indices()) );
}

// Unit test for assigning a group of bits to an arbitrary value
void
bigwhole_group_bit_assign_unit_test
(
)
{
    using boost::math::big_whole;
    using std::size_t;

    typedef std::valarray<size_t>  va_size_t;
    typedef std::set<size_t>       st_size_t;

    // no segment is zero
    big_whole  x( 255 );

    BOOST_CHECK_EQUAL( 255u, x.to_uintmax() );

    x.bits_assign( 3, 5, big_whole(2) );
    BOOST_CHECK_EQUAL( 215u, x.to_uintmax() );

    // old middle part is zero
    x.assign( 199 );
    BOOST_CHECK_EQUAL( 199u, x.to_uintmax() );

    x.bits_assign( 3, 5, big_whole(7) );
    BOOST_CHECK_EQUAL( 255u, x.to_uintmax() );

    // new-value bits above the range are ignored
    x.assign( 199 );
    BOOST_CHECK_EQUAL( 199u, x.to_uintmax() );

    x.bits_assign( 3, 4, big_whole(14) );
    BOOST_CHECK_EQUAL( 215u, x.to_uintmax() );

    // new middle part is zero
    x.assign( 255 );
    BOOST_CHECK_EQUAL( 255u, x.to_uintmax() );

    x.bits_assign( 3, 5, big_whole() );
    BOOST_CHECK_EQUAL( 199u, x.to_uintmax() );

    // change the lowest bits
    x.assign( 255 );
    BOOST_CHECK_EQUAL( 255u, x.to_uintmax() );

    x.bits_assign( 0, 2, big_whole(33) );
    BOOST_CHECK_EQUAL( 249u, x.to_uintmax() );

    // change the lowest bits to zero
    x.assign( 255 );
    BOOST_CHECK_EQUAL( 255u, x.to_uintmax() );

    x.bits_assign( 0, 3, big_whole() );
    BOOST_CHECK_EQUAL( 240u, x.to_uintmax() );

    // change the highest bits
    x.assign( 255 );
    BOOST_CHECK_EQUAL( 255u, x.to_uintmax() );

    x.bits_assign( 5, 7, big_whole(5) );
    BOOST_CHECK_EQUAL( 191u, x.to_uintmax() );

    // change the highest bits to zero
    x.assign( 255 );
    BOOST_CHECK_EQUAL( 255u, x.to_uintmax() );

    x.bits_assign( 4, 7, big_whole() );
    BOOST_CHECK_EQUAL( 15u, x.to_uintmax() );

    // turn zero into a nonzero
    x.reset();
    BOOST_CHECK_EQUAL( 0u, x.to_uintmax() );

    x.bits_assign( 8, 9, big_whole(35) );
    BOOST_CHECK_EQUAL( 768u, x.to_uintmax() );

    // keep zero a zero
    x.reset();
    BOOST_CHECK_EQUAL( 0u, x.to_uintmax() );

    x.bits_assign( 11, 45, big_whole() );
    BOOST_CHECK_EQUAL( 0u, x.to_uintmax() );

    // multiple words
    size_t const     t11_old_i[] = { 0,  2 * wlimits_type::digits + 1 };
    size_t const     t11_old_s = sizeof( t11_old_i ) / sizeof( t11_old_i[0] );
    va_size_t const  t11_old( t11_old_i, t11_old_s );
    st_size_t        t11_new( t11_old_i, t11_old_i + t11_old_s );

    x.reconfigure( t11_old );
    BOOST_CHECK( equal_valarrays(t11_old, x.to_bit_indices()) );
    BOOST_CHECK( equal_valarrays(set_to_valarray( t11_new ), t11_old) );

    x.bits_assign( wlimits_type::digits - 1, wlimits_type::digits + 2, big_whole(7) );
    insert_value_range( t11_new, wlimits_type::digits - 1, wlimits_type::digits + 1 );
    BOOST_CHECK( equal_valarrays(x.to_bit_indices(), set_to_valarray( t11_new )) );
}

// Unit test for checking if an object is even/odd
void
bigwhole_is_even_unit_test
(
)
{
    using boost::math::big_whole;
    using std::size_t;

    typedef std::valarray<size_t>  va_size_t;

    // zero
    big_whole  a;

    BOOST_CHECK_EQUAL( 0u, a.to_uintmax() );
    BOOST_CHECK( a.is_even() );

    // non-zero odd
    a.assign( 5 );
    BOOST_CHECK_EQUAL( 5u, a.to_uintmax() );
    BOOST_CHECK( !a.is_even() );

    // non-zero even
    a.assign( 18 );
    BOOST_CHECK_EQUAL( 18u, a.to_uintmax() );
    BOOST_CHECK( a.is_even() );

    // multi-word odd
    size_t const     t3_i[] = { 0, wlimits_type::digits + 1 };
    size_t const     t3_s = sizeof( t3_i ) / sizeof( t3_i[0] );
    va_size_t const  t3( t3_i, t3_s );

    a.reconfigure( t3 );
    BOOST_CHECK( equal_valarrays(t3, a.to_bit_indices()) );
    BOOST_CHECK( !a.is_even() );

    // multi-word even
    size_t const     t4_i = 2 * wlimits_type::digits + 3;
    va_size_t const  t4( t4_i, 1 );

    a.reconfigure( t4 );
    BOOST_CHECK( equal_valarrays(t4, a.to_bit_indices()) );
    BOOST_CHECK( a.is_even() );
}

// Unit test for comparisons
void
bigwhole_compare_unit_test
(
)
{
    using boost::math::big_whole;

    typedef std::valarray<std::size_t>  va_size_t;

    big_whole const  z;
    big_whole const  a( 7 );
    big_whole const  b( 100 );
    big_whole const  c( va_size_t(2 * wlimits_type::digits + 3, 1) );

    // "compare" directly
    BOOST_CHECK( z.compare(z) == 0 );
    BOOST_CHECK( z.compare(a) < 0 );
    BOOST_CHECK( z.compare(b) < 0 );
    BOOST_CHECK( z.compare(c) < 0 );

    BOOST_CHECK( a.compare(z) > 0 );
    BOOST_CHECK( a.compare(a) == 0 );
    BOOST_CHECK( a.compare(b) < 0 );
    BOOST_CHECK( a.compare(c) < 0 );

    BOOST_CHECK( b.compare(z) > 0 );
    BOOST_CHECK( b.compare(a) > 0 );
    BOOST_CHECK( b.compare(b) == 0 );
    BOOST_CHECK( b.compare(c) < 0 );

    BOOST_CHECK( c.compare(z) > 0 );
    BOOST_CHECK( c.compare(a) > 0 );
    BOOST_CHECK( c.compare(b) > 0 );
    BOOST_CHECK( c.compare(c) == 0 );

    // use ==
    BOOST_CHECK( z == z ); BOOST_CHECK( !(z == a) ); BOOST_CHECK( !(z == b) ); BOOST_CHECK( !(z == c) );
    BOOST_CHECK( !(a == z) ); BOOST_CHECK( a == a ); BOOST_CHECK( !(a == b) ); BOOST_CHECK( !(a == c) );
    BOOST_CHECK( !(b == z) ); BOOST_CHECK( !(b == a) ); BOOST_CHECK( b == b ); BOOST_CHECK( !(b == c) );
    BOOST_CHECK( !(c == z) ); BOOST_CHECK( !(c == a) ); BOOST_CHECK( !(c == b) ); BOOST_CHECK( c == c );

    // use !=
    BOOST_CHECK( !(z != z) ); BOOST_CHECK( z != a ); BOOST_CHECK( z != b ); BOOST_CHECK( z != c );
    BOOST_CHECK( a != z ); BOOST_CHECK( !(a != a) ); BOOST_CHECK( a != b ); BOOST_CHECK( a != c );
    BOOST_CHECK( b != z ); BOOST_CHECK( b != a ); BOOST_CHECK( !(b != b) ); BOOST_CHECK( b != c );
    BOOST_CHECK( c != z ); BOOST_CHECK( c != a ); BOOST_CHECK( c != b ); BOOST_CHECK( !(c != c) );

    // use <
    BOOST_CHECK( !(z < z) ); BOOST_CHECK( z < a ); BOOST_CHECK( z < b ); BOOST_CHECK( z < c );
    BOOST_CHECK( !(a < z) ); BOOST_CHECK( !(a < a) ); BOOST_CHECK( a < b ); BOOST_CHECK( a < c );
    BOOST_CHECK( !(b < z) ); BOOST_CHECK( !(b < a) ); BOOST_CHECK( !(b < b) ); BOOST_CHECK( b < c );
    BOOST_CHECK( !(c < z) ); BOOST_CHECK( !(c < a) ); BOOST_CHECK( !(c < b) ); BOOST_CHECK( !(c < c) );

    // use >
    BOOST_CHECK( !(z > z) ); BOOST_CHECK( !(z > a) ); BOOST_CHECK( !(z > b) ); BOOST_CHECK( !(z > c) );
    BOOST_CHECK( a > z ); BOOST_CHECK( !(a > a) ); BOOST_CHECK( !(a > b) ); BOOST_CHECK( !(a > c) );
    BOOST_CHECK( b > z ); BOOST_CHECK( b > a ); BOOST_CHECK( !(b > b) ); BOOST_CHECK( !(b > c) );
    BOOST_CHECK( c > z ); BOOST_CHECK( c > a ); BOOST_CHECK( c > b ); BOOST_CHECK( !(c > c) );

    // use <=
    BOOST_CHECK( z <= z ); BOOST_CHECK( z <= a ); BOOST_CHECK( z <= b ); BOOST_CHECK( z <= c );
    BOOST_CHECK( !(a <= z) ); BOOST_CHECK( a <= a ); BOOST_CHECK( a <= b ); BOOST_CHECK( a <= c );
    BOOST_CHECK( !(b <= z) ); BOOST_CHECK( !(b <= a) ); BOOST_CHECK( b <= b ); BOOST_CHECK( b <= c );
    BOOST_CHECK( !(c <= z) ); BOOST_CHECK( !(c <= a) ); BOOST_CHECK( !(c <= b) ); BOOST_CHECK( c <= c );

    // use >=
    BOOST_CHECK( z >= z ); BOOST_CHECK( !(z >= a) ); BOOST_CHECK( !(z >= b) ); BOOST_CHECK( !(z >= c) );
    BOOST_CHECK( a >= z ); BOOST_CHECK( a >= a ); BOOST_CHECK( !(a >= b) ); BOOST_CHECK( !(a >= c) );
    BOOST_CHECK( b >= z ); BOOST_CHECK( b >= a ); BOOST_CHECK( b >= b ); BOOST_CHECK( !(b >= c) );
    BOOST_CHECK( c >= z ); BOOST_CHECK( c >= a ); BOOST_CHECK( c >= b ); BOOST_CHECK( c >= c );
}

// Unit test for bitwise-and
void
bigwhole_bitwise_and_unit_test
(
)
{
    using boost::math::big_whole;
    using std::size_t;

    typedef std::valarray<size_t>  va_size_t;

    size_t const     di[] = { wlimits_type::digits + 1, 2 * wlimits_type::digits + 3 };
    size_t const     ds = sizeof( di ) / sizeof( di[0] );
    va_size_t const  dv( di, ds );
    size_t const     ei[] = { 1, wlimits_type::digits - 1, wlimits_type::digits, wlimits_type::digits + 1 };
    size_t const     es = sizeof( ei ) / sizeof( ei[0] );
    va_size_t const  ev( ei, es );

    big_whole const  z;
    big_whole const  a( 7 );
    big_whole const  b( 12 );
    big_whole const  c( 100 );
    big_whole const  d( dv );
    big_whole const  e( ev );

    // self-and
    BOOST_CHECK_EQUAL( z, z & z );
    BOOST_CHECK_EQUAL( a, a & a );
    BOOST_CHECK_EQUAL( b, b & b );
    BOOST_CHECK_EQUAL( c, c & c );
    BOOST_CHECK_EQUAL( d, d & d );
    BOOST_CHECK_EQUAL( e, e & e );

    // zero-and
    BOOST_CHECK_EQUAL( z, z & a ); BOOST_CHECK_EQUAL( z, a & z );
    BOOST_CHECK_EQUAL( z, z & b ); BOOST_CHECK_EQUAL( z, b & z );
    BOOST_CHECK_EQUAL( z, z & c ); BOOST_CHECK_EQUAL( z, c & z );
    BOOST_CHECK_EQUAL( z, z & d ); BOOST_CHECK_EQUAL( z, d & z );
    BOOST_CHECK_EQUAL( z, z & e ); BOOST_CHECK_EQUAL( z, e & z );

    // various combinations
    word_type const  ab = 4, ac = 4, ad = 0, ae = 2;
    word_type const  ba = ab, bc = 4, bd = 0, be = 0;
    word_type const  ca = ac, cb = bc, cd = 0, ce = 0;
    word_type const  da = ad, db = bd, dc = cd;
    word_type const  ea = ae, eb = be, ec = ce;
    big_whole const  de( va_size_t(wlimits_type::digits + 1, 1) ), ed( de );

    BOOST_CHECK_EQUAL( ab, a & b ); BOOST_CHECK_EQUAL( ba, b & a );
    BOOST_CHECK_EQUAL( ac, a & c ); BOOST_CHECK_EQUAL( ca, c & a );
    BOOST_CHECK_EQUAL( ad, a & d ); BOOST_CHECK_EQUAL( da, d & a );
    BOOST_CHECK_EQUAL( ae, a & e ); BOOST_CHECK_EQUAL( ea, e & a );

    BOOST_CHECK_EQUAL( bc, b & c ); BOOST_CHECK_EQUAL( cb, c & b );
    BOOST_CHECK_EQUAL( bd, b & d ); BOOST_CHECK_EQUAL( db, d & b );
    BOOST_CHECK_EQUAL( be, b & e ); BOOST_CHECK_EQUAL( eb, e & b );

    BOOST_CHECK_EQUAL( cd, c & d ); BOOST_CHECK_EQUAL( dc, d & c );
    BOOST_CHECK_EQUAL( ce, c & e ); BOOST_CHECK_EQUAL( ec, e & c );

    BOOST_CHECK_EQUAL( de, d & e ); BOOST_CHECK_EQUAL( ed, e & d );
}

// Unit test for bitwise-or
void
bigwhole_bitwise_or_unit_test
(
)
{
    using boost::math::big_whole;
    using std::size_t;

    typedef std::valarray<size_t>  va_size_t;

    size_t const     di[] = { wlimits_type::digits + 1, 2 * wlimits_type::digits + 3 };
    size_t const     ds = sizeof( di ) / sizeof( di[0] );
    va_size_t const  dv( di, ds );
    size_t const     ei[] = { 1, wlimits_type::digits - 1, wlimits_type::digits, wlimits_type::digits + 1 };
    size_t const     es = sizeof( ei ) / sizeof( ei[0] );
    va_size_t const  ev( ei, es );

    big_whole const  z;
    big_whole const  a( 7 );
    big_whole const  b( 12 );
    big_whole const  c( 100 );
    big_whole const  d( dv );
    big_whole const  e( ev );

    // self-or
    BOOST_CHECK_EQUAL( z, z | z );
    BOOST_CHECK_EQUAL( a, a | a );
    BOOST_CHECK_EQUAL( b, b | b );
    BOOST_CHECK_EQUAL( c, c | c );
    BOOST_CHECK_EQUAL( d, d | d );
    BOOST_CHECK_EQUAL( e, e | e );

    // zero-or
    BOOST_CHECK_EQUAL( a, z | a ); BOOST_CHECK_EQUAL( a, a | z );
    BOOST_CHECK_EQUAL( b, z | b ); BOOST_CHECK_EQUAL( b, b | z );
    BOOST_CHECK_EQUAL( c, z | c ); BOOST_CHECK_EQUAL( c, c | z );
    BOOST_CHECK_EQUAL( d, z | d ); BOOST_CHECK_EQUAL( d, d | z );
    BOOST_CHECK_EQUAL( e, z | e ); BOOST_CHECK_EQUAL( e, e | z );

    // various combinations
    size_t const  ed_i[] = { 1, wlimits_type::digits - 1, wlimits_type::digits, wlimits_type::digits + 1, 2 * wlimits_type::digits + 3 };
    size_t const  ed_s = sizeof( ed_i ) / sizeof( ed_i[0] );
    size_t const  ec_i[] = { 1, 2, 5, 6, wlimits_type::digits - 1, wlimits_type::digits, wlimits_type::digits + 1 };
    size_t const  ec_s = sizeof( ec_i ) / sizeof( ec_i[0] );
    size_t const  eb_i[] = { 1, 2, 3, wlimits_type::digits - 1, wlimits_type::digits, wlimits_type::digits + 1 };
    size_t const  eb_s = sizeof( eb_i ) / sizeof( eb_i[0] );
    size_t const  ea_i[] = { 0, 1, 2, wlimits_type::digits - 1, wlimits_type::digits, wlimits_type::digits + 1 };
    size_t const  ea_s = sizeof( ea_i ) / sizeof( ea_i[0] );
    size_t const  dc_i[] = { 2, 5, 6, wlimits_type::digits + 1, 2 * wlimits_type::digits + 3 };
    size_t const  dc_s = sizeof( dc_i ) / sizeof( dc_i[0] );
    size_t const  db_i[] = { 2, 3, wlimits_type::digits + 1, 2 * wlimits_type::digits + 3 };
    size_t const  db_s = sizeof( db_i ) / sizeof( db_i[0] );
    size_t const  da_i[] = { 0, 1, 2, wlimits_type::digits + 1, 2 * wlimits_type::digits + 3 };
    size_t const  da_s = sizeof( da_i ) / sizeof( da_i[0] );

    va_size_t const  ed_v( ed_i, ed_s ), ec_v( ec_i, ec_s ), eb_v( eb_i, eb_s );
    va_size_t const  ea_v( ea_i, ea_s ), dc_v( dc_i, dc_s ), db_v( db_i, db_s );
    va_size_t const  da_v( da_i, da_s );

    big_whole const  ed( ed_v ), ec( ec_v ), eb( eb_v ), ea( ea_v );
    big_whole const  de( ed ), dc( dc_v ), db( db_v ), da( da_v );
    big_whole const  ce( ec ), cd( dc ), be( eb ), bd( db ), ae( ea ), ad( da );
    word_type const  cb = 108, ca = 103, bc = cb, ba = 15, ac = ca, ab = ba;

    BOOST_CHECK_EQUAL( ab, a | b ); BOOST_CHECK_EQUAL( ba, b | a );
    BOOST_CHECK_EQUAL( ac, a | c ); BOOST_CHECK_EQUAL( ca, c | a );
    BOOST_CHECK_EQUAL( ad, a | d ); BOOST_CHECK_EQUAL( da, d | a );
    BOOST_CHECK_EQUAL( ae, a | e ); BOOST_CHECK_EQUAL( ea, e | a );

    BOOST_CHECK_EQUAL( bc, b | c ); BOOST_CHECK_EQUAL( cb, c | b );
    BOOST_CHECK_EQUAL( bd, b | d ); BOOST_CHECK_EQUAL( db, d | b );
    BOOST_CHECK_EQUAL( be, b | e ); BOOST_CHECK_EQUAL( eb, e | b );

    BOOST_CHECK_EQUAL( cd, c | d ); BOOST_CHECK_EQUAL( dc, d | c );
    BOOST_CHECK_EQUAL( ce, c | e ); BOOST_CHECK_EQUAL( ec, e | c );

    BOOST_CHECK_EQUAL( de, d | e ); BOOST_CHECK_EQUAL( ed, e | d );
}

// Unit test for bitwise-xor
void
bigwhole_bitwise_xor_unit_test
(
)
{
    using boost::math::big_whole;
    using std::size_t;

    typedef std::valarray<size_t>  va_size_t;

    size_t const     di[] = { wlimits_type::digits + 1, 2 * wlimits_type::digits + 3 };
    size_t const     ds = sizeof( di ) / sizeof( di[0] );
    va_size_t const  dv( di, ds );
    size_t const     ei[] = { 1, wlimits_type::digits - 1, wlimits_type::digits, wlimits_type::digits + 1 };
    size_t const     es = sizeof( ei ) / sizeof( ei[0] );
    va_size_t const  ev( ei, es );

    big_whole const  z;
    big_whole const  a( 7 );
    big_whole const  b( 12 );
    big_whole const  c( 100 );
    big_whole const  d( dv );
    big_whole const  e( ev );

    // self-or
    BOOST_CHECK_EQUAL( z, z ^ z );
    BOOST_CHECK_EQUAL( z, a ^ a );
    BOOST_CHECK_EQUAL( z, b ^ b );
    BOOST_CHECK_EQUAL( z, c ^ c );
    BOOST_CHECK_EQUAL( z, d ^ d );
    BOOST_CHECK_EQUAL( z, e ^ e );

    // zero-or
    BOOST_CHECK_EQUAL( a, z ^ a ); BOOST_CHECK_EQUAL( a, a ^ z );
    BOOST_CHECK_EQUAL( b, z ^ b ); BOOST_CHECK_EQUAL( b, b ^ z );
    BOOST_CHECK_EQUAL( c, z ^ c ); BOOST_CHECK_EQUAL( c, c ^ z );
    BOOST_CHECK_EQUAL( d, z ^ d ); BOOST_CHECK_EQUAL( d, d ^ z );
    BOOST_CHECK_EQUAL( e, z ^ e ); BOOST_CHECK_EQUAL( e, e ^ z );

    // various combinations
    size_t const  ed_i[] = { 1, wlimits_type::digits - 1, wlimits_type::digits, 2 * wlimits_type::digits + 3 };
    size_t const  ed_s = sizeof( ed_i ) / sizeof( ed_i[0] );
    size_t const  ec_i[] = { 1, 2, 5, 6, wlimits_type::digits - 1, wlimits_type::digits, wlimits_type::digits + 1 };
    size_t const  ec_s = sizeof( ec_i ) / sizeof( ec_i[0] );
    size_t const  eb_i[] = { 1, 2, 3, wlimits_type::digits - 1, wlimits_type::digits, wlimits_type::digits + 1 };
    size_t const  eb_s = sizeof( eb_i ) / sizeof( eb_i[0] );
    size_t const  ea_i[] = { 0, 2, wlimits_type::digits - 1, wlimits_type::digits, wlimits_type::digits + 1 };
    size_t const  ea_s = sizeof( ea_i ) / sizeof( ea_i[0] );
    size_t const  dc_i[] = { 2, 5, 6, wlimits_type::digits + 1, 2 * wlimits_type::digits + 3 };
    size_t const  dc_s = sizeof( dc_i ) / sizeof( dc_i[0] );
    size_t const  db_i[] = { 2, 3, wlimits_type::digits + 1, 2 * wlimits_type::digits + 3 };
    size_t const  db_s = sizeof( db_i ) / sizeof( db_i[0] );
    size_t const  da_i[] = { 0, 1, 2, wlimits_type::digits + 1, 2 * wlimits_type::digits + 3 };
    size_t const  da_s = sizeof( da_i ) / sizeof( da_i[0] );

    va_size_t const  ed_v( ed_i, ed_s ), ec_v( ec_i, ec_s ), eb_v( eb_i, eb_s );
    va_size_t const  ea_v( ea_i, ea_s ), dc_v( dc_i, dc_s ), db_v( db_i, db_s );
    va_size_t const  da_v( da_i, da_s );

    big_whole const  ed( ed_v ), ec( ec_v ), eb( eb_v ), ea( ea_v );
    big_whole const  de( ed ), dc( dc_v ), db( db_v ), da( da_v );
    big_whole const  ce( ec ), cd( dc ), be( eb ), bd( db ), ae( ea ), ad( da );
    word_type const  cb = 104, ca = 99, bc = cb, ba = 11, ac = ca, ab = ba;

    BOOST_CHECK_EQUAL( ab, a ^ b ); BOOST_CHECK_EQUAL( ba, b ^ a );
    BOOST_CHECK_EQUAL( ac, a ^ c ); BOOST_CHECK_EQUAL( ca, c ^ a );
    BOOST_CHECK_EQUAL( ad, a ^ d ); BOOST_CHECK_EQUAL( da, d ^ a );
    BOOST_CHECK_EQUAL( ae, a ^ e ); BOOST_CHECK_EQUAL( ea, e ^ a );

    BOOST_CHECK_EQUAL( bc, b ^ c ); BOOST_CHECK_EQUAL( cb, c ^ b );
    BOOST_CHECK_EQUAL( bd, b ^ d ); BOOST_CHECK_EQUAL( db, d ^ b );
    BOOST_CHECK_EQUAL( be, b ^ e ); BOOST_CHECK_EQUAL( eb, e ^ b );

    BOOST_CHECK_EQUAL( cd, c ^ d ); BOOST_CHECK_EQUAL( dc, d ^ c );
    BOOST_CHECK_EQUAL( ce, c ^ e ); BOOST_CHECK_EQUAL( ec, e ^ c );

    BOOST_CHECK_EQUAL( de, d ^ e ); BOOST_CHECK_EQUAL( ed, e ^ d );
}

// Unit test for left-shift (value-decreasing)
void
bigwhole_left_shift_unit_test
(
)
{
    using boost::math::big_whole;
    using std::size_t;
    using std::length_error;

    typedef std::valarray<size_t>  va_size_t;

    typedef std::numeric_limits<size_t>  size_limits;

    size_t const     ei[] = { wlimits_type::digits + 1, 2 * wlimits_type::digits + 3 };
    size_t const     es = sizeof( ei ) / sizeof( ei[0] );
    va_size_t const  ev( ei, es );
    size_t const     e1i[] = { 1, wlimits_type::digits + 3 };
    size_t const     e1s = sizeof( e1i ) / sizeof( e1i[0] );
    va_size_t const  e1v( e1i, e1s );
    size_t const     fi[] = { 1, wlimits_type::digits - 1, wlimits_type::digits, wlimits_type::digits + 2 };
    size_t const     fs = sizeof( fi ) / sizeof( fi[0] );
    va_size_t const  fv( fi, fs );

    big_whole const  z;
    big_whole const  a( 1 );
    big_whole const  b( 5 );
    big_whole const  c( 10 );
    big_whole const  d( 100 );
    big_whole const  e( ev ), e1( e1v );
    big_whole const  f( fv );

    // shifting zero
    BOOST_CHECK( !(z >> z) );
    BOOST_CHECK( !(z >> a) );
    BOOST_CHECK( !(z >> d) );

    // zero-shift
    BOOST_CHECK_EQUAL( z, z >> z );
    BOOST_CHECK_EQUAL( a, a >> z );
    BOOST_CHECK_EQUAL( b, b >> z );
    BOOST_CHECK_EQUAL( c, c >> z );
    BOOST_CHECK_EQUAL( d, d >> z );
    BOOST_CHECK_EQUAL( e, e >> z );
    BOOST_CHECK_EQUAL( f, f >> z );

    // single-word shifts
    BOOST_CHECK_EQUAL( z, a >> a );

    BOOST_CHECK_EQUAL( z, b >> 3 );
    BOOST_CHECK_EQUAL( a, b >> 2 );
    BOOST_CHECK_EQUAL( 2, b >> 1 );

    BOOST_CHECK_EQUAL( b, c >> 1 );
    BOOST_CHECK_EQUAL( 2, c >> 2 );
    BOOST_CHECK_EQUAL( a, c >> 3 );
    BOOST_CHECK_EQUAL( z, c >> 4 );

    BOOST_CHECK_EQUAL( 50, d >> 1 );
    BOOST_CHECK_EQUAL( 25, d >> 2 );
    BOOST_CHECK_EQUAL( 12, d >> 3 );
    BOOST_CHECK_EQUAL( 6, d >> 4 );
    BOOST_CHECK_EQUAL( 3, d >> 5 );
    BOOST_CHECK_EQUAL( a, d >> 6 );
    BOOST_CHECK_EQUAL( z, d >> 7 );

    // over shifts
    BOOST_CHECK_EQUAL( z, a >> 4 );
    BOOST_CHECK_EQUAL( z, b >> 10 );
    BOOST_CHECK_EQUAL( z, c >> 25 );
    BOOST_CHECK_EQUAL( z, d >> 1000 );
    BOOST_CHECK_EQUAL( z, e >> (3 * wlimits_type::digits) );
    BOOST_CHECK_EQUAL( z, f >> (2 * wlimits_type::digits) );

    // multi-word shifts
    BOOST_CHECK_EQUAL( e1, e >> wlimits_type::digits );
    BOOST_CHECK_EQUAL( 8, e >> (2 * wlimits_type::digits) );
    BOOST_CHECK_EQUAL( 4, e1 >> (wlimits_type::digits + 1) );
    BOOST_CHECK_EQUAL( b, f >> wlimits_type::digits );
    BOOST_CHECK_EQUAL( 11, f >> (wlimits_type::digits - 1) );

    // exceptions (hope we don't get any out-of-memory exceptions)
    va_size_t const  si( size_limits::digits, 1 );
    big_whole const  s( si );
    big_whole const  t( size_limits::max() );  // == s - 1
    big_whole const  u = s | 1;  // == s + 1

    BOOST_CHECK_NO_THROW( z >> s );
    BOOST_CHECK_THROW( a >> s, length_error );
    BOOST_CHECK_THROW( f >> s, length_error );
    BOOST_CHECK_NO_THROW( z >> t );
    BOOST_CHECK_NO_THROW( a >> t );
    BOOST_CHECK_NO_THROW( f >> t );
    BOOST_CHECK_NO_THROW( z >> u );
    BOOST_CHECK_THROW( a >> u, length_error );
    BOOST_CHECK_THROW( f >> u, length_error );
}

// Unit test for right-shift (value-increasing)
void
bigwhole_right_shift_unit_test
(
)
{
    using boost::math::big_whole;
    using std::size_t;
    using std::length_error;
    using std::overflow_error;

    typedef std::valarray<size_t>  va_size_t;

    typedef std::numeric_limits<size_t>  size_limits;

    size_t const     ei[] = { wlimits_type::digits + 1, 2 * wlimits_type::digits + 3 };
    size_t const     es = sizeof( ei ) / sizeof( ei[0] );
    va_size_t const  ev( ei, es );
    size_t const     e1i[] = { 1, wlimits_type::digits + 3 };
    size_t const     e1s = sizeof( e1i ) / sizeof( e1i[0] );
    va_size_t const  e1v( e1i, e1s );
    size_t const     fi[] = { 1, wlimits_type::digits - 1, wlimits_type::digits, wlimits_type::digits + 2 };
    size_t const     fs = sizeof( fi ) / sizeof( fi[0] );
    va_size_t const  fv( fi, fs );

    big_whole const  z;
    big_whole const  a( 1 );
    big_whole const  b( 5 );
    big_whole const  c( 10 );
    big_whole const  d( 100 );
    big_whole const  e( ev ), e1( e1v );
    big_whole const  f( fv );

    // shifting zero
    BOOST_CHECK( !(z << z) );
    BOOST_CHECK( !(z << a) );
    BOOST_CHECK( !(z << d) );

    // zero-shift
    BOOST_CHECK_EQUAL( z, z << z );
    BOOST_CHECK_EQUAL( a, a << z );
    BOOST_CHECK_EQUAL( b, b << z );
    BOOST_CHECK_EQUAL( c, c << z );
    BOOST_CHECK_EQUAL( d, d << z );
    BOOST_CHECK_EQUAL( e, e << z );
    BOOST_CHECK_EQUAL( f, f << z );

    // single-word shifts
    BOOST_CHECK_EQUAL( 2, a << a );
    BOOST_CHECK_EQUAL( 32, a << b );

    BOOST_CHECK_EQUAL( c, b << a );
    BOOST_CHECK_EQUAL( 20, b << 2 );
    BOOST_CHECK_EQUAL( 160, b << b );

    BOOST_CHECK_EQUAL( 20, c << a );
    BOOST_CHECK_EQUAL( 320, c << b );
    BOOST_CHECK_EQUAL( 10240, c << c );

    BOOST_CHECK_EQUAL( 200, d << a );
    BOOST_CHECK_EQUAL( 3200, d << b );
    BOOST_CHECK_EQUAL( 102400ul, d << c );

    // multi-word shifts
    BOOST_CHECK_EQUAL( e, e1 << wlimits_type::digits );
    BOOST_CHECK_EQUAL( 2 | (( a | c ) << ( wlimits_type::digits - 1 )), f );

    // exceptions (hope we don't get any out-of-memory exceptions)
    va_size_t const  si( size_limits::digits, 1 );
    big_whole const  s( si );
    big_whole const  t( size_limits::max() );
    big_whole const  u = s | 1;

    BOOST_CHECK_NO_THROW( z << s );
    BOOST_CHECK_THROW( a << s, length_error );
    BOOST_CHECK_THROW( f << s, length_error );

    #if PRIVATE_USE_BIG_MEMORY
    // This needs a lot of memory, so it shouldn't be regularly
    // tested.  Worse, what happens if a std::bad_alloc occurs?
    BOOST_CHECK_NO_THROW( a << t );
    #endif

    BOOST_CHECK_THROW( 2 << t, overflow_error );
    BOOST_CHECK_THROW( f << t, overflow_error );

    BOOST_CHECK_NO_THROW( z << u );
    BOOST_CHECK_THROW( a << u, length_error );
    BOOST_CHECK_THROW( f << u, length_error );
}

// Unit test for bitwise-and, bitwise-complementing the second value
void
bigwhole_bitwise_and_not_unit_test
(
)
{
    using boost::math::big_whole;
    using std::size_t;
    using boost::math::and_not;

    typedef std::valarray<size_t>  va_size_t;

    size_t const     di[] = { wlimits_type::digits + 1, 2 * wlimits_type::digits + 3 };
    size_t const     ds = sizeof( di ) / sizeof( di[0] );
    va_size_t const  dv( di, ds );
    size_t const     ei[] = { 1, wlimits_type::digits - 1, wlimits_type::digits, wlimits_type::digits + 1 };
    size_t const     es = sizeof( ei ) / sizeof( ei[0] );
    va_size_t const  ev( ei, es );

    big_whole const  z;
    big_whole const  a( 7 );
    big_whole const  b( 12 );
    big_whole const  c( 100 );
    big_whole const  d( dv );
    big_whole const  e( ev );

    // self-and
    BOOST_CHECK_EQUAL( z, and_not(z, z) );
    BOOST_CHECK_EQUAL( z, and_not(a, a) );
    BOOST_CHECK_EQUAL( z, and_not(b, b) );
    BOOST_CHECK_EQUAL( z, and_not(c, c) );
    BOOST_CHECK_EQUAL( z, and_not(d, d) );
    BOOST_CHECK_EQUAL( z, and_not(e, e) );

    // zero-and
    BOOST_CHECK_EQUAL( z, and_not(z, a) ); BOOST_CHECK_EQUAL( a, and_not(a, z) );
    BOOST_CHECK_EQUAL( z, and_not(z, b) ); BOOST_CHECK_EQUAL( b, and_not(b, z) );
    BOOST_CHECK_EQUAL( z, and_not(z, c) ); BOOST_CHECK_EQUAL( c, and_not(c, z) );
    BOOST_CHECK_EQUAL( z, and_not(z, d) ); BOOST_CHECK_EQUAL( d, and_not(d, z) );
    BOOST_CHECK_EQUAL( z, and_not(z, e) ); BOOST_CHECK_EQUAL( e, and_not(e, z) );

    // various combinations
    big_whole const  ab = 3, ac = 3, ad = a, ae = 5;
    big_whole const  ba = 8, bc = 8, bd = b, be = b;
    big_whole const  ca = 96, cb = 96, cd = c, ce = c;
    big_whole const  da = d, db = d, dc = d, de( va_size_t(di + 1, ds - 1) );
    big_whole const  ea( va_size_t(ei + 1, es - 1) ), eb = e, ec = e, ed( va_size_t(ei, es - 1) );

    BOOST_CHECK_EQUAL( ab, and_not(a, b) ); BOOST_CHECK_EQUAL( ba, and_not(b, a) );
    BOOST_CHECK_EQUAL( ac, and_not(a, c) ); BOOST_CHECK_EQUAL( ca, and_not(c, a) );
    BOOST_CHECK_EQUAL( ad, and_not(a, d) ); BOOST_CHECK_EQUAL( da, and_not(d, a) );
    BOOST_CHECK_EQUAL( ae, and_not(a, e) ); BOOST_CHECK_EQUAL( ea, and_not(e, a) );

    BOOST_CHECK_EQUAL( bc, and_not(b, c) ); BOOST_CHECK_EQUAL( cb, and_not(c, b) );
    BOOST_CHECK_EQUAL( bd, and_not(b, d) ); BOOST_CHECK_EQUAL( db, and_not(d, b) );
    BOOST_CHECK_EQUAL( be, and_not(b, e) ); BOOST_CHECK_EQUAL( eb, and_not(e, b) );

    BOOST_CHECK_EQUAL( cd, and_not(c, d) ); BOOST_CHECK_EQUAL( dc, and_not(d, c) );
    BOOST_CHECK_EQUAL( ce, and_not(c, e) ); BOOST_CHECK_EQUAL( ec, and_not(e, c) );

    BOOST_CHECK_EQUAL( de, and_not(d, e) ); BOOST_CHECK_EQUAL( ed, and_not(e, d) );
}

// Unit test for unary + and - operators
void
bigwhole_unary_plus_minus_unit_test
(
)
{
    using boost::math::big_whole;
    using std::range_error;

    big_whole const  z;
    big_whole const  a = !z;
    big_whole const  b = a << ( wlimits_type::digits + 3 );
    big_whole const  c = a | b;

    // plus
    BOOST_CHECK_EQUAL( z, +z );
    BOOST_CHECK_EQUAL( a, +a );
    BOOST_CHECK_EQUAL( b, +b );
    BOOST_CHECK_EQUAL( c, +c );

    // minus (hope we don't get any out-of-memory exceptions)
    BOOST_CHECK_NO_THROW( -z );
    BOOST_CHECK_EQUAL( z, -z );

    BOOST_CHECK_THROW( -a, range_error );
    BOOST_CHECK_THROW( -b, range_error );
    BOOST_CHECK_THROW( -c, range_error );
}

// Unit test for ++ and -- operators
void
bigwhole_double_plus_minus_unit_test
(
)
{
    using boost::math::big_whole;
    using std::range_error;
    using boost::math::and_not;

    typedef std::valarray<std::size_t>  va_size_t;

    big_whole const  z = 0;
    big_whole const  aa = wlimits_type::max();
    big_whole const  bb( va_size_t(wlimits_type::digits, 1) );

    // pre-increment
    big_whole  a = z;
    big_whole  b = aa;

    BOOST_CHECK_EQUAL( 1, ++a );
    BOOST_CHECK_EQUAL( 2, ++a );

    BOOST_CHECK_EQUAL( bb, ++b );
    BOOST_CHECK_EQUAL( bb | 1, ++b );

    // post-increment
    a.assign( 0 );
    b.assign( aa );

    BOOST_CHECK_EQUAL( 0, a++ );
    BOOST_CHECK_EQUAL( 1, a++ );
    BOOST_CHECK_EQUAL( 2, a );

    BOOST_CHECK_EQUAL( aa, b++ );
    BOOST_CHECK_EQUAL( bb, b++ );
    BOOST_CHECK_EQUAL( bb | 1, b );

    // decrement errors
    a.assign( 0 );

    BOOST_CHECK_THROW( --a, range_error );
    BOOST_CHECK( !a );
    BOOST_CHECK_THROW( a--, range_error );
    BOOST_CHECK( !a );

    // pre-decrement
    a.assign( 4 );
    b.assign( bb | 1 );

    BOOST_CHECK_EQUAL( 3, --a );
    BOOST_CHECK_EQUAL( 2, --a );
    BOOST_CHECK_EQUAL( 1, --a );
    BOOST_CHECK_EQUAL( 0, --a );

    BOOST_CHECK_EQUAL( bb, --b );
    BOOST_CHECK_EQUAL( aa, --b );
    BOOST_CHECK_EQUAL( and_not(aa, 1), --b );

    // post-decrement
    a.assign( 4 );
    b.assign( bb | 1 );

    BOOST_CHECK_EQUAL( 4, a-- );
    BOOST_CHECK_EQUAL( 3, a-- );
    BOOST_CHECK_EQUAL( 2, a-- );
    BOOST_CHECK_EQUAL( 1, a-- );
    BOOST_CHECK_EQUAL( 0, a );

    BOOST_CHECK_EQUAL( bb | 1, b-- );
    BOOST_CHECK_EQUAL( bb, b-- );
    BOOST_CHECK_EQUAL( aa, b-- );
    BOOST_CHECK_EQUAL( and_not(aa, 1), b );
}

// Unit test for the abs and sgn functions
void
bigwhole_abs_sgn_unit_test
(
)
{
    using boost::math::big_whole;

    typedef std::valarray<std::size_t>  va_size_t;

    big_whole const  z = 0;
    big_whole const  o = 1;
    big_whole const  a = wlimits_type::max();
    big_whole const  b( va_size_t(wlimits_type::digits, 1) );
    big_whole const  c = b | o;
    big_whole const  d = c << ( 2 * wlimits_type::digits + 3 );

    BOOST_CHECK_EQUAL( z, abs(z) );
    BOOST_CHECK_EQUAL( 0, sgn(z) );

    BOOST_CHECK_EQUAL( o, abs(o) );
    BOOST_CHECK_EQUAL( +1, sgn(o) );

    BOOST_CHECK_EQUAL( a, abs(a) );
    BOOST_CHECK_EQUAL( +1, sgn(a) );

    BOOST_CHECK_EQUAL( b, abs(b) );
    BOOST_CHECK_EQUAL( +1, sgn(b) );

    BOOST_CHECK_EQUAL( c, abs(c) );
    BOOST_CHECK_EQUAL( +1, sgn(c) );

    BOOST_CHECK_EQUAL( d, abs(d) );
    BOOST_CHECK_EQUAL( +1, sgn(d) );
}

// Unit test for binary + and - operators
void
bigwhole_binary_plus_minus_unit_test
(
)
{
    using boost::math::big_whole;
    using std::range_error;

    typedef std::valarray<std::size_t>  va_size_t;

    int const        wd = wlimits_type::digits;
    big_whole const  z;
    big_whole const  a = !z;
    big_whole const  b = wlimits_type::max();
    big_whole const  c( va_size_t(wd, 1) );
    big_whole const  d = a | c;
    big_whole const  e = a << ( wd + 3 );
    big_whole const  f = d << ( 2 * wd + 3 );

    // same-add
    BOOST_CHECK_EQUAL( z, z + z );
    BOOST_CHECK_EQUAL( a << 1, a + a );
    BOOST_CHECK_EQUAL( b << 1, b + b );
    BOOST_CHECK_EQUAL( c << 1, c + c );
    BOOST_CHECK_EQUAL( d << 1, d + d );
    BOOST_CHECK_EQUAL( e << 1, e + e );
    BOOST_CHECK_EQUAL( f << 1, f + f );

    // same-subtract
    BOOST_CHECK_EQUAL( z, z - z );
    BOOST_CHECK_EQUAL( z, a - a );
    BOOST_CHECK_EQUAL( z, b - b );
    BOOST_CHECK_EQUAL( z, c - c );
    BOOST_CHECK_EQUAL( z, d - d );
    BOOST_CHECK_EQUAL( z, e - e );
    BOOST_CHECK_EQUAL( z, f - f );

    // zero-add
    BOOST_CHECK_EQUAL( a, z + a );  BOOST_CHECK_EQUAL( a, a + z );
    BOOST_CHECK_EQUAL( b, z + b );  BOOST_CHECK_EQUAL( b, b + z );
    BOOST_CHECK_EQUAL( c, z + c );  BOOST_CHECK_EQUAL( c, c + z );
    BOOST_CHECK_EQUAL( d, z + d );  BOOST_CHECK_EQUAL( d, d + z );
    BOOST_CHECK_EQUAL( e, z + e );  BOOST_CHECK_EQUAL( e, e + z );
    BOOST_CHECK_EQUAL( f, z + f );  BOOST_CHECK_EQUAL( f, f + z );

    // zero-subtract
    BOOST_CHECK_THROW( z - a, range_error );  BOOST_CHECK_EQUAL( a, a - z );
    BOOST_CHECK_THROW( z - b, range_error );  BOOST_CHECK_EQUAL( b, b - z );
    BOOST_CHECK_THROW( z - c, range_error );  BOOST_CHECK_EQUAL( c, c - z );
    BOOST_CHECK_THROW( z - d, range_error );  BOOST_CHECK_EQUAL( d, d - z );
    BOOST_CHECK_THROW( z - e, range_error );  BOOST_CHECK_EQUAL( e, e - z );
    BOOST_CHECK_THROW( z - f, range_error );  BOOST_CHECK_EQUAL( f, f - z );

    // more bad subtractions
    BOOST_CHECK_THROW( a - b, range_error );  BOOST_CHECK_THROW( a - c, range_error );
    BOOST_CHECK_THROW( a - d, range_error );  BOOST_CHECK_THROW( a - e, range_error );
    BOOST_CHECK_THROW( a - f, range_error );  BOOST_CHECK_THROW( b - c, range_error );
    BOOST_CHECK_THROW( b - d, range_error );  BOOST_CHECK_THROW( b - e, range_error );
    BOOST_CHECK_THROW( b - f, range_error );  BOOST_CHECK_THROW( c - d, range_error );
    BOOST_CHECK_THROW( c - e, range_error );  BOOST_CHECK_THROW( c - f, range_error );
    BOOST_CHECK_THROW( d - e, range_error );  BOOST_CHECK_THROW( d - f, range_error );
    BOOST_CHECK_THROW( e - f, range_error );

    // additions
    big_whole const  ab1 = c, ac1 = d, ad1 = 3 ^ d, ae1 = a | e, af1 = a | f;
    big_whole const  ba1 = ab1, bc1 = b | c, bd1 = c << 1, be1 = b | e, bf1 = b | f;
    big_whole const  ca1 = ac1, cb1 = bc1, cd1 = 3 ^ (d << 1), ce1 = c | e, cf1 = c | f;
    big_whole const  da1 = ad1, db1 = bd1, dc1 = cd1, de1 = d | e, df1 = d | f;
    big_whole const  ea1 = ae1, eb1 = be1, ec1 = ce1, ed1 = de1, ef1 = e | f;
    big_whole const  fa1 = af1, fb1 = bf1, fc1 = cf1, fd1 = df1, fe1 = ef1;

    BOOST_CHECK_EQUAL( ab1, a + b );  BOOST_CHECK_EQUAL( ba1, b + a );
    BOOST_CHECK_EQUAL( ac1, a + c );  BOOST_CHECK_EQUAL( ca1, c + a );
    BOOST_CHECK_EQUAL( ad1, a + d );  BOOST_CHECK_EQUAL( da1, d + a );
    BOOST_CHECK_EQUAL( ae1, a + e );  BOOST_CHECK_EQUAL( ea1, e + a );
    BOOST_CHECK_EQUAL( af1, a + f );  BOOST_CHECK_EQUAL( fa1, f + a );

    BOOST_CHECK_EQUAL( bc1, b + c );  BOOST_CHECK_EQUAL( cb1, c + b );
    BOOST_CHECK_EQUAL( bd1, b + d );  BOOST_CHECK_EQUAL( db1, d + b );
    BOOST_CHECK_EQUAL( be1, b + e );  BOOST_CHECK_EQUAL( eb1, e + b );
    BOOST_CHECK_EQUAL( bf1, b + f );  BOOST_CHECK_EQUAL( fb1, f + b );

    BOOST_CHECK_EQUAL( cd1, c + d );  BOOST_CHECK_EQUAL( dc1, d + c );
    BOOST_CHECK_EQUAL( ce1, c + e );  BOOST_CHECK_EQUAL( ec1, e + c );
    BOOST_CHECK_EQUAL( cf1, c + f );  BOOST_CHECK_EQUAL( fc1, f + c );

    BOOST_CHECK_EQUAL( de1, d + e );  BOOST_CHECK_EQUAL( ed1, e + d );
    BOOST_CHECK_EQUAL( df1, d + f );  BOOST_CHECK_EQUAL( fd1, f + d );

    BOOST_CHECK_EQUAL( ef1, e + f );  BOOST_CHECK_EQUAL( fe1, f + e );

    // subtractions
    big_whole const  six = 6, seven = 7, eight = 8, one = a, maxd = b;

    big_whole const  ba2 = a ^ b;
    big_whole const  cb2 = a, ca2 = b;
    big_whole const  dc2 = a, db2 = a << 1, da2 = c;
    big_whole const  ed2 = b | ( six << wd );
    big_whole const  ec2 = seven << wd, eb2 = ec2 | a;
    big_whole const  ea2 = b | ( seven << wd );
    big_whole const  fe2 = ( (( (eight << wd) | seven ) << wd) | (maxd ^ seven) ) << wd;
    big_whole const  fd2 = ( (( (( eight << wd ) | seven) << wd ) | (maxd ^ one)) << wd ) | maxd;
    big_whole const  fc2 = ( (( (eight << wd) | seven ) << wd) | maxd ) << wd;
    big_whole const  fb2 = ( (( (( eight << wd ) | seven) << wd ) | maxd) << wd ) | one;
    big_whole const  fa2 = ( (( (( eight << wd ) | seven) << wd ) | maxd) << wd ) | maxd;

    BOOST_CHECK_EQUAL( ba2, b - a );  BOOST_CHECK_EQUAL( ca2, c - a );
    BOOST_CHECK_EQUAL( cb2, c - b );  BOOST_CHECK_EQUAL( da2, d - a );
    BOOST_CHECK_EQUAL( db2, d - b );  BOOST_CHECK_EQUAL( dc2, d - c );
    BOOST_CHECK_EQUAL( ea2, e - a );  BOOST_CHECK_EQUAL( eb2, e - b );
    BOOST_CHECK_EQUAL( ec2, e - c );  BOOST_CHECK_EQUAL( ed2, e - d );
    BOOST_CHECK_EQUAL( fa2, f - a );  BOOST_CHECK_EQUAL( fb2, f - b ); 
    BOOST_CHECK_EQUAL( fc2, f - c );  BOOST_CHECK_EQUAL( fd2, f - d );
    BOOST_CHECK_EQUAL( fe2, f - e );
}

// Unit test for intersects
void
bigwhole_intersects_unit_test
(
)
{
    using boost::math::big_whole;
    using std::size_t;

    typedef std::valarray<size_t>  va_size_t;

    size_t const     di[] = { wlimits_type::digits + 1, 2 * wlimits_type::digits + 3 };
    size_t const     ds = sizeof( di ) / sizeof( di[0] );
    va_size_t const  dv( di, ds );
    size_t const     ei[] = { 1, wlimits_type::digits - 1, wlimits_type::digits, wlimits_type::digits + 1 };
    size_t const     es = sizeof( ei ) / sizeof( ei[0] );
    va_size_t const  ev( ei, es );

    big_whole const  z;
    big_whole const  a( 7 );
    big_whole const  b( 12 );
    big_whole const  c( 100 );
    big_whole const  d( dv );
    big_whole const  e( ev );

    // self-intersects
    BOOST_CHECK( !z.intersects(z) );
    BOOST_CHECK(  a.intersects(a) );
    BOOST_CHECK(  b.intersects(b) );
    BOOST_CHECK(  c.intersects(c) );
    BOOST_CHECK(  d.intersects(d) );
    BOOST_CHECK(  e.intersects(e) );

    // zero-intersects
    BOOST_CHECK( !z.intersects(a) ); BOOST_CHECK( !a.intersects(z) );
    BOOST_CHECK( !z.intersects(b) ); BOOST_CHECK( !b.intersects(z) );
    BOOST_CHECK( !z.intersects(c) ); BOOST_CHECK( !c.intersects(z) );
    BOOST_CHECK( !z.intersects(d) ); BOOST_CHECK( !d.intersects(z) );
    BOOST_CHECK( !z.intersects(e) ); BOOST_CHECK( !e.intersects(z) );

    // various combinations
    bool const  ab = true, ac = true, ad = false, ae = true;
    bool const  ba = ab, bc = true, bd = false, be = false;
    bool const  ca = ac, cb = bc, cd = false, ce = false;
    bool const  da = ad, db = bd, dc = cd, de = true;
    bool const  ea = ae, eb = be, ec = ce, ed = de;

    BOOST_CHECK_EQUAL( ab, a.intersects(b) ); BOOST_CHECK_EQUAL( ba, b.intersects(a) );
    BOOST_CHECK_EQUAL( ac, a.intersects(c) ); BOOST_CHECK_EQUAL( ca, c.intersects(a) );
    BOOST_CHECK_EQUAL( ad, a.intersects(d) ); BOOST_CHECK_EQUAL( da, d.intersects(a) );
    BOOST_CHECK_EQUAL( ae, a.intersects(e) ); BOOST_CHECK_EQUAL( ea, e.intersects(a) );

    BOOST_CHECK_EQUAL( bc, b.intersects(c) ); BOOST_CHECK_EQUAL( cb, c.intersects(b) );
    BOOST_CHECK_EQUAL( bd, b.intersects(d) ); BOOST_CHECK_EQUAL( db, d.intersects(b) );
    BOOST_CHECK_EQUAL( be, b.intersects(e) ); BOOST_CHECK_EQUAL( eb, e.intersects(b) );

    BOOST_CHECK_EQUAL( cd, c.intersects(d) ); BOOST_CHECK_EQUAL( dc, d.intersects(c) );
    BOOST_CHECK_EQUAL( ce, c.intersects(e) ); BOOST_CHECK_EQUAL( ec, e.intersects(c) );

    BOOST_CHECK_EQUAL( de, d.intersects(e) ); BOOST_CHECK_EQUAL( ed, e.intersects(d) );
}

// Unit test for bit searching
void
bigwhole_bit_search_unit_test
(
)
{
    using boost::math::big_whole;
    using std::size_t;

    typedef std::valarray<size_t>  va_size_t;

    typedef std::numeric_limits<size_t>  size_limits;

    size_t const     wd = wlimits_type::digits;
    size_t const     ei[] = { wd + 1, 2 * wd + 3 };
    size_t const     es = sizeof( ei ) / sizeof( ei[0] );
    va_size_t const  ev( ei, es );

    big_whole const  z;
    big_whole const  a( 1 );
    big_whole const  b( 5 );
    big_whole const  e( ev );

    // limit violations
    BOOST_CHECK_EQUAL( size_limits::min(), z.next_set_bit(size_limits::max()) );
    BOOST_CHECK_EQUAL( size_limits::min(), a.next_set_bit(size_limits::max()) );
    BOOST_CHECK_EQUAL( size_limits::min(), b.next_set_bit(size_limits::max()) );
    BOOST_CHECK_EQUAL( size_limits::min(), e.next_set_bit(size_limits::max()) );

    BOOST_CHECK_EQUAL( size_limits::min(), z.next_reset_bit(size_limits::max()) );
    BOOST_CHECK_EQUAL( size_limits::min(), a.next_reset_bit(size_limits::max()) );
    BOOST_CHECK_EQUAL( size_limits::min(), b.next_reset_bit(size_limits::max()) );
    BOOST_CHECK_EQUAL( size_limits::min(), e.next_reset_bit(size_limits::max()) );

    BOOST_CHECK_EQUAL( size_limits::max(), z.previous_set_bit(size_limits::min()) );
    BOOST_CHECK_EQUAL( size_limits::max(), a.previous_set_bit(size_limits::min()) );
    BOOST_CHECK_EQUAL( size_limits::max(), b.previous_set_bit(size_limits::min()) );
    BOOST_CHECK_EQUAL( size_limits::max(), e.previous_set_bit(size_limits::min()) );

    BOOST_CHECK_EQUAL( size_limits::max(), z.previous_reset_bit(size_limits::min()) );
    BOOST_CHECK_EQUAL( size_limits::max(), a.previous_reset_bit(size_limits::min()) );
    BOOST_CHECK_EQUAL( size_limits::max(), b.previous_reset_bit(size_limits::min()) );
    BOOST_CHECK_EQUAL( size_limits::max(), e.previous_reset_bit(size_limits::min()) );

    // checking zero
    BOOST_CHECK_EQUAL( 1u, z.next_reset_bit(0) );
    BOOST_CHECK_EQUAL( size_limits::min(), z.next_set_bit(0) );
    BOOST_CHECK_EQUAL( 23u, z.previous_reset_bit(24) );
    BOOST_CHECK_EQUAL( size_limits::max(), z.previous_set_bit(24) );

    // more checks
    BOOST_CHECK_EQUAL( 1u, a.next_reset_bit(0) );
    BOOST_CHECK_EQUAL( size_limits::min(), a.next_set_bit(0) );
    BOOST_CHECK_EQUAL( 23u, a.previous_reset_bit(24) );
    BOOST_CHECK_EQUAL( 0u, a.previous_set_bit(24) );

    BOOST_CHECK_EQUAL( 1u, b.next_reset_bit(0) );
    BOOST_CHECK_EQUAL( 2u, b.next_set_bit(0) );
    BOOST_CHECK_EQUAL( 3u, b.next_reset_bit(1) );
    BOOST_CHECK_EQUAL( 2u, b.next_set_bit(1) );
    BOOST_CHECK_EQUAL( 3u, b.next_reset_bit(2) );
    BOOST_CHECK_EQUAL( size_limits::min(), b.next_set_bit(2) );
    BOOST_CHECK_EQUAL( 4u, b.next_reset_bit(3) );
    BOOST_CHECK_EQUAL( size_limits::min(), b.next_set_bit(3) );

    BOOST_CHECK_EQUAL( 3u, b.previous_reset_bit(4) );
    BOOST_CHECK_EQUAL( 2u, b.previous_set_bit(4) );
    BOOST_CHECK_EQUAL( 1u, b.previous_reset_bit(3) );
    BOOST_CHECK_EQUAL( 2u, b.previous_set_bit(3) );
    BOOST_CHECK_EQUAL( 1u, b.previous_reset_bit(2) );
    BOOST_CHECK_EQUAL( 0u, b.previous_set_bit(2) );
    BOOST_CHECK_EQUAL( size_limits::max(), b.previous_reset_bit(1) );
    BOOST_CHECK_EQUAL( 0u, b.previous_set_bit(1) );

    BOOST_CHECK_EQUAL( 1u, e.next_reset_bit(0) );
    BOOST_CHECK_EQUAL( wd + 1u, e.next_set_bit(0) );

    BOOST_CHECK_EQUAL( size_limits::max(), e.previous_set_bit(wd + 1) );
    BOOST_CHECK_EQUAL( wd, e.previous_reset_bit(wd + 1) );
    BOOST_CHECK_EQUAL( wd + 2, e.next_reset_bit(wd + 1) );
    BOOST_CHECK_EQUAL( 2 * wd + 3, e.next_set_bit(wd + 1) );

    BOOST_CHECK_EQUAL( wd + 1, e.previous_set_bit(3 * wd / 2 + 4) );
    BOOST_CHECK_EQUAL( 3 * wd / 2 + 3, e.previous_reset_bit(3 * wd / 2 + 4) );
    BOOST_CHECK_EQUAL( 3 * wd / 2 + 5, e.next_reset_bit(3 * wd / 2 + 4) );
    BOOST_CHECK_EQUAL( 2 * wd + 3, e.next_set_bit(3 * wd / 2 + 4) );

    BOOST_CHECK_EQUAL( wd + 1, e.previous_set_bit(2 * wd + 3) );
    BOOST_CHECK_EQUAL( 2 * wd + 2, e.previous_reset_bit(2 * wd + 3) );
    BOOST_CHECK_EQUAL( 2 * wd + 4, e.next_reset_bit(2 * wd + 3) );
    BOOST_CHECK_EQUAL( size_limits::min(), e.next_set_bit(2 * wd + 3) );

    BOOST_CHECK_EQUAL( 2 * wd + 3, e.previous_set_bit(size_limits::max()) );
    BOOST_CHECK_EQUAL( size_limits::max() - 1, e.previous_reset_bit(size_limits::max()) );
}

// Unit test for the scale factor
void
bigwhole_scale_unit_test
(
)
{
    using boost::math::big_whole;
    using std::size_t;

    typedef std::valarray<size_t>  va_size_t;

    size_t const     wd = wlimits_type::digits;
    size_t const     ei[] = { wd + 1, 2 * wd + 3 };
    size_t const     es = sizeof( ei ) / sizeof( ei[0] );
    va_size_t const  ev( ei, es );
    size_t const     fi[] = { 1, wd - 1, wd, wd + 2 };
    size_t const     fs = sizeof( fi ) / sizeof( fi[0] );
    va_size_t const  fv( fi, fs );

    big_whole const  z;
    big_whole const  a( 1 );
    big_whole const  b( 5 );
    big_whole const  c( 10 );
    big_whole const  d( 100 );
    big_whole const  e( ev );
    big_whole const  f( fv );

    BOOST_CHECK_EQUAL( 0u, z.scale() );
    BOOST_CHECK_EQUAL( 0u, a.scale() );
    BOOST_CHECK_EQUAL( 0u, b.scale() );
    BOOST_CHECK_EQUAL( 1u, c.scale() );
    BOOST_CHECK_EQUAL( 2u, d.scale() );
    BOOST_CHECK_EQUAL( wd + 1u, e.scale() );
    BOOST_CHECK_EQUAL( 1u, (e >> wd).scale() );
    BOOST_CHECK_EQUAL( wd + 5u, (e << 4).scale() );
    BOOST_CHECK_EQUAL( 1u, f.scale() );
}


// Unit test program
boost::unit_test_framework::test_suite *
init_unit_test_suite
(
    int         ,   // "argc" is unused
    char *      []  // "argv" is unused
)
{
    boost::unit_test_framework::test_suite *  test
     = BOOST_TEST_SUITE( "big_whole test" );

    test->add( BOOST_TEST_CASE(basic_bigwhole_unit_test) );
    test->add( BOOST_TEST_CASE(bigwhole_multi_bit_check_unit_test) );
    test->add( BOOST_TEST_CASE(bigwhole_reverse_unit_test) );

    test->add( BOOST_TEST_CASE(bigwhole_all_bit_reset_unit_test) );
    test->add( BOOST_TEST_CASE(bigwhole_single_bit_reset_unit_test) );
    test->add( BOOST_TEST_CASE(bigwhole_group_bit_reset_unit_test) );
    test->add( BOOST_TEST_CASE(bigwhole_single_bit_set_unit_test) );
    test->add( BOOST_TEST_CASE(bigwhole_group_bit_set_unit_test) );
    test->add( BOOST_TEST_CASE(bigwhole_single_bit_flip_unit_test) );
    test->add( BOOST_TEST_CASE(bigwhole_group_bit_flip_unit_test) );
    test->add( BOOST_TEST_CASE(bigwhole_group_bit_assign_unit_test) );

    test->add( BOOST_TEST_CASE(bigwhole_is_even_unit_test) );
    test->add( BOOST_TEST_CASE(bigwhole_compare_unit_test) );

    test->add( BOOST_TEST_CASE(bigwhole_bitwise_and_unit_test) );
    test->add( BOOST_TEST_CASE(bigwhole_bitwise_or_unit_test) );
    test->add( BOOST_TEST_CASE(bigwhole_bitwise_xor_unit_test) );

    test->add( BOOST_TEST_CASE(bigwhole_left_shift_unit_test) );
    test->add( BOOST_TEST_CASE(bigwhole_right_shift_unit_test) );

    test->add( BOOST_TEST_CASE(bigwhole_bitwise_and_not_unit_test) );

    test->add( BOOST_TEST_CASE(bigwhole_unary_plus_minus_unit_test) );
    test->add( BOOST_TEST_CASE(bigwhole_double_plus_minus_unit_test) );
    test->add( BOOST_TEST_CASE(bigwhole_abs_sgn_unit_test) );
    test->add( BOOST_TEST_CASE(bigwhole_binary_plus_minus_unit_test) );

    test->add( BOOST_TEST_CASE(bigwhole_intersects_unit_test) );
    test->add( BOOST_TEST_CASE(bigwhole_bit_search_unit_test) );
    test->add( BOOST_TEST_CASE(bigwhole_scale_unit_test) );

    return test;
}
