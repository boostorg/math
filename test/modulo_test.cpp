//  Boost modular arithmetic test program file  ------------------------------//

//  Copyright 2002 Daryle Walker.  Use, modification, and distribution are
//  subject to the Boost Software License, Version 1.0.  (See accompanying file
//  LICENSE_1_0.txt or a copy at <http://www.boost.org/LICENSE_1_0.txt>.)

//  See <http://www.boost.org/libs/math/> for the library's home page.

//  Revision History
//  20 Sep 2002  Initial version (Daryle Walker)


#include <boost/config.hpp>        // for BOOST_NO_MEMBER_TEMPLATES
#include <boost/cstdlib.hpp>       // for boost::exit_success
#include <boost/math/modulo.hpp>   // for boost::math::modulo, etc.
#include <boost/test/minimal.hpp>  // for main, BOOST_TEST

#include <iomanip>   // for std::setw
#include <ios>       // for std::ios_base
#include <iostream>  // for std::cout
#include <ostream>   // for std::ostream, std::endl
#include <sstream>   // for std::ostringstream, etc.
#include <string>    // for std::string


// Control macro
#ifndef CONTROL_BACKWARD_CONVERSION_TEST
#define CONTROL_BACKWARD_CONVERSION_TEST  0
#endif

// Macros to compact code
#define PRIVATE_PRINT_UNARY_TABLE( Modulus, Op, Stream ) do { \
 typedef ::boost::math::modulo< (Modulus) >  modulo_type;     \
 (Stream) << "Unary Operation mod-" << (Modulus)              \
  << " (" << #Op << " Value):\n";                             \
 for ( unsigned i = 0u ; i < (Modulus) ; ++i )                \
     (Stream) << ' ' << ::std::setw(6) << i;                  \
 (Stream) << '\n';                                            \
 for ( unsigned j = 0u ; j < (Modulus) ; ++j )                \
     (Stream) << ' ' << ::std::setw(6) << ( Op modulo_type(j) ); \
 (Stream) << '\n' << ::std::endl; } while (false)

#define PRIVATE_PRINT_ALL_UNARY_TABLES( Modulus, Stream ) do { \
 PRIVATE_PRINT_UNARY_TABLE( (Modulus), +, (Stream) );          \
 PRIVATE_PRINT_UNARY_TABLE( (Modulus), -, (Stream) );          \
 PRIVATE_PRINT_UNARY_TABLE( (Modulus), ~, (Stream) );          \
 PRIVATE_PRINT_UNARY_TABLE( (Modulus), (bool), (Stream) );     \
 PRIVATE_PRINT_UNARY_TABLE( (Modulus), !, (Stream) );          \
 PRIVATE_PRINT_UNARY_TABLE( (Modulus), ++, (Stream) );         \
 PRIVATE_PRINT_UNARY_TABLE( (Modulus), --, (Stream) ); } while (false)

#define PRIVATE_PRINT_BINARY_TABLE( Modulus, Op, Stream ) do {  \
 typedef ::boost::math::modulo< (Modulus) >  modulo_type;       \
 (Stream) << "Binary Operation mod-" << (Modulus)               \
  << " (Left " << #Op << " Top):\n";                            \
 (Stream) << ::std::setw(6) << ' ';                             \
 for ( unsigned i = 0u ; i < (Modulus) ; ++i )                  \
     (Stream) << ' ' << ::std::setw(6) << i;                    \
 (Stream) << '\n';                                              \
 for ( unsigned j = 0u ; j < (Modulus) ; ++j ) {                \
     (Stream) << ::std::setw(6) << j;                           \
     for ( unsigned k = 0u ; k < (Modulus) ; ++k ) {            \
         ::std::ostringstream  oss;                             \
         ::std::string         answer;                          \
         oss.copyfmt( (Stream) );                               \
         try { oss << ( modulo_type(j) Op modulo_type(k) );     \
          answer = oss.str(); } catch (...) { answer.assign( "-" ); } \
         (Stream) << ' ' << ::std::setw(6) << answer;           \
     } (Stream) << '\n';                                        \
 } (Stream) << ::std::endl; } while (false)

#define PRIVATE_PRINT_ALL_BINARY_TABLES( Modulus, Stream ) do { \
 PRIVATE_PRINT_BINARY_TABLE( (Modulus), +, (Stream) );          \
 PRIVATE_PRINT_BINARY_TABLE( (Modulus), -, (Stream) );          \
 PRIVATE_PRINT_BINARY_TABLE( (Modulus), *, (Stream) );          \
 PRIVATE_PRINT_BINARY_TABLE( (Modulus), /, (Stream) );          \
 PRIVATE_PRINT_BINARY_TABLE( (Modulus), &, (Stream) );          \
 PRIVATE_PRINT_BINARY_TABLE( (Modulus), |, (Stream) );          \
 PRIVATE_PRINT_BINARY_TABLE( (Modulus), ^, (Stream) );          \
 PRIVATE_PRINT_BINARY_TABLE( (Modulus), &&, (Stream) );         \
 PRIVATE_PRINT_BINARY_TABLE( (Modulus), ||, (Stream) );         \
 PRIVATE_PRINT_BINARY_TABLE( (Modulus), ==, (Stream) );         \
 PRIVATE_PRINT_BINARY_TABLE( (Modulus), !=, (Stream) );         \
 PRIVATE_PRINT_BINARY_TABLE( (Modulus), >=, (Stream) );         \
 PRIVATE_PRINT_BINARY_TABLE( (Modulus), <=, (Stream) );         \
 PRIVATE_PRINT_BINARY_TABLE( (Modulus), >, (Stream) );          \
 PRIVATE_PRINT_BINARY_TABLE( (Modulus), <, (Stream) ); } while (false)

#define PRIVATE_PRINT_BINARY_FUNC_TABLE( Modulus, Func, Stream ) do { \
 typedef ::boost::math::modulo< (Modulus) >  modulo_type; \
 (Stream) << "Binary Operation mod-" << (Modulus)       \
  << " " << #Func << "(Left, Top):\n";                  \
 (Stream) << ::std::setw(6) << ' ';                     \
 for ( unsigned i = 0u ; i < (Modulus) ; ++i )          \
     (Stream) << ' ' << ::std::setw(6) << i;            \
 (Stream) << '\n';                                      \
 for ( unsigned j = 0u ; j < (Modulus) ; ++j ) {        \
     (Stream) << ::std::setw(6) << j;                   \
     for ( unsigned k = 0u ; k < (Modulus) ; ++k ) {    \
         ::std::ostringstream  oss;                     \
         ::std::string         answer;                  \
         oss.copyfmt( (Stream) );                       \
         try { oss << (Func)( modulo_type(j), k );      \
          answer = oss.str(); } catch (...) { answer.assign( "-" ); } \
         (Stream) << ' ' << ::std::setw(6) << answer;   \
     } (Stream) << '\n';                                \
 } (Stream) << ::std::endl; } while (false)

#define PRIVATE_PRINT_ALL_TABLES( Modulus, Stream ) do { \
 PRIVATE_PRINT_ALL_UNARY_TABLES( (Modulus), (Stream) );  \
 PRIVATE_PRINT_ALL_BINARY_TABLES( (Modulus), (Stream) ); } while (false)


// Main testing function
int
test_main
(
    int         ,   // "argc" is unused
    char *      []  // "argv" is unused
)
{
    using boost::math::modulo;
    using std::cout;
    using boost::math::pow;
    using std::endl;

    typedef modulo<6>  modulo6_t;
    typedef modulo<3>  modulo3_t;

    // Try out most of the operations
    cout.setf( std::ios_base::boolalpha );
    PRIVATE_PRINT_ALL_TABLES( 7, cout );
    PRIVATE_PRINT_BINARY_FUNC_TABLE( 7, pow, cout );

    // Test input
    std::istringstream  as( "[3]" ), bs( "[2]" ), cs( "-4" );
    modulo6_t           a, b, c;

    cout << "Doing input tests." << endl;
    BOOST_TEST( as >> a );
    BOOST_TEST( residue(a) == 3 );
    BOOST_TEST( bs >> b );
    BOOST_TEST( residue(b) == 2 );
    BOOST_TEST( !(cs >> c) );
    BOOST_TEST( (cs.clear(), cs.str( "[ty" ), !( cs >> c )) );
    BOOST_TEST( (cs.clear(), cs.str( "[5?" ), !( cs >> c )) );
    BOOST_TEST( (cs.clear(), cs.str( "[1]r0]]" ), ( cs >> c )) );
    BOOST_TEST( residue(c) == 1 );

    // Test post-(in/de)crement
    cout << "Doing post-(in/de)crement tests." << endl;
    BOOST_TEST( (c++, residue( c ) == 2) );
    BOOST_TEST( (c--, residue( c ) == 1) );

    // Test swap
    cout << "Doing swap test." << endl;
    swap( a, b );
    BOOST_TEST( (residue( a ) == 2) && (residue( b ) == 3) );

    #ifndef BOOST_NO_MEMBER_TEMPLATES
    // Test cross-conversions
    cout << "Doing cross-conversion tests." << endl;
    BOOST_TEST( residue(modulo3_t( modulo6_t(0) )) == 0 );
    BOOST_TEST( residue(modulo3_t( modulo6_t(1) )) == 1 );
    BOOST_TEST( residue(modulo3_t( modulo6_t(2) )) == 2 );
    BOOST_TEST( residue(modulo3_t( modulo6_t(3) )) == 0 );
    BOOST_TEST( residue(modulo3_t( modulo6_t(4) )) == 1 );
    BOOST_TEST( residue(modulo3_t( modulo6_t(5) )) == 2 );

    #if CONTROL_BACKWARD_CONVERSION_TEST
    // Should cause an error #if compiled
    BOOST_TEST( residue(modulo6_t( modulo3_t(5) )) == 2 );
    #endif
    #endif

    // Test quasi-traditional notation
    cout << "Doing cross-comparison tests." << endl;
    BOOST_TEST( 17 == modulo<5>(2) );
    BOOST_TEST( 5 != modulo<2>(0) );
    BOOST_TEST( modulo<5>(3) != 17 );
    BOOST_TEST( modulo<2>(1) == 5 );

    // Test assignment
    cout << "Doing assignment tests." << endl;
    BOOST_TEST( (c.assign( 4 ), residue( c ) == 4) );
    BOOST_TEST( (c.assign( 17 ), residue( c ) == 5) );

    // Test inversion status
    cout << "Doing invertibilty tests." << endl;
    BOOST_TEST( !modulo6_t(0).is_invertible() );
    BOOST_TEST( modulo6_t(1).is_invertible() );
    BOOST_TEST( !modulo6_t(2).is_invertible() );
    BOOST_TEST( !modulo6_t(3).is_invertible() );
    BOOST_TEST( !modulo6_t(4).is_invertible() );
    BOOST_TEST( modulo6_t(5).is_invertible() );
    BOOST_TEST( !modulo3_t(0).is_invertible() );
    BOOST_TEST( modulo3_t(1).is_invertible() );
    BOOST_TEST( modulo3_t(2).is_invertible() );

    // Test Chinese remainder theorem
    cout << "Doing Chinese Remainder Theorem tests." << endl;
    BOOST_TEST( boost::math::chinese_remainder(modulo<7>( 3 ), modulo<13>( 7 ))
     == modulo<91>(59) );

    return boost::exit_success;
}
