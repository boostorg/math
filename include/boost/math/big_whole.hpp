//  Boost math/big_whole.hpp header file  ------------------------------------//

//  Copyright 2004 Daryle Walker.  Use, modification, and distribution are
//  subject to the Boost Software License, Version 1.0.  (See accompanying file
//  LICENSE_1_0.txt or a copy at <http://www.boost.org/LICENSE_1_0.txt>.)

//  See <http://www.boost.org/libs/math/> for the library's home page.

#ifndef BOOST_MATH_BIG_WHOLE_HPP
#define BOOST_MATH_BIG_WHOLE_HPP

#include <boost/math_fwd.hpp>             // self include
#include <boost/math/big_whole_core.hpp>  // for boost::math::big_whole

#include <cstddef>   // for std::size_t
#include <istream>   // for std::basic_istream
#include <ostream>   // for std::basic_ostream
#include <string>    // for std::basic_string
#include <valarray>  // for std::valarray


namespace boost
{
namespace math
{


//  Class/template/function/operator forward declarations  -------------------//

template < typename Ch, class Tr, class Al >
    std::basic_string<Ch, Tr, Al>  bigwhole_to_bitstring( big_whole const &w,
     Ch zero_char, Ch one_char );
template < typename Ch, class Tr, class Al >
    void  copy_bigwhole_to_bitstring( std::basic_string<Ch, Tr, Al> &s,
     big_whole const &w, Ch zero_char, Ch one_char );

template < typename Ch, class Tr >
    std::basic_ostream<Ch, Tr> &  operator <<( std::basic_ostream<Ch, Tr> &os,
     big_whole const &w );

template < typename Ch, class Tr >
    std::basic_istream<Ch, Tr> &  operator >>( std::basic_istream<Ch, Tr> &is,
     big_whole &w );


//  Arbitrary-length whole-number bit-string printer definitions  ------------//

template < typename Ch, class Tr, class Al >
std::basic_string<Ch, Tr, Al>
bigwhole_to_bitstring
(
    big_whole const &  w,
    Ch                 zero_char,
    Ch                 one_char
)
{
    using std::size_t;

    std::valarray<size_t> const    sz = w.to_bit_indices();
    size_t const                   s = sz.size();
    size_t const                   char_count = 1 + ( s ? sz.max() : 0 );
    size_t const                   max_index = char_count - 1;
    std::basic_string<Ch, Tr, Al>  str( char_count, zero_char );

    for ( size_t  i = 0 ; i < s ; ++i )
    {
        str[ max_index - sz[i] ] = one_char;
    }

    return str;
}

template < typename Ch, class Tr, class Al >
inline
void
copy_bigwhole_to_bitstring
(
    std::basic_string<Ch, Tr, Al> &  s,
    big_whole const &                w,
    Ch                               zero_char,
    Ch                               one_char
)
{
    s = bigwhole_to_bitstring<Ch, Tr, Al>( w, zero_char, one_char );
}


//  Arbitrary-length whole-number streaming operator definitions  ------------//

template < typename Ch, class Tr >
inline
std::basic_ostream<Ch, Tr> &
operator <<
(
    std::basic_ostream<Ch, Tr> &  os,
    big_whole const &             w
)
{
    // NOTE: this solution is Temporary!

    std::basic_string<Ch, Tr>  str;

    copy_bigwhole_to_bitstring( str, w, os.widen('0'), os.widen('1') );
    return os.write( str.data(), str.length() );
}

template < typename Ch, class Tr >
inline
std::basic_istream<Ch, Tr> &
operator >>
(
    std::basic_istream<Ch, Tr> &  is,
    big_whole &                   w
)
{
    return is;  // FILL IN LATER!
}


}  // namespace math
}  // namespace boost


#endif  // BOOST_MATH_BIG_WHOLE_HPP
