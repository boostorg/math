//  Boost math/big_whole.hpp header file  ------------------------------------//

//  Copyright 2004 Daryle Walker.  Use, modification, and distribution are
//  subject to the Boost Software License, Version 1.0.  (See accompanying file
//  LICENSE_1_0.txt or a copy at <http://www.boost.org/LICENSE_1_0.txt>.)

//  See <http://www.boost.org/libs/math/> for the library's home page.

#ifndef BOOST_MATH_BIG_WHOLE_HPP
#define BOOST_MATH_BIG_WHOLE_HPP

#include <boost/math_fwd.hpp>             // self include
#include <boost/math/big_whole_core.hpp>  // for boost::math::big_whole

#include <istream>  // for std::basic_istream
#include <ostream>  // for std::basic_ostream
#include <string>   // for std::basic_string


namespace boost
{
namespace math
{


//  Class/template/function/operator forward declarations  -------------------//

template < typename Ch, class Tr >
    std::basic_ostream<Ch, Tr> &  operator <<( std::basic_ostream<Ch, Tr> &os,
     big_whole const &w );

template < typename Ch, class Tr >
    std::basic_istream<Ch, Tr> &  operator >>( std::basic_istream<Ch, Tr> &is,
     big_whole &w );


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

    using std::size_t;

    std::valarray<size_t> const  sz = w.to_bit_indices();

    if ( size_t const  s = sz.size() )
    {
        size_t const               w_len = w.length();
        std::basic_string<Ch, Tr>  str( w_len, os.widen('0') );
        Ch const                   one_char = os.widen( '1' );

        for ( size_t  i = 0 ; i < s ; ++i )
        {
            str.at( w_len - sz[i] - 1 ) = one_char;
        }

        return os.write( str.data(), str.length() );
    }
    else
    {
        return os.put( os.widen('0') );
    }
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
