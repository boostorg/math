//  Boost math/modulo.hpp header file  ---------------------------------------//

//  (C) Copyright Daryle Walker 2002.  Permission to copy, use, modify, sell and
//  distribute this software is granted provided this copyright notice appears
//  in all copies.  This software is provided "as is" without express or implied
//  warranty, and with no claim as to its suitability for any purpose.

//  See http://www.boost.org for updates, documentation, and revision history. 

#ifndef BOOST_MATH_MODULO_HPP
#define BOOST_MATH_MODULO_HPP

#include <boost/math/modulo_core.hpp>  // for boost::math::modulo

#include <ios>         // for std::ios_base, etc.
#include <istream>     // for std::basic_istream
#include <ostream>     // for std::basic_ostream
#include <sstream>     // for std::basic_ostringstream
#include <string>      // for std::basic_string


namespace boost
{
namespace math
{


//  Forward declarations  ----------------------------------------------------//

template < unsigned long Modulus, typename Ch, class Tr >
    std::basic_istream<Ch, Tr> &  operator >>( std::basic_istream<Ch, Tr> &bis,
     modulo<Modulus> &x );

template < unsigned long Modulus, typename Ch, class Tr >
    std::basic_ostream<Ch, Tr> &  operator <<( std::basic_ostream<Ch, Tr> &bos,
     modulo<Modulus> const &x );


//  Modular arithmetic I/O function definitions  -----------------------------//

template < unsigned long Modulus, typename Ch, class Tr >
std::basic_istream<Ch, Tr> &
operator >>
(
    std::basic_istream<Ch, Tr> &  bis,
    modulo<Modulus> &             x
)
{
    using std::ios_base;

    typedef modulo<Modulus>  modulo_type;

    try
    {
        Ch  next;

        if ( (bis >> next) && Tr::eq(next, bis.widen( '[' )) )
        {
            typename modulo_type::value_type  residue;

            if ( bis >> residue )
            {
                if ( (bis >> next) && Tr::eq(next, bis.widen( ']' )) )
                {
                    modulo_type  temp( residue );

                    x.swap( temp );
                }
                else
                {
                    bis.setstate( ios_base::failbit );
                }
            }
            else
            {
                bis.setstate( ios_base::failbit );
            }
        }
        else
        {
            bis.setstate( ios_base::failbit );
        }
    }
    catch ( ios_base::failure & )
    {
        throw;
    }
    catch ( ... )
    {
        bis.setstate( ios_base::failbit );
    }

    return bis;
}

template < unsigned long Modulus, typename Ch, class Tr >
std::basic_ostream<Ch, Tr> &
operator <<
(
    std::basic_ostream<Ch, Tr> &  bos,
    modulo<Modulus> const &       x
)
{
    using std::ios_base;

    // Simulate most of the original stream's configuraion
    std::basic_ostringstream<Ch, Tr>  boss;

    boss.copyfmt( bos );
    boss.setf( ios_base::dec, ios_base::basefield );
    boss.unsetf( ios_base::showpos );
    boss.width( 0 );

    // Print the formatted value to scratch space
    boss << '[' << residue( x ) << ']';

    // Sent the scratch text to the real output
    if ( boss.fail() )
    {
        bos.setstate( ios_base::failbit );
    }
    else
    {
        bos << boss.str();
    }
    return bos;
}


}  // namespace math
}  // namespace boost


#endif  // BOOST_MATH_MODULO_HPP
