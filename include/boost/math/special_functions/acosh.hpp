//    boost asinh.hpp header file

//  (C) Copyright Eric Ford 2001 & Hubert Holin.
//  Distributed under the Boost Software License, Version 1.0. (See
//  accompanying file LICENSE_1_0.txt or copy at
//  http://www.boost.org/LICENSE_1_0.txt)

// See http://www.boost.org for updates, documentation, and revision history.

#ifndef BOOST_ACOSH_HPP
#define BOOST_ACOSH_HPP


#include <cmath>
#include <boost/config.hpp>
#include <boost/math/tools/precision.hpp>
#include <boost/math/tools/error_handling.hpp>

// This is the inverse of the hyperbolic cosine function.

namespace boost
{
    namespace math
    {
#if defined(__GNUC__) && (__GNUC__ < 3)
        // gcc 2.x ignores function scope using declarations,
        // put them in the scope of the enclosing namespace instead:
        
        using    ::std::abs;
        using    ::std::sqrt;
        using    ::std::log;
        
        using    ::std::numeric_limits;
#endif
        
        template<typename T>
        inline T    acosh(const T x)
        {
            using    ::std::abs;
            using    ::std::sqrt;
            using    ::std::log;
            
            T const    one = static_cast<T>(1);
            T const    two = static_cast<T>(2);
            
            static T const    taylor_2_bound = sqrt(tools::epsilon<T>());
            static T const    taylor_n_bound = sqrt(taylor_2_bound);
            static T const    upper_taylor_2_bound = one/taylor_2_bound;
            
            if(x < one)
            {
               return tools::domain_error<T>(
                  BOOST_CURRENT_FUNCTION,
                  "acosh requires x >= 1, but got x = %1%.", x);
            }
            else if    (x >= taylor_n_bound)
            {
                if    (x > upper_taylor_2_bound)
                {
                    // approximation by laurent series in 1/x at 0+ order from -1 to 0
                    return( log( x*two) );
                }
                else
                {
                    return( log( x + sqrt(x*x-one) ) );
                }
            }
            else
            {
                T    y = sqrt(x-one);
                
                // approximation by taylor series in y at 0 up to order 2
                T    result = y;
                
                if    (y >= taylor_2_bound)
                {
                    T    y3 = y*y*y;
                    
                    // approximation by taylor series in y at 0 up to order 4
                    result -= y3/static_cast<T>(12);
                }
                
                return(sqrt(static_cast<T>(2))*result);
            }
        }
    }
}

#endif /* BOOST_ACOSH_HPP */


