//    boost atanh.hpp header file

//  (C) Copyright Hubert Holin 2001.
//  Distributed under the Boost Software License, Version 1.0. (See
//  accompanying file LICENSE_1_0.txt or copy at
//  http://www.boost.org/LICENSE_1_0.txt)

// See http://www.boost.org for updates, documentation, and revision history.

#ifndef BOOST_ATANH_HPP
#define BOOST_ATANH_HPP


#include <cmath>
#include <boost/config.hpp>
#include <boost/math/tools/precision.hpp>
#include <boost/math/tools/error_handling.hpp>

// This is the inverse of the hyperbolic tangent function.

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
        
        // This is the main fare
        
        template<typename T>
        inline T    atanh(const T x)
        {
            using    ::std::abs;
            using    ::std::sqrt;
            using    ::std::log;
            
            using    ::std::numeric_limits;
            
            T const            one = static_cast<T>(1);
            T const            two = static_cast<T>(2);
            
            static T const    taylor_2_bound = sqrt(tools::epsilon<T>());
            static T const    taylor_n_bound = sqrt(taylor_2_bound);
            
            if        (x < -one)
            {
               return tools::domain_error<T>(
                  BOOST_CURRENT_FUNCTION,
                  "atanh requires x >= -1, but got x = %1%.", x);
            }
            else if    (x < -one + tools::epsilon<T>())
            {
               // -Infinity:
               return -tools::overflow_error<T>(
                  BOOST_CURRENT_FUNCTION);
            }
            else if    (x > one - tools::epsilon<T>())
            {
               // Infinity:
               return -tools::overflow_error<T>(
                  BOOST_CURRENT_FUNCTION);
            }
            else if    (x > +one)
            {
               return tools::domain_error<T>(
                  BOOST_CURRENT_FUNCTION,
                  "atanh requires x <= 1, but got x = %1%.", x);
            }
            else if    (abs(x) >= taylor_n_bound)
            {
                return(log( (one + x) / (one - x) ) / two);
            }
            else
            {
                // approximation by taylor series in x at 0 up to order 2
                T    result = x;
                
                if    (abs(x) >= taylor_2_bound)
                {
                    T    x3 = x*x*x;
                    
                    // approximation by taylor series in x at 0 up to order 4
                    result += x3/static_cast<T>(3);
                }
                
                return(result);
            }
        }

    }
}

#endif /* BOOST_ATANH_HPP */

