//  boost sinhc.hpp header file

//  (C) Copyright Hubert Holin 2001. Permission to copy, use, modify, sell and
//  distribute this software is granted provided this copyright notice appears
//  in all copies. This software is provided "as is" without express or implied
//  warranty, and with no claim as to its suitability for any purpose.

// See http://www.boost.org for updates, documentation, and revision history.

#ifndef BOOST_SINHC_HPP
#define BOOST_SINHC_HPP


#include <cmath>
#include <limits>
#include <string>
#include <stdexcept>


// These are the the "Hyperbolic Sinus Cardinal" functions.

namespace boost
{
    // This is the "Hyperbolic Sinus Cardinal" of index Pi.
    
    template<typename T>
    inline T    sinhc_pi(const T x)
    {
        using    ::std::abs;
        using    ::std::sinh;
        using    ::std::sqrt;
        
        using    ::std::numeric_limits;
        
        static T const    e1 = numeric_limits<T>::epsilon();
        static T const    e2 = sqrt(e1);
        static T const    e3 = sqrt(e2);
        
        if    (abs(x) > e3)
        {
            return(sinh(x)/x);
        }
        else
        {
            T    result = static_cast<T>(1);
            
            if    (abs(x) > e1)
            {
                T    x2 = x*x;
                
                result += x2/static_cast<T>(6);
                
                if    (abs(x) > e2)
                {
                    result += (x2*x2)/static_cast<T>(120);
                }
            }
            
            return(result);
        }
    }
    
    
    template<typename T, template<typename> class U>
    inline U<T>    sinhc_pi(const U<T> x)
    {
        using    ::std::abs;
        using    ::std::sinh;
        using    ::std::sqrt;
        
        using    ::std::numeric_limits;
        
        static T const    e1 = numeric_limits<T>::epsilon();
        static T const    e2 = sqrt(e1);
        static T const    e3 = sqrt(e2);
        
        if    (abs(x) > e3)
        {
            return(sinh(x)/x);
        }
        else
        {
            U<T>    result = static_cast< U<T> >(1);
            
            if    (abs(x) > e1)
            {
                U<T>    x2 = x*x;
                
                result += x2/static_cast<T>(6);
                
                if    (abs(x) > e2)
                {
                    result += (x2*x2)/static_cast<T>(120);
                }
            }
            
            return(result);
        }
    }
}

#endif /* BOOST_SINHC_HPP */
