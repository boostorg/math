//  boost atanh.hpp header file

//  (C) Copyright Hubert Holin 2001. Permission to copy, use, modify, sell and
//  distribute this software is granted provided this copyright notice appears
//  in all copies. This software is provided "as is" without express or implied
//  warranty, and with no claim as to its suitability for any purpose.

// See http://www.boost.org for updates, documentation, and revision history.

#ifndef BOOST_ATANH_HPP
#define BOOST_ATANH_HPP


#include <cmath>
#include <limits>
#include <string>
#include <stdexcept>


// This is the inverse of the hyperbolic tangent function.

namespace boost
{
    // These are implementation details (for main fare see below)
    
    namespace detail
    {
        template    <
                        typename T,
                        bool InfinitySupported
                    >
        struct atanh_helper1_t
        {
            static T    get_pos_infinity()
            {
                return(+std::numeric_limits<T>::infinity());
            }
            
            static T    get_neg_infinity()
            {
                return(-std::numeric_limits<T>::infinity());
            }
        };  // boost::detail::atanh_helper1_t
        
        
        template<typename T>
        struct atanh_helper1_t<T, false>
        {
            static T    get_pos_infinity()
            {
                ::std::string        error_reporting("Argument to atanh is +1 (result: +Infinity)!");
                ::std::out_of_range  bad_argument(error_reporting);
                
                throw(bad_argument);
            }
            
            static T    get_neg_infinity()
            {
                ::std::string        error_reporting("Argument to atanh is -1 (result: -Infinity)!");
                ::std::out_of_range  bad_argument(error_reporting);
                
                throw(bad_argument);
            }
        };  // boost::detail::atanh_helper1_t
        
        
        template    <
                        typename T,
                        bool QuietNanSupported
                    >
        struct atanh_helper2_t
        {
            static T    get_pos_NaN()
            {
                return(+std::numeric_limits<T>::quiet_NaN());
            }
            
            static T    get_neg_NaN()
            {
                return(-std::numeric_limits<T>::quiet_NaN());
            }
        };  // boost::detail::atanh_helper2_t
        
        
        template<typename T>
        struct atanh_helper2_t<T, false>
        {
            static T    get_pos_NaN()
            {
                ::std::string        error_reporting("Argument to atanh is strictly greater than +1!");
                ::std::domain_error  bad_argument(error_reporting);
                
                throw(bad_argument);
            }
            
            static T    get_neg_NaN()
            {
                ::std::string        error_reporting("Argument to atanh is strictly smaller than -1!");
                ::std::domain_error  bad_argument(error_reporting);
                
                throw(bad_argument);
            }
        };  // boost::detail::atanh_helper2_t
        
    }  // boost::detail
    
    
    // This is the main fare
    
    template<typename T>
    inline T    atanh(const T x)
    {
        using    ::std::abs;
        using    ::std::sqrt;
        
        using    ::std::numeric_limits;
        
        typedef    detail::atanh_helper1_t<T, std::numeric_limits<T>::has_infinity>    helper1_type;
        typedef    detail::atanh_helper2_t<T, std::numeric_limits<T>::has_quiet_NaN>   helper2_type;
        
        static T const    e0 = sqrt(numeric_limits<T>::epsilon());
        
        T const    one = static_cast<T>(1);
        T const    two = static_cast<T>(2);
        
        if        (x < -one)
        {
            return(helper2_type::get_neg_NaN());
        }
        else if    (x == -one)
        {
            return(helper1_type::get_neg_infinity());
        }
        else if    (x == +one)
        {
            return(helper1_type::get_pos_infinity());
        }
        else if    (x > +one)
        {
            return(helper2_type::get_pos_NaN());
        }
        else if    (abs(x) < e0)
        {
            return(x);
        }
        else        // -one < x <= -e0 or +e0 <= x < +one
        {
            using     ::std::log;
            
            return(log( (one + x) / (one - x) ) / two);
        }
    }
}

#endif /* BOOST_ATANH_HPP */
