// test file for special_functions.hpp

//  (C) Copyright Hubert Holin 2001. Permission to copy, use, modify, sell and
//  distribute this software is granted provided this copyright notice appears
//  in all copies. This software is provided "as is" without express or implied
//  warranty, and with no claim as to its suitability for any purpose.


#include <iostream>
#include <iomanip>
#include <complex>

#include <boost/math/special_functions/atanh.hpp>
#include <boost/math/special_functions/sinc.hpp>
#include <boost/math/special_functions/sinhc.hpp>


#define BOOST_INCLUDE_MAIN  // for testing, include rather than link
#include <boost/test/test_tools.hpp>


int    test_main(int, char *[])

{
    using ::std::abs;
    using ::std::tanh;
    using ::std::log;
    using ::std::sinh;
    using ::std::sin;
        
    using ::std::numeric_limits;
    
    using ::boost::math::atanh;
    using ::boost::math::sinc_pi;
    using ::boost::math::sinhc_pi;
    
    
    // tests for evaluation by humans
     
    for    (int i = 0; i <= 100; i++)
    {
        float        xf =
            tanh(static_cast<float>(i-50)/static_cast<float>(5));
        
        double       xd =
            tanh(static_cast<double>(i-50)/static_cast<double>(5));
        
        long double  xl =
            tanh(static_cast<long double>(i-50)/static_cast<long double>(5));
        
        if    (
                std::numeric_limits<float>::has_infinity &&
                std::numeric_limits<double>::has_infinity &&
                std::numeric_limits<long double>::has_infinity
            )
        {
            ::std::cout << ::std::setw(15)
                        << atanh<float>(xf)
                        << ::std::setw(15)
                        << atanh<double>(xd)
                        << ::std::setw(15)
                        << atanh<long double>(xl)
                        << ::std::endl;
        }
        else
        {
            if    (
                    (abs(xf-static_cast<float>(1)) <
                        numeric_limits<float>::epsilon())||
                    (abs(xf+static_cast<float>(1)) <
                        numeric_limits<float>::epsilon())||
                    (abs(xf-static_cast<double>(1)) <
                        numeric_limits<double>::epsilon())||
                    (abs(xf+static_cast<double>(1)) <
                        numeric_limits<double>::epsilon())||
                    (abs(xf-static_cast<long double>(1)) <
                        numeric_limits<long double>::epsilon())||
                    (abs(xf+static_cast<long double>(1)) <
                        numeric_limits<long double>::epsilon())
                )
            {
                ::std::cout << "Platform's numerics may lack precision."
                            << ::std::endl;
            }
            else
            {
                ::std::cout << ::std::setw(15)
                            << atanh<float>(xf)
                            << ::std::setw(15)
                            << atanh<double>(xd)
                            << ::std::setw(15)
                            << atanh<long double>(xl)
                            << ::std::endl;
            }
        }
    }
    
    ::std::cout << ::std::endl;
    
    for    (int i = 0; i <= 100; i++)
    {
        ::std::cout << ::std::setw(15)
                    << sinc_pi<float>(static_cast<float>(i-50)/
                                                static_cast<float>(50))
                    << ::std::setw(15)
                    << sinc_pi<double>(static_cast<double>(i-50)/
                                                static_cast<double>(50))
                    << ::std::setw(15)
                    << sinc_pi<long double>(static_cast<long double>(i-50)/
                                                static_cast<long double>(50))
                    << ::std::endl;
    }
    
    ::std::cout << ::std::endl;
    
    for    (int i = 0; i <= 100; i++)
    {
        ::std::cout << ::std::setw(15)
                    << sinhc_pi<float>(static_cast<float>(i-50)/
                                                static_cast<float>(50))
                    << ::std::setw(15)
                    << sinhc_pi<double>(static_cast<double>(i-50)/
                                                static_cast<double>(50))
                    << ::std::setw(15)
                    << sinhc_pi<long double>(static_cast<long double>(i-50)/
                                                static_cast<long double>(50))
                    << ::std::endl;
    }
    
    ::std::cout << ::std::endl;
    
    
    // tests for evaluation by scripts
    
#define    BOOST_SPECIAL_FUNCTIONS_TEST_ATANH(type)         \
                                                            \
    BOOST_TEST(abs(atanh<type>(static_cast<type>(0))) <=    \
        numeric_limits<type>::epsilon());                   \
                                                            \
    BOOST_TEST(abs(atanh<type>(static_cast<type>(3)/5)-     \
        log(static_cast<type>(2))) <=                       \
        numeric_limits<type>::epsilon());                   \
                                                            \
    BOOST_TEST(abs(atanh<type>(static_cast<type>(-3)/5)+    \
        log(static_cast<type>(2))) <=                       \
        numeric_limits<type>::epsilon());
    
    
#define    BOOST_SPECIAL_FUNCTIONS_TEST_SINC_PI(type)           \
                                                                \
    BOOST_TEST(abs(sinc_pi<type>(static_cast<type>(0))-         \
        static_cast<type>(1)) <=                                \
        numeric_limits<type>::epsilon());                       \
                                                                \
    BOOST_TEST(abs(sinc_pi<type>(::std::complex<type>(0, 1))-   \
        ::std::complex<type>(sinh(static_cast<type>(1)))) <=    \
        numeric_limits<type>::epsilon());
    
    
#define    BOOST_SPECIAL_FUNCTIONS_TEST_SINHC_PI(type)          \
                                                                \
    BOOST_TEST(abs(sinhc_pi<type>(static_cast<type>(0))-        \
        static_cast<type>(1)) <=                                \
        numeric_limits<type>::epsilon());                       \
                                                                \
    BOOST_TEST(abs(sinhc_pi<type>(::std::complex<type>(0, 1))-  \
        ::std::complex<type>(sin(static_cast<type>(1)))) <=     \
        numeric_limits<type>::epsilon());
    
    
#define    BOOST_SPECIAL_FUNCTIONS_TEST(type)   \
    BOOST_SPECIAL_FUNCTIONS_TEST_ATANH(type)    \
    BOOST_SPECIAL_FUNCTIONS_TEST_SINC_PI(type)  \
    BOOST_SPECIAL_FUNCTIONS_TEST_SINHC_PI(type)
    
    
    BOOST_SPECIAL_FUNCTIONS_TEST(float)
    BOOST_SPECIAL_FUNCTIONS_TEST(double)
    BOOST_SPECIAL_FUNCTIONS_TEST(long double)
    
    
#undef    BOOST_SPECIAL_FUNCTIONS_TEST

#undef    BOOST_SPECIAL_FUNCTIONS_TEST_ATANH
#undef    BOOST_SPECIAL_FUNCTIONS_TEST_SINC_PI
#undef    BOOST_SPECIAL_FUNCTIONS_TEST_SINHC_PI
    
    return(::boost::exit_success);
}
