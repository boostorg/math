// test file for special_functions.hpp

//  (C) Copyright Hubert Holin 2001. Permission to copy, use, modify, sell and
//  distribute this software is granted provided this copyright notice appears
//  in all copies. This software is provided "as is" without express or implied
//  warranty, and with no claim as to its suitability for any purpose.


#include <iostream>
#include <iomanip>

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
        
    using ::std::numeric_limits;
    
    using ::boost::atanh;
    using ::boost::sinc_pi;
    using ::boost::sinhc_pi;
    
    
    // tests for evaluation by humans
     
    for    (int i = 0; i <= 100; i++)
    {
        ::std::cout << ::std::setw(15)
                    << atanh<float>(tanh(static_cast<float>(i-50)/
                                                static_cast<float>(5)))
                    << ::std::setw(15)
                    << atanh<double>(tanh(static_cast<double>(i-50)/
                                                static_cast<double>(5)))
                    << ::std::setw(15)
                    << atanh<long double>(tanh(static_cast<long double>(i-50)/
                                                static_cast<long double>(5)))
                    << ::std::endl;
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
    
    BOOST_TEST(abs(atanh<float>(static_cast<float>(0))) <=
        numeric_limits<float>::epsilon());
    
    BOOST_TEST(abs(atanh<double>(static_cast<double>(0))) <=
        numeric_limits<double>::epsilon());
    
    BOOST_TEST(abs(atanh<long double>(static_cast<long double>(0))) <=
        numeric_limits<long double>::epsilon());
    
    
    BOOST_TEST(abs(atanh<float>(static_cast<float>(3)/5)-
        log(static_cast<float>(2))) <=
        numeric_limits<float>::epsilon());
    
    BOOST_TEST(abs(atanh<double>(static_cast<double>(3)/5)-
        log(static_cast<double>(2))) <=
        numeric_limits<double>::epsilon());
    
    BOOST_TEST(abs(atanh<long double>(static_cast<long double>(3)/5)-
        log(static_cast<long double>(2))) <=
        numeric_limits<long double>::epsilon());
    
    
    BOOST_TEST(abs(atanh<float>(static_cast<float>(-3)/5)+
        log(static_cast<float>(2))) <=
        numeric_limits<float>::epsilon());
    
    BOOST_TEST(abs(atanh<double>(static_cast<double>(-3)/5)+
        log(static_cast<double>(2))) <=
        numeric_limits<double>::epsilon());
    
    BOOST_TEST(abs(atanh<long double>(static_cast<long double>(-3)/5)+
        log(static_cast<long double>(2))) <=
        numeric_limits<long double>::epsilon());
    
    
    BOOST_TEST(abs(sinc_pi<float>(static_cast<float>(0))-
        static_cast<float>(1)) <=
        numeric_limits<float>::epsilon());
    
    BOOST_TEST(abs(sinc_pi<double>(static_cast<double>(0))-
        static_cast<double>(1)) <=
        numeric_limits<double>::epsilon());
    
    BOOST_TEST(abs(sinc_pi<long double>(static_cast<long double>(0))-
        static_cast<long double>(1)) <=
        numeric_limits<long double>::epsilon());
    
    
    BOOST_TEST(abs(sinhc_pi<float>(static_cast<float>(0))-
        static_cast<float>(1)) <=
        numeric_limits<float>::epsilon());
    
    BOOST_TEST(abs(sinhc_pi<double>(static_cast<double>(0))-
        static_cast<double>(1)) <=
        numeric_limits<double>::epsilon());
    
    BOOST_TEST(abs(sinhc_pi<long double>(static_cast<long double>(0))-
        static_cast<long double>(1)) <=
        numeric_limits<long double>::epsilon());
    
    
    return(::boost::exit_success);
}
