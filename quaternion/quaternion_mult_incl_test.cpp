// test file for quaternion.hpp

//  (C) Copyright Hubert Holin 2001. Permission to copy, use, modify, sell and
//  distribute this software is granted provided this copyright notice appears
//  in all copies. This software is provided "as is" without express or implied
//  warranty, and with no claim as to its suitability for any purpose.


#define BOOST_INCLUDE_MAIN  // for testing, include rather than link
#include <boost/test/test_tools.hpp>


#include <boost/config.hpp>

#include "quaternion_mi1.h"
#include "quaternion_mi2.h"

int    test_main(int, char *[])

{
    
    quaternion_mi1();
    
    quaternion_mi2();
    
    return(::boost::exit_success);
}
