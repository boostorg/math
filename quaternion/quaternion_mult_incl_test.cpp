// test file for quaternion.hpp

//  (C) Copyright Hubert Holin 2001.
//  Distributed under the Boost Software License, Version 1.0. (See
//  accompanying file LICENSE_1_0.txt or copy at
//  http://www.boost.org/LICENSE_1_0.txt)


#include <boost/test/unit_test.hpp>
#include <boost/test/unit_test_suite_ex.hpp>

#include "quaternion_mi1.h"
#include "quaternion_mi2.h"


boost::unit_test_framework::test_suite *    init_unit_test_suite(int, char *[])
{
    ::boost::unit_test_framework::unit_test_log::instance().
            set_log_threshold_level_by_name("messages");
    
    boost::unit_test_framework::test_suite *    test =
        BOOST_TEST_SUITE("quaternion_multiple_inclusion_test");
    
    test->add(BOOST_TEST_CASE(&quaternion_mi1));
    test->add(BOOST_TEST_CASE(&quaternion_mi2));
    
    return(test);
}

