//  (C) Copyright John Maddock 2024.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include "math_unit_test.hpp"
#include <cfenv>
#include <iostream>
#include <boost/math/special_functions/gamma.hpp>


int main()
{
   CHECK_EQUAL(boost::math::tgamma(-200.5), 0.0); // triggers internal exception handling
   CHECK_EQUAL(boost::math::gamma_p(500.125, 1e-50), 0.0); // triggers internal exception handling
   return boost::math::test::report_errors();
}
