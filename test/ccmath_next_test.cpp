//  (C) Copyright John Maddock 2008 - 2022.
//  (C) Copyright Matt Borland 2022.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <iostream>
#include <iomanip>
#include <limits>
#include <boost/math/ccmath/next.hpp>

template <typename T>
void test_next()
{

}

int main(void)
{
    test_next<float>();
    test_next<double>();

    #ifndef BOOST_MATH_NO_LONG_DOUBLE_MATH_FUNCTIONS
    test_next<long double>();
    #endif
}
