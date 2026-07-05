//  Copyright Matt Borland 2026
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// Regression test for https://github.com/boostorg/math/issues/1412
//
// Part 2: Test all of #include <boost/math/special_functions.hpp> as that is what libc++ is going to use

#define BOOST_MATH_NO_EXCEPTIONS

// Claim the <iomanip> include guard so any later #include <iomanip> expands to nothing and std::setprecision is unavailable.
// The original reproducer with _SSTREAM broke <complex> which is outside our purview.
#define _LIBCPP_IOMANIP        // libc++
#define _GLIBCXX_IOMANIP 1     // libstdc++

#include "compile_test/instantiate.hpp"

int main(int argc, char* [])
{
    if (argc > 10000)
    {
        instantiate(double(0));
        instantiate_mixed(double(0));
    }
}
