//  Copyright Matt Borland 2026
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// Regression test for https://github.com/boostorg/math/issues/1412

#define BOOST_MATH_NO_EXCEPTIONS

// Claim the <iomanip> include guard so any later #include <iomanip> expands to nothing and std::setprecision is unavailable.
// The original reproducer with _SSTREAM broke <complex> which iso outside our purview
#define _LIBCPP_IOMANIP        // libc++
#define _GLIBCXX_IOMANIP 1     // libstdc++

#include <boost/math/policies/error_handling.hpp>
#include <boost/math/special_functions/laguerre.hpp>

int main()
{
    // Exercise a special function so the no-exceptions header path is
    // instantiated, not just parsed.
    return static_cast<int>(boost::math::laguerre(1u, 0u, 0.5));
}
