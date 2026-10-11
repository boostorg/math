//  (C) Copyright Matt Borland 2026.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// A std-free module consumer: it names no standard library entity, so it builds
// on toolchains (GCC/libstdc++) that cannot yet mix textual std includes with
// import std in a consumer. Instantiating these functions makes the compiler
// resolve the namespace-scope tables/constants they use, which is where GCC
// enforces the module exposure rule that clang does not diagnose.

#ifndef BOOST_MATH_BUILD_MODULE
#include <boost/math/special_functions/lambert_w.hpp>
#else
import boost.math;
#endif

int main()
{
    int status {0};

    // lambert_w0/wm1 index the namespace-scope lookup tables at run time; that
    // odr-use is what makes GCC reject the tables if they carry internal linkage.
    const double w0 {boost::math::lambert_w0(1.0)};
    const double wm1 {boost::math::lambert_wm1(-0.1)};
    status += (w0 > 0.5 && w0 < 0.6) ? 0 : 1;
    status += (wm1 < -3.0) ? 0 : 1;

    return status;
}
