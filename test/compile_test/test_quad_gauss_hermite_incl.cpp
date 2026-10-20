//  Copyright Nick Thompson 2026.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// Basic sanity check that header
// #includes all the files that it needs to.
//
#include <boost/math/quadrature/gauss_hermite.hpp>
//
// Note this header includes no other headers, this is
// important if this test is to be meaningful:
//
#include "test_compile_result.hpp"

void compile_and_link_test()
{
    auto integrand = [](double x) { return x; };
    check_result<double>(boost::math::quadrature::gauss_hermite<double, 7>::integrate(integrand));
    check_result<double>(boost::math::quadrature::gauss_hermite<double, 4>::integrate(integrand));
}
