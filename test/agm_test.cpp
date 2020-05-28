/*
 * Copyright Nick Thompson, 2019
 * Use, modification and distribution are subject to the
 * Boost Software License, Version 1.0. (See accompanying file
 * LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
 */

#include "math_unit_test.hpp"
#include <numeric>
#include <utility>
#include <random>
#include <boost/lexical_cast.hpp>
#include <boost/core/demangle.hpp>
#include <boost/math/tools/agm.hpp>
#ifdef BOOST_HAS_FLOAT128
#include <boost/multiprecision/float128.hpp>
using boost::multiprecision::float128;
#endif

using boost::math::tools::agm;

template<class Real>
void test_gauss_constant()
{
    // http://oeis.org/A014549/constant
    Real G_expected = boost::lexical_cast<Real>(".83462684167407318628142973279904680899399301349034700244982737010368199270952641186969116035127532412906785");

    Real G_computed = 1/agm(sqrt(Real(2)), Real(1));
    if(!CHECK_ULP_CLOSE(G_expected, G_computed, 2)) {
        std::cerr << "  Gauss constant not computed correctly.\n";
    }
}

template<typename Real>
void test_scaling()
{
    Real a = 2;
    Real g = 1;
    Real scale = 7;
    Real expected = agm(scale*a, scale*g);
    Real computed = scale*agm(a, g);
    if(!CHECK_ULP_CLOSE(expected, computed, 2)) {
        std::cerr << "  Scaling property agm(kx,ky) = k*agm(x, y) is violated.\n";
    }
}



int main()
{
    test_gauss_constant<float>();
    test_gauss_constant<double>();
    test_gauss_constant<long double>();
    return boost::math::test::report_errors();
}
