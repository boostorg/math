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
#include <boost/math/concepts/real_concept.hpp>
#include <boost/math/tools/agm.hpp>
#ifdef BOOST_HAS_FLOAT128
#include <boost/multiprecision/float128.hpp>
using boost::multiprecision::float128;
#endif

using boost::math::tools::agm;

template<class Real, class Policy>
void test_gauss_constant(const Policy &pol)
{
    // http://oeis.org/A014549/constant
    Real G_expected = boost::lexical_cast<Real>(".83462684167407318628142973279904680899399301349034700244982737010368199270952641186969116035127532412906785");

    Real G_computed = 1/agm(sqrt(Real(2)), Real(1), pol);
    if(!CHECK_ULP_CLOSE(G_expected, G_computed, 2)) {
        std::cerr << "  Gauss constant not computed correctly.\n";
    }
}

template<typename Real, class Policy>
void test_scaling(const Policy &pol)
{
    Real a = 2;
    Real g = 1;
    Real scale = 7;
    Real expected = agm(scale*a, scale*g, pol);
    Real computed = scale*agm(a, g, pol);
    if(!CHECK_ULP_CLOSE(expected, computed, 2)) {
        std::cerr << "  Scaling property agm(kx,ky) = k*agm(x, y) is violated.\n";
    }

    expected = 0;
    computed = agm(a, Real(0), pol);
    if(!CHECK_ULP_CLOSE(expected, computed, 0)) {
        std::cerr << "  agm(a, 0) != 0.\n";
    }

    computed = agm(Real(0), Real(0), pol);
    if(!CHECK_ULP_CLOSE(expected, computed, 0)) {
        std::cerr << "  agm(0, 0) != 0.\n";
    }

    expected = 1;
    computed = agm(Real(1), Real(1), pol);
    if(!CHECK_ULP_CLOSE(expected, computed, 0)) {
        std::cerr << "  agm(1, 1) != 1.\n";
    }

    expected = 7;
    computed = agm(Real(7), Real(7), pol);
    if(!CHECK_ULP_CLOSE(expected, computed, 0)) {
        std::cerr << "  agm(7, 7) != 7.\n";
    }

    // Properties I found at: https://mathworld.wolfram.com/Arithmetic-GeometricMean.html
    // agm(x,y) = agm((x+y)/2, sqrt(xy))
    expected = agm(Real(3), Real(1), pol);
    computed = agm(Real(2), sqrt(Real(3)), pol);
    if(!CHECK_ULP_CLOSE(expected, computed, 0)) {
        std::cerr << "  agm(x, y) != agm((x+y)/2, sqrt(xy)).\n";
    }

    //computed = agm(std::numeric_limits<Real>::infinity(), Real(7));
    //std::cout << "Computed at infinity = " << computed << "\n";

    for (Real x = 0; x < 1; x += Real(1)/128) {
        expected = agm(Real(1), sqrt(1-x*x), pol);
        computed = agm(1+x, 1-x, pol);
        if(!CHECK_ULP_CLOSE(expected, computed, 0)) {
            std::cerr << "  agm(1, sqrt(1-x^2) != agm(1+x,1-x).\n";
        }
    }
}


int main()
{
    typedef boost::math::policies::policy<> pol;
    pol test_policy = pol();

    test_gauss_constant<float, pol>(test_policy);
    test_gauss_constant<double, pol>(test_policy);
    test_gauss_constant<long double, pol>(test_policy);

    test_scaling<float, pol>(test_policy);
    test_scaling<double, pol>(test_policy);
    test_scaling<long double, pol>(test_policy);
    test_scaling<boost::math::concepts::real_concept, pol>(test_policy);

    #ifdef BOOST_HAS_FLOAT128
    test_gauss_constant<float128, pol>(test_policy);
    test_scaling<float128, pol>(test_policy);
    #endif
    return boost::math::test::report_errors();
}
