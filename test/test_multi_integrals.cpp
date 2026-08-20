/*
 * Copyright Jacob Hass, 2026
 * Use, modification and distribution are subject to the
 * Boost Software License, Version 1.0. (See accompanying file
 * LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
 */

#define BOOST_TEST_MODULE multi_integrals

#include <algorithm>
#include <type_traits>
#include <vector>
#include <boost/math/tools/config.hpp>
#include <boost/test/included/unit_test.hpp>
#include <boost/test/tools/floating_point_comparison.hpp>
#include <boost/math/tools/test_value.hpp>
#include <boost/math/concepts/real_concept.hpp>
#include <boost/math/quadrature/multi_integrals.hpp>
#include <boost/math/quadrature/tanh_sinh.hpp>
#include <boost/math/quadrature/trapezoidal.hpp>
#include <boost/math/quadrature/gauss_kronrod.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/multiprecision/cpp_bin_float.hpp>

#if __has_include(<stdfloat>)
#  include <stdfloat>
#endif

using boost::math::quadrature::integrateND;
using boost::math::quadrature::tanh_sinh;
using boost::math::quadrature::gauss_kronrod;

template <typename RealType>
void test_catalan()
{
    auto f = [](const std::vector<RealType>& coords) {
        return 1 / (1 + coords[0] * coords[0] * coords[1] * coords[1]);
    };

    std::vector<RealType> lower_bounds = {0.0, 0.0};
    std::vector<RealType> upper_bounds = {1.0, 1.0}; 
    
    typedef std::function<RealType(RealType)> func;
    
    auto integrate = [](func f, RealType a, RealType b)
    {
        tanh_sinh<RealType> integrator;
        return integrator.integrate(f, a, b);
    }; 

    RealType result = integrateND(f, lower_bounds, upper_bounds, integrate);
    RealType tol = boost::math::tools::epsilon<RealType>();
    BOOST_CHECK_CLOSE_FRACTION(result, boost::math::constants::catalan<RealType>(), 4 * tol);
}

template <typename RealType>
void test_euler()
{
    BOOST_MATH_STD_USING

    auto f = [](const std::vector<RealType>& coords) {
        RealType x = coords[0];
        RealType y = coords[1];
        return exp( -x * y ) * (1 / (1 + x) - exp(-x));
    };

    std::vector<RealType> lower_bounds = {0.0, 0.0};
    std::vector<RealType> upper_bounds = {std::numeric_limits<RealType>::infinity(), std::numeric_limits<RealType>::infinity()}; 
    
    typedef std::function<RealType(RealType)> func;
    
    auto integrate = [](func f, RealType a, RealType b)
    {
        return gauss_kronrod<RealType, 61>::integrate(f, a, b);
    };

    RealType result = integrateND(f, lower_bounds, upper_bounds, integrate);
    RealType tol = boost::math::tools::epsilon<RealType>();
    // Long double precision is not good compared to other precisions
    RealType check_tol = std::is_same<long double, RealType>::value ? 3000 * tol : 4 * tol;
    BOOST_CHECK_CLOSE_FRACTION(result, boost::math::constants::euler<RealType>(), check_tol);
}

template <typename RealType>
void test_square()
{
    BOOST_MATH_STD_USING 

    auto f = [](const std::vector<RealType>& coords)
    {
        return abs(coords[0] + coords[1] + coords[2]);
    };

    typedef std::function<RealType(RealType)> func;
    auto integrate = [](func f, RealType a, RealType b)
    {
        return gauss_kronrod<RealType, 61>::integrate(f, a, b);
    }; 

    std::vector<RealType> lower_bounds = {-1.0, -1.0, -1.0};
    std::vector<RealType> upper_bounds = {1.0, 1.0, 1.0};

    RealType result = integrateND(f, lower_bounds, upper_bounds, integrate);
    RealType tol = boost::math::tools::epsilon<RealType>();
    RealType multFactor;

    // Does not seem to perform well on this integral
    if(std::is_same<RealType, float>::value){
        multFactor = 200;
    }
    else if (std::is_same<RealType, double>::value){
        multFactor = 4 * pow(10, 6);
    }
    else if (std::is_same<RealType, long double>::value){
        multFactor = 2 * pow(10, 8);
    }
    
    BOOST_CHECK_CLOSE_FRACTION(result, RealType(13) / RealType(2), multFactor * tol);
}

BOOST_AUTO_TEST_CASE(multi_integrals)
{
    // Catalan constant tests
    test_catalan<float>();
    test_catalan<double>();
    #ifndef BOOST_MATH_NO_LONG_DOUBLE_MATH_FUNCTIONS
        test_catalan<long double>();
    #endif
    #ifndef BOOST_MATH_NO_REAL_CONCEPT_TESTS
        test_catalan<boost::math::concepts::real_concept>();
    #endif

    // Euler constant tests
    test_euler<float>();
    test_euler<double>();
    #ifndef BOOST_MATH_NO_LONG_DOUBLE_MATH_FUNCTIONS
        test_euler<long double>();
    #endif

    // Square integral tests
    test_square<float>();
    test_square<double>();
    #ifndef BOOST_MATH_NO_LONG_DOUBLE_MATH_FUNCTIONS
        test_square<long double>();
    #endif
}