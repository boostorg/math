/*
 * Copyright Nick Thompson, 2017
 * Use, modification and distribution are subject to the
 * Boost Software License, Version 1.0. (See accompanying file
 * LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
 */
#define BOOST_TEST_MODULE catmull_rom_test

#include <array>
#include <boost/cstdfloat.hpp>
#include <boost/type_index.hpp>
#include <boost/test/included/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/math/interpolators/catmull_rom.hpp>
#include <boost/multiprecision/cpp_bin_float.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>

using boost::multiprecision::cpp_bin_float_50;
using boost::math::catmull_rom;

template<class Real>
void test_alpha_distance()
{
    Real tol = std::numeric_limits<Real>::epsilon();
    std::array<Real, 3> v1 = {0,0,0};
    std::array<Real, 3> v2 = {1,0,0};
    Real alpha = 0.5;
    Real d = boost::math::detail::alpha_distance<Real, std::array<Real, 3>, 3>(v1, v2, alpha);
    BOOST_CHECK_CLOSE_FRACTION(d, 1, tol);

    d = boost::math::detail::alpha_distance<Real, std::array<Real, 3>, 3>(v1, v2, 0);
    BOOST_CHECK_CLOSE_FRACTION(d, 1, tol);

    d = boost::math::detail::alpha_distance<Real, std::array<Real, 3>, 3>(v1, v2, 1);
    BOOST_CHECK_CLOSE_FRACTION(d, 1, tol);


    v2[0] = 2;
    d = boost::math::detail::alpha_distance<Real, std::array<Real, 3>, 3>(v1, v2, alpha);
    BOOST_CHECK_CLOSE_FRACTION(d, pow(2, (Real)1/ (Real) 2), tol);

    d = boost::math::detail::alpha_distance<Real, std::array<Real, 3>, 3>(v1, v2, 0);
    BOOST_CHECK_CLOSE_FRACTION(d, 1, tol);

    d = boost::math::detail::alpha_distance<Real, std::array<Real, 3>, 3>(v1, v2, 1);
    BOOST_CHECK_CLOSE_FRACTION(d, 2, tol);
}


template<class Real>
void test_linear()
{
    Real tol = 10*std::numeric_limits<Real>::epsilon();
    std::vector<std::array<Real, 3>> v(4);
    v[0] = {0,0,0};
    v[1] = {1,0,0};
    v[2] = {2,0,0};
    v[3] = {3,0,0};
    catmull_rom<Real, std::array<Real, 3>, 3> cr(v.data(), v.size());

    // Test that the interpolation condition is obeyed:
    BOOST_CHECK_CLOSE_FRACTION(cr.max_parameter(), 3, tol);
    auto p0 = cr(0.0);
    BOOST_CHECK_SMALL(p0[0], tol);
    BOOST_CHECK_SMALL(p0[1], tol);
    BOOST_CHECK_SMALL(p0[2], tol);
    auto p1 = cr(1.0);
    BOOST_CHECK_CLOSE_FRACTION(p1[0], 1, tol);
    BOOST_CHECK_SMALL(p1[1], tol);
    BOOST_CHECK_SMALL(p1[2], tol);

    auto p2 = cr(2.0);
    BOOST_CHECK_CLOSE_FRACTION(p2[0], 2, tol);
    BOOST_CHECK_SMALL(p2[1], tol);
    BOOST_CHECK_SMALL(p2[2], tol);


    auto p3 = cr(3.0);
    BOOST_CHECK_CLOSE_FRACTION(p3[0], 3, tol);
    BOOST_CHECK_SMALL(p3[1], tol);
    BOOST_CHECK_SMALL(p3[2], tol);

    Real s = cr.parameter_at_point(0);
    BOOST_CHECK_SMALL(s, tol);

    s = cr.parameter_at_point(1);
    BOOST_CHECK_CLOSE_FRACTION(s, 1, tol);

    s = cr.parameter_at_point(2);
    BOOST_CHECK_CLOSE_FRACTION(s, 2, tol);

    s = cr.parameter_at_point(3);
    BOOST_CHECK_CLOSE_FRACTION(s, 3, tol);

    // Test that the function is linear on the interval [1,2]:
    for (double s = 1; s < 2; s += 0.01)
    {
        auto p = cr(s);
        BOOST_CHECK_CLOSE_FRACTION(p[0], s, tol);
        BOOST_CHECK_SMALL(p[1], tol);
        BOOST_CHECK_SMALL(p[2], tol);

        auto tangent = cr.prime(s);
        // TODO: Fix the tangent vector!
        //BOOST_CHECK_CLOSE_FRACTION(p[0], 1, tol);
        BOOST_CHECK_SMALL(p[1], tol);
        BOOST_CHECK_SMALL(p[2], tol);
    }

}

template<class Real>
void test_affine_invariance()
{
    std::cout << "Testing that the Catmull-Rom spline is affine invariant.\n";

}

template<class Real>
void test_data_representations()
{
    std::cout << "Testing that the Catmull-Rom spline works with multiple data representations.\n";
}

BOOST_AUTO_TEST_CASE(catmull_rom_test)
{
    test_alpha_distance<float>();
    test_alpha_distance<double>();
    test_alpha_distance<long double>();
    test_alpha_distance<cpp_bin_float_50>();

    test_linear<float>();
    test_linear<double>();
    test_linear<long double>();
    test_linear<cpp_bin_float_50>();
    /*
    test_affine_invariance<float>();
    test_data_representations<float>();*/
}
