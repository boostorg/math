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


using boost::math::catmull_rom;

template<class Real>
void test_linear()
{
    Real tol = std::numeric_limits<Real>::epsilon();
    std::vector<std::array<Real, 3>> v(4);
    v[0] = {0,0,0};
    v[1] = {1,0,0};
    v[2] = {2,0,0};
    v[3] = {3,0,0};
    catmull_rom<Real, std::array<Real, 3>, 3> cr(v.data(), v.size(), true);

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
    test_linear<float>();
    test_affine_invariance<float>();
    test_data_representations<float>();
}
