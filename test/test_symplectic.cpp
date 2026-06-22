/*
 * Copyright Nick Thompson, 2017
 * Use, modification and distribution are subject to the
 * Boost Software License, Version 1.0. (See accompanying file
 * LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
 */
#define BOOST_TEST_MODULE symplectic_quadrature

#include <valarray>
#include <boost/math/tools/config.hpp>
#include <boost/test/included/unit_test.hpp>
#include <boost/test/tools/floating_point_comparison.hpp>
#include <boost/math/tools/test_value.hpp>
#include <boost/math/concepts/real_concept.hpp>
#include <boost/math/quadrature/symplectic.hpp>
#ifdef BOOST_HAS_FLOAT128
#include <boost/multiprecision/complex128.hpp>
#endif

#if __has_include(<stdfloat>)
#  include <stdfloat>
#endif

using boost::math::quadrature::integrate_hamiltonian;

template <class Real>
Real singleton_dHdp(Real p) { return p; }
template <class Real>
Real singleton_dHdq(Real q) { return q; }

template <class Real>
std::valarray<Real> vector_dHdp(std::valarray<Real> p) { return p; }
template <class Real>
std::valarray<Real> vector_dHdq(std::valarray<Real> q) { return q; }

template<class RealType>
void test_input_singleton_vs_vector()
{
    RealType dt = 0.05;
    unsigned int steps = 3;
    RealType q0 = 1;
    RealType p0 = 0;


    std::vector<RealType> p;
    std::vector<RealType> q;

    std::tie(p, q) = boost::math::quadrature::integrate_hamiltonian(p0, q0, dt, steps, singleton_dHdp<RealType>, singleton_dHdq<RealType>);

    std::valarray<RealType> q0_vector = {1};
    std::valarray<RealType> p0_vector = {0};
    std::vector<std::valarray<RealType> > p_vector;
    std::vector<std::valarray<RealType> > q_vector;

    std::tie(p_vector, q_vector) = boost::math::quadrature::integrate_hamiltonian(p0_vector, q0_vector, dt, steps, vector_dHdp<RealType>, vector_dHdq<RealType>);
    
    BOOST_CHECK_EQUAL(p[0], p_vector[0][0]);
    BOOST_CHECK_EQUAL(p[1], p_vector[1][0]);
    BOOST_CHECK_EQUAL(p[2], p_vector[2][0]);

}

BOOST_AUTO_TEST_CASE(symplectic_quadrature)
{
    test_input_singleton_vs_vector<double>();
}
