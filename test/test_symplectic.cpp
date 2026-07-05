/*
 * Copyright Jacob Hass, 2026
 * Use, modification and distribution are subject to the
 * Boost Software License, Version 1.0. (See accompanying file
 * LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
 */
#define BOOST_TEST_MODULE symplectic_quadrature

#include <algorithm>
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

// Equations of motion for simple harmonic oscillator
template <class Real>
std::vector<Real> vector_dHdp(std::vector<Real> p) { return p; }
template <class Real>
std::vector<Real> vector_dHdq(std::vector<Real> q) { return q; }

template <class RealType>
void test_invalid_parameters()
{
    std::vector<RealType> q0 = {1};
    std::vector<RealType> p0 = {0};
    // Negative timestep
    BOOST_CHECK_THROW(boost::math::quadrature::integrate_hamiltonian(q0, p0, -0.1, 10, vector_dHdp<RealType>, vector_dHdq<RealType>), std::domain_error);

    // Method not in {'Y6', 'Y4', 'Y2'}
    BOOST_CHECK_THROW(boost::math::quadrature::integrate_hamiltonian(q0, p0, 0.1, 10, vector_dHdp<RealType>, vector_dHdq<RealType>, "InvalidMethod"), std::out_of_range);
}

/* Test if SHO energy fluctuations are below a given tolerance*/
template <class RealType>
void test_harmonic_oscillator(const RealType tol, const std::string method)
{
    BOOST_MATH_STD_USING

    RealType dt = 0.05;
    RealType t_end = 100;
    unsigned int steps = t_end / dt;

    std::vector<RealType> q0 = {1};
    std::vector<RealType> p0 = {0};

    std::vector<std::vector<RealType> > p;
    std::vector<std::vector<RealType> > q;

    std::tie(p, q) = boost::math::quadrature::integrate_hamiltonian(p0, q0, dt, steps, vector_dHdp<RealType>, vector_dHdp<RealType>, method);

    RealType p_val;
    RealType q_val;
    std::vector<RealType> abs_energy_error(p.size());
    for (unsigned i=0; i < p.size(); i++)
    {
        p_val = p[i][0];
        q_val = q[i][0];

        abs_energy_error[i] = std::abs(std::pow(p_val, 2) + std::pow(q_val, 2) - 1);
    }

    RealType max_error = *std::max_element(std::begin(abs_energy_error), std::end(abs_energy_error));
    BOOST_CHECK_LE(max_error, tol);
}

BOOST_AUTO_TEST_CASE(symplectic_quadrature)
{
    test_invalid_parameters<double>();
    test_harmonic_oscillator<double>(1e-10, "Y6");
    test_harmonic_oscillator<double>(7e-4, "Y4");
    test_harmonic_oscillator<double>(7e-4, "Y2");
}
