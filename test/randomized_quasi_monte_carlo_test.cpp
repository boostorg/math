/*
 *  (C) Copyright Nick Thompson 2018.
 *  Use, modification and distribution are subject to the
 *  Boost Software License, Version 1.0. (See accompanying file
 *  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
 *
 */

#define BOOST_TEST_MODULE quasi_monte_carlo_test
#include <iostream>
#include <boost/math/constants/constants.hpp>
#include <boost/test/included/unit_test.hpp>
#include <boost/math/quadrature/randomized_quasi_monte_carlo.hpp>

using boost::math::constants::pi;
using boost::math::quadrature::randomized_quasi_monte_carlo;

template<class Real>
void test_pi()
{
    std::cout << "Testing pi is calculated correctly using quasi-Monte-Carlo\n";
    auto g = [](std::vector<Real> const & x)->Real
    {
        Real r = x[0]*x[0]+x[1]*x[1];
        if (r <= 1)
        {
            return 4;
        }
        return 0;
    };

    std::vector<std::pair<Real, Real>> bounds{{0, 1}, {0, 1}};
    randomized_quasi_monte_carlo<Real, decltype(g)> mc(g, bounds, (Real) 0.0005);

    auto task = mc.integrate();
    Real pi_estimated = task.get();
    if (abs(pi_estimated - pi<Real>())/pi<Real>() > 0.005)
    {
        std::cout << "Error in estimation of pi too high, function calls: " << mc.calls() << "\n";
        BOOST_CHECK_CLOSE_FRACTION(pi_estimated, pi<Real>(), 0.005);
    }

}

template<class Real>
void test_constant()
{
    std::cout << "Testing constants are integrated correctly using randomized quasi-Monte-Carlo\n";
    auto g = [](std::vector<Real> const & x)->Real
    {
      return 1;
    };

    std::vector<std::pair<Real, Real>> bounds{{0, 1}, {0, 1}};
    randomized_quasi_monte_carlo<Real, decltype(g)> mc(g, bounds, (Real) 0.0001);

    auto task = mc.integrate();
    Real one = task.get();
    BOOST_CHECK_CLOSE_FRACTION(one, 1, 0.001);
    BOOST_CHECK_SMALL(mc.current_error_estimate(), std::numeric_limits<Real>::epsilon());
    BOOST_CHECK(mc.calls() > 1000);
}

template<class Real>
void test_nan()
{
    std::cout << "Testing that a reasonable action is performed by the randomized quasi-Monte-Carlo integrator when singularities at hit.\n";
    auto g = [](std::vector<Real> const & x)->Real
    {
      return (Real) 1/ (Real) 0;
    };

    std::vector<std::pair<Real, Real>> bounds{{0, 1}, {0, 1}};
    randomized_quasi_monte_carlo<Real, decltype(g)> mc(g, bounds, (Real) 0.0001);

    auto task = mc.integrate();
    Real result = task.get();
    // I think this is reasonable, but should it throw an exception?
    BOOST_CHECK(std::isnan(result));
}


template<class Real>
void test_exception_from_integrand()
{
    std::cout << "Testing that a reasonable action is performed by the quasi-Monte-Carlo integrator when the integrand throws an exception.\n";
    auto g = [](std::vector<Real> const & x)->Real
    {
        if (x[0] > 0.5 && x[0] < 0.5001)
        {
            throw std::domain_error("You have done something wrong.\n");
        }
        return (Real) 1;
    };

    std::vector<std::pair<Real, Real>> bounds{{0, 1}, {0, 1}};
    randomized_quasi_monte_carlo<Real, decltype(g)> mc(g, bounds, (Real) 0.0001);

    auto task = mc.integrate();
    bool caught_exception = false;
    try
    {
      Real result = task.get();
      // Get rid of unused variable warning:
      std::ostream cnull(0);
      cnull << result;
    }
    catch(std::exception const & e)
    {
        caught_exception = true;
    }
    BOOST_CHECK(caught_exception);
}


template<class Real, size_t dimension>
void test_product()
{
    std::cout << "Testing that product functions are integrated correctly by quasi-Monte-Carlo\n";
    auto g = [&](std::vector<Real> const & x)->Real
    {
        Real y = 1;
        for (size_t i = 0; i < x.size(); ++i)
        {
            y *= 2*x[i];
        }
        return y;
    };

    std::vector<std::pair<Real, Real>> bounds(dimension);
    for (size_t i = 0; i < dimension; ++i)
    {
        bounds[i] = std::make_pair<Real, Real>(0, 1);
    }
    randomized_quasi_monte_carlo<Real, decltype(g)> mc(g, bounds, (Real) 0.001);

    auto task = mc.integrate();
    Real y = task.get();
    BOOST_CHECK_CLOSE_FRACTION(y, 1, 0.01);
}

template<class Real>
void test_upper_bound_infinite()
{
    std::cout << "Testing that infinite upper bounds are integrated correctly by quasi-Monte-Carlo\n";
    auto g = [](std::vector<Real> const & x)->Real
    {
        return 1.0/(x[0]*x[0] + 1.0);
    };

    std::vector<std::pair<Real, Real>> bounds(1);
    for (size_t i = 0; i < bounds.size(); ++i)
    {
        bounds[i] = std::make_pair<Real, Real>(0, std::numeric_limits<Real>::infinity());
    }
    randomized_quasi_monte_carlo<Real, decltype(g)> mc(g, bounds, (Real) 0.001);

    auto task = mc.integrate();
    Real y = task.get();
    BOOST_CHECK_CLOSE_FRACTION(y, M_PI/2, 0.01);
}

template<class Real>
void test_lower_bound_infinite()
{
    std::cout << "Testing that infinite lower bounds are integrated correctly by randomized quasi-Monte-Carlo\n";
    auto g = [](std::vector<Real> const & x)->Real
    {
        return 1.0/(x[0]*x[0] + 1.0);
    };

    std::vector<std::pair<Real, Real>> bounds(1);
    for (size_t i = 0; i < bounds.size(); ++i)
    {
        bounds[i] = std::make_pair<Real, Real>(-std::numeric_limits<Real>::infinity(), 0);
    }
    randomized_quasi_monte_carlo<Real, decltype(g)> mc(g, bounds, (Real) 0.001);

    auto task = mc.integrate();
    Real y = task.get();
    BOOST_CHECK_CLOSE_FRACTION(y, M_PI/2, 0.01);
}

template<class Real>
void test_double_infinite()
{
    std::cout << "Testing that double infinite bounds are integrated correctly by randomized quasi-Monte-Carlo\n";
    auto g = [](std::vector<Real> const & x)->Real
    {
        return 1.0/(x[0]*x[0] + 1.0);
    };

    std::vector<std::pair<Real, Real>> bounds(1);
    for (size_t i = 0; i < bounds.size(); ++i)
    {
        bounds[i] = std::make_pair<Real, Real>(-std::numeric_limits<Real>::infinity(), std::numeric_limits<Real>::infinity());
    }
    randomized_quasi_monte_carlo<Real, decltype(g)> mc(g, bounds, (Real) 0.001);

    auto task = mc.integrate();
    Real y = task.get();
    BOOST_CHECK_CLOSE_FRACTION(y, M_PI, 0.01);
}

template<class Real, size_t dimension>
void test_radovic()
{
    // See: Generalized Halton Sequences in 2008: A Comparative Study, function g1:
    auto g = [](std::vector<Real> const & x)->Real
    {
        using std::abs;
        Real alpha = 0.01;
        Real z = 1;
        for (size_t i = 0; i < dimension; ++i)
        {
            z *= (abs(4*x[i]-2) + alpha)/(1+alpha);
        }
        return z;
    };

    std::vector<std::pair<Real, Real>> bounds(dimension);
    for (size_t i = 0; i < bounds.size(); ++i)
    {
        bounds[i] = std::make_pair<Real, Real>(0, 1);
    }
    randomized_quasi_monte_carlo<Real, decltype(g)> rqmc(g, bounds, (Real) 0.001);

    auto task = rqmc.integrate();
    Real y = task.get();
    BOOST_CHECK_CLOSE_FRACTION(y, 1, 0.01);
}


BOOST_AUTO_TEST_CASE(randomized_quasi_monte_carlo_test)
{
    test_nan<float>();
    test_pi<float>();
    test_pi<double>();
    test_pi<long double>();
    test_constant<float>();
    test_constant<double>();
    test_constant<long double>();
    test_exception_from_integrand<float>();
    test_product<float, 1>();
    test_product<float, 2>();
    test_product<float, 3>();
    test_product<float, 4>();
    test_product<float, 5>();
    test_product<float, 6>();
    test_product<double, 1>();
    test_product<double, 2>();
    test_product<double, 3>();
    test_product<double, 4>();
    test_upper_bound_infinite<float>();
    test_upper_bound_infinite<double>();
    test_lower_bound_infinite<float>();
    test_lower_bound_infinite<double>();
    test_double_infinite<float>();
    test_double_infinite<double>();
    test_radovic<float, 1>();
    test_radovic<float, 2>();
    test_radovic<float, 3>();
    test_radovic<double, 1>();
    test_radovic<double, 2>();
    test_radovic<double, 3>();
    test_radovic<double, 4>();
    test_radovic<double, 5>();
}
