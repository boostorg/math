//  (C) Copyright Nick Thompson 2018.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_TOOLS_BIVARIATE_STATISTICS_HPP
#define BOOST_MATH_TOOLS_BIVARIATE_STATISTICS_HPP

#include <iterator>
#include <tuple>
#include <boost/type_traits.hpp>
#include <boost/assert.hpp>
#include <boost/multiprecision/detail/number_base.hpp>


namespace boost{ namespace math{ namespace tools {

template<class Container>
auto
means_and_population_covariance(Container const & u, Container const & v)
{
    using Real = typename Container::value_type;
    using std::size;
    BOOST_ASSERT_MSG(size(u) == size(v), "The size of each vector must be the same to compute covariance.");
    BOOST_ASSERT_MSG(size(u) > 0, "Computing covariance requires at least one sample.");

    // See Equation III.9 of "Numerically Stable, Single-Pass, Parallel Statistics Algorithms", Bennet et al.
    Real cov = 0;
    Real mu_u = u[0];
    Real mu_v = v[0];

    for(size_t i = 1; i < size(u); ++i)
    {
        Real u_tmp = (u[i] - mu_u)/(i+1);
        Real v_tmp = v[i] - mu_v;
        cov += i*u_tmp*v_tmp;
        mu_u = mu_u + u_tmp;
        mu_v = mu_v + v_tmp/(i+1);
    }

    return std::make_tuple(mu_u, mu_v, cov/size(u));
}

template<class Container>
auto
population_covariance(Container const & u, Container const & v)
{
    auto [mu_u, mu_v, cov] = boost::math::tools::means_and_population_covariance(u, v);
    return cov;
}

}}}
#endif
