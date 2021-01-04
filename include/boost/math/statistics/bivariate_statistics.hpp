//  (C) Copyright Nick Thompson 2018.
//  (C) Copyright Matt Borland 2021.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_STATISTICS_BIVARIATE_STATISTICS_HPP
#define BOOST_MATH_STATISTICS_BIVARIATE_STATISTICS_HPP

#include <iterator>
#include <tuple>
#include <type_traits>
#include <cmath>
#include <cstddef>
#include <boost/assert.hpp>


namespace boost{ namespace math{ namespace statistics { namespace detail {

template<class ReturnType, class Container>
ReturnType means_and_covariance_impl(Container const & u, Container const & v)
{
    using Real = typename std::tuple_element<0, ReturnType>::type;

    BOOST_ASSERT_MSG(u.size() == v.size(), "The size of each vector must be the same to compute covariance.");
    BOOST_ASSERT_MSG(u.size() > 0, "Computing covariance requires at least one sample.");

    // See Equation III.9 of "Numerically Stable, Single-Pass, Parallel Statistics Algorithms", Bennet et al.
    Real cov = 0;
    Real mu_u = u[0];
    Real mu_v = v[0];

    for(std::size_t i = 1; i < u.size(); ++i)
    {
        Real u_tmp = (u[i] - mu_u)/(i+1);
        Real v_tmp = v[i] - mu_v;
        cov += i*u_tmp*v_tmp;
        mu_u = mu_u + u_tmp;
        mu_v = mu_v + v_tmp/(i+1);
    }

    return std::make_tuple(mu_u, mu_v, cov/u.size());
}

template<class ReturnType, class Container>
ReturnType correlation_coefficient_impl(Container const & u, Container const & v)
{
    using Real = ReturnType;
    using std::sqrt;
    BOOST_ASSERT_MSG(u.size() == v.size(), "The size of each vector must be the same to compute covariance.");
    BOOST_ASSERT_MSG(u.size() > 0, "Computing covariance requires at least two samples.");

    Real cov = 0;
    Real mu_u = u[0];
    Real mu_v = v[0];
    Real Qu = 0;
    Real Qv = 0;

    for(std::size_t i = 1; i < u.size(); ++i)
    {
        Real u_tmp = u[i] - mu_u;
        Real v_tmp = v[i] - mu_v;
        Qu = Qu + (i*u_tmp*u_tmp)/(i+1);
        Qv = Qv + (i*v_tmp*v_tmp)/(i+1);
        cov += i*u_tmp*v_tmp/(i+1);
        mu_u = mu_u + u_tmp/(i+1);
        mu_v = mu_v + v_tmp/(i+1);
    }

    // If both datasets are constant, then they are perfectly correlated.
    if (Qu == 0 && Qv == 0)
    {
        return Real(1);
    }
    // If one dataset is constant and the other isn't, then they have no correlation:
    if (Qu == 0 || Qv == 0)
    {
        return Real(0);
    }

    // Make sure rho in [-1, 1], even in the presence of numerical noise.
    Real rho = cov/sqrt(Qu*Qv);
    if (rho > 1) {
        rho = 1;
    }
    if (rho < -1) {
        rho = -1;
    }
    return rho;
}
} // namespace detail

template<typename Container, typename Real = typename Container::value_type, typename std::enable_if<std::is_integral<Real>::value, bool>::type = true>
inline auto means_and_covariance(Container const & u, Container const & v) -> std::tuple<double, double, double>
{
    return detail::means_and_covariance_impl<std::tuple<double, double, double>>(u, v);
}

template<typename Container, typename Real = typename Container::value_type, typename std::enable_if<!std::is_integral<Real>::value, bool>::type = true>
inline auto means_and_covariance(Container const & u, Container const & v) -> std::tuple<Real, Real, Real>
{
    return detail::means_and_covariance_impl<std::tuple<Real, Real, Real>>(u, v);
}

template<typename Container, typename Real = typename Container::value_type, typename std::enable_if<std::is_integral<Real>::value, bool>::type = true>
inline double covariance(Container const & u, Container const & v)
{
    std::tuple<double, double, double> temp = detail::means_and_covariance_impl<std::tuple<double, double, double>>(u, v);
    return std::get<2>(temp);
}

template<typename Container, typename Real = typename Container::value_type, typename std::enable_if<!std::is_integral<Real>::value, bool>::type = true>
inline Real covariance(Container const & u, Container const & v)
{
    std::tuple<Real, Real, Real> temp = detail::means_and_covariance_impl<std::tuple<Real, Real, Real>>(u, v);
    return std::get<2>(temp);
}

template<typename Container, typename Real = typename Container::value_type, typename std::enable_if<std::is_integral<Real>::value, bool>::type = true>
inline double correlation_coefficient(Container const & u, Container const & v)
{
    return detail::correlation_coefficient_impl<double>(u, v);
}

template<typename Container, typename Real = typename Container::value_type, typename std::enable_if<!std::is_integral<Real>::value, bool>::type = true>
inline Real correlation_coefficient(Container const & u, Container const & v)
{
    return detail::correlation_coefficient_impl<Real>(u, v);
}

}}} // namespace boost::math::statistics
#endif
