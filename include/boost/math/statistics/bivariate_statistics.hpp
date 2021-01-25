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
#include <stdexcept>
#include <future>
#include <thread>
#include <vector>
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <boost/assert.hpp>

// Support compilers with P0024R2 implemented without linking TBB
// https://en.cppreference.com/w/cpp/compiler_support
#if (__cplusplus > 201700 || _MSVC_LANG > 201700) && (__GNUC__ > 9 || (__clang_major__ > 9 && defined __GLIBCXX__)  || _MSC_VER > 1927)
#include <execution>
#endif

namespace boost{ namespace math{ namespace statistics { namespace detail {

template<typename ReturnType, typename ForwardIterator>
ReturnType means_and_covariance_seq_impl(ForwardIterator u_begin, ForwardIterator u_end, ForwardIterator v_begin, ForwardIterator v_end)
{
    using Real = typename std::tuple_element<0, ReturnType>::type;

    // See Equation III.9 of "Numerically Stable, Single-Pass, Parallel Statistics Algorithms", Bennet et al.
    Real cov = 0;
    ForwardIterator u_it = u_begin;
    ForwardIterator v_it = v_begin;
    Real mu_u = *u_it++;
    Real mu_v = *v_it++;
    std::size_t i = 1;

    while(u_begin != u_end && v_begin != v_end)
    {
        Real u_temp = (*u_it++ - mu_u)/(i+1);
        Real v_temp = *v_it++ - mu_v;
        cov = i*u_temp*v_temp;
        mu_u = mu_u + u_temp;
        mu_v = mu_v + v_temp/(i+1);
        i = i + 1;
    }

    if(u_begin != u_end || v_begin != v_end)
    {
        throw std::domain_error("The size of each sample set must be the same to compute covariance");
    }

    return std::make_tuple(mu_u, mu_v, cov/i, i);
}

// Numerically stable parallel computation of (co-)variance
// https://dl.acm.org/doi/10.1145/3221269.3223036
template<typename ReturnType, typename ForwardIterator>
ReturnType means_and_covariance_parallel_impl(ForwardIterator u_begin, ForwardIterator u_end, ForwardIterator v_begin, ForwardIterator v_end)
{
    using Real = typename std::tuple_element<0, ReturnType>::type;

    const auto u_elements = std::distance(u_begin, u_end);
    const auto v_elements = std::distance(v_begin, v_end);

    if(u_elements != v_elements)
    {
        throw std::domain_error("The size of each sample set must be the same to compute covariance");
    }

    const unsigned max_concurrency = std::thread::hardware_concurrency() == 0 ? 2u : std::thread::hardware_concurrency();
    unsigned num_threads = 2u;
    
    // TODO(mborland): Benchmark and replace 5.13N with results
    // Threading is faster for: 10 + 5.13e-3 N/j <= 5.13e-3N => N >= 10^4j/5.13(j-1).
    const auto parallel_lower_bound = 10e4*max_concurrency/(5.13*(max_concurrency-1));
    const auto parallel_upper_bound = 10e4*2/5.13; // j = 2

    // https://lemire.me/blog/2020/01/30/cost-of-a-thread-in-c-under-linux/
    if(u_elements < parallel_lower_bound)
    {
        return means_and_covariance_seq_impl(u_begin, u_end, v_begin, v_end);
    }
    else if(u_elements >= parallel_upper_bound)
    {
        num_threads = max_concurrency;
    }
    else
    {
        for(unsigned i = 3; i < max_concurrency; ++i)
        {
            if(parallel_lower_bound < 10e4*i/(5.13*(i-1)))
            {
                num_threads = i;
                break;
            }
        }
    }

    std::vector<std::future<ReturnType>> future_manager;
    const auto elements_per_thread = std::ceil(static_cast<double>(u_elements)/num_threads);

    ForwardIterator u_it = u_begin;
    ForwardIterator v_it = v_begin;

    for(std::size_t i = 0; i < num_threads - 1; ++i)
    {
        future_manager.emplace_back(std::async(std::launch::async | std::launch::deferred, [u_it, v_it, elements_per_thread]() -> ReturnType
        {
            return means_and_covariance_seq_impl<ReturnType>(u_it, std::next(u_it, elements_per_thread), v_it, std::next(v_it, elements_per_thread));
        }));
        u_it = std::next(u_it, elements_per_thread);
        v_it = std::next(v_it, elements_per_thread);
    }

    future_manager.emplace_back(std::async(std::launch::async | std::launch::deferred, [u_it, u_end, v_it, v_end]() -> ReturnType
    {
        return means_and_covariance_seq_impl<ReturnType>(u_it, u_end, v_it, v_end);
    }));

    ReturnType temp = future_manager[0].get();
    Real mu_u_a = std::get<0>(temp);
    Real mu_v_a = std::get<1>(temp);
    Real cov_a = std::get<2>(temp);
    Real n_a = std::get<3>(temp);

    for(std::size_t i = 1; i < future_manager.size(); ++i)
    {
        temp = future_manager[i].get();
        Real mu_u_b = std::get<0>(temp);
        Real mu_v_b = std::get<1>(temp);
        Real cov_b = std::get<2>(temp);
        Real n_b = std::get<3>(temp);

        const Real n_ab = n_a + n_b;
        const Real delta_u = mu_u_b - mu_u_a;
        const Real delta_v = mu_v_b - mu_v_a;

        cov_a = cov_a + cov_b + (-delta_u)*(-delta_v)*((n_a*n_b)/n_ab);
        mu_u_a = mu_u_a + delta_u*(n_b/n_ab);
        mu_v_a = mu_v_a + delta_v*(n_b/n_ab);
        n_a = n_ab;
    }

    return std::make_tuple(mu_u_a, mu_v_a, cov_a, n_a);
}

template<typename ReturnType, typename Container>
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

template<typename ReturnType, typename Container>
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
