//  (C) Copyright Nick Thompson 2019.
//  (C) Copyright Matt Borland 2021.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_STATISTICS_T_TEST_HPP
#define BOOST_MATH_STATISTICS_T_TEST_HPP

#include <cmath>
#include <iterator>
#include <utility>
#include <type_traits>
#include <boost/math/distributions/students_t.hpp>
#include <boost/math/statistics/univariate_statistics.hpp>

namespace boost { namespace math { namespace statistics { namespace detail {

template<typename ReturnType, typename T>
ReturnType one_sample_t_test_impl(T sample_mean, T sample_variance, T num_samples, T assumed_mean) 
{
    using Real = typename std::tuple_element<0, ReturnType>::type;
    using std::sqrt;
    typedef boost::math::policies::policy<
          boost::math::policies::promote_float<false>,
          boost::math::policies::promote_double<false> >
          no_promote_policy;

    Real test_statistic = (sample_mean - assumed_mean)/sqrt(sample_variance/num_samples);
    auto student = boost::math::students_t_distribution<Real, no_promote_policy>(num_samples - 1);
    Real pvalue;
    if (test_statistic > 0) {
        pvalue = 2*boost::math::cdf<Real>(student, -test_statistic);;
    }
    else {
        pvalue = 2*boost::math::cdf<Real>(student, test_statistic);
    }
    return std::make_pair(test_statistic, pvalue);
}

template<typename ReturnType, typename ForwardIterator>
ReturnType one_sample_t_test_impl(ForwardIterator begin, ForwardIterator end, typename std::iterator_traits<ForwardIterator>::value_type assumed_mean) 
{
    using Real = typename std::tuple_element<0, ReturnType>::type;
    std::pair<Real, Real> temp = mean_and_sample_variance(begin, end);
    Real mu = std::get<0>(temp);
    Real s_sq = std::get<1>(temp);
    return one_sample_t_test_impl<ReturnType>(mu, s_sq, Real(std::distance(begin, end)), Real(assumed_mean));
}
} // namespace detail

template<typename Real, typename std::enable_if<std::is_integral<Real>::value, bool>::type = true>
inline auto one_sample_t_test(Real sample_mean, Real sample_variance, Real num_samples, Real assumed_mean) -> std::pair<double, double>
{
    return detail::one_sample_t_test_impl<std::pair<double, double>>(sample_mean, sample_variance, num_samples, assumed_mean);
}

template<typename Real, typename std::enable_if<!std::is_integral<Real>::value, bool>::type = true>
inline auto one_sample_t_test(Real sample_mean, Real sample_variance, Real num_samples, Real assumed_mean) -> std::pair<Real, Real>
{
    return detail::one_sample_t_test_impl<std::pair<Real, Real>>(sample_mean, sample_variance, num_samples, assumed_mean);
}

template<typename ForwardIterator, typename Real = typename std::iterator_traits<ForwardIterator>::value_type, 
         typename std::enable_if<std::is_integral<Real>::value, bool>::type = true>
inline auto one_sample_t_test(ForwardIterator begin, ForwardIterator end, Real assumed_mean) -> std::pair<double, double>
{
    return detail::one_sample_t_test_impl<std::pair<double, double>>(begin, end, assumed_mean);
}

template<typename ForwardIterator, typename Real = typename std::iterator_traits<ForwardIterator>::value_type, 
         typename std::enable_if<!std::is_integral<Real>::value, bool>::type = true>
inline auto one_sample_t_test(ForwardIterator begin, ForwardIterator end, Real assumed_mean) -> std::pair<Real, Real>
{
    return detail::one_sample_t_test_impl<std::pair<Real, Real>>(begin, end, assumed_mean);
}

template<typename Container, typename Real = typename Container::value_type,
         typename std::enable_if<std::is_integral<Real>::value, bool>::type = true>
inline auto one_sample_t_test(Container const & v, Real assumed_mean) -> std::pair<double, double>
{
    return detail::one_sample_t_test_impl<std::pair<double, double>>(std::begin(v), std::end(v), assumed_mean);
}

template<typename Container, typename Real = typename Container::value_type,
         typename std::enable_if<!std::is_integral<Real>::value, bool>::type = true>
inline auto one_sample_t_test(Container const & v, Real assumed_mean) -> std::pair<Real, Real>
{
    return detail::one_sample_t_test_impl<std::pair<Real, Real>>(std::begin(v), std::end(v), assumed_mean);
}

}}} // namespace boost::math::statistics
#endif
