//  (C) Copyright Nick Thompson 2019.
//  (C) Copyright Matt Borland 2021.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_STATISTICS_T_TEST_HPP
#define BOOST_MATH_STATISTICS_T_TEST_HPP

#include <cmath>
#include <cstddef>
#include <iterator>
#include <utility>
#include <type_traits>
#include <vector>
#include <stdexcept>
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

template<typename ReturnType, typename T, typename U = std::size_t>
ReturnType two_sample_t_test_impl(T mean_1, T variance_1, U size_1, T mean_2, T variance_2, U size_2)
{
    using Real = typename std::tuple_element<0, ReturnType>::type;
    using no_promote_policy = boost::math::policies::policy<boost::math::policies::promote_float<false>, boost::math::policies::promote_double<false>>;
    using std::sqrt;

    // TODO(mborland): Check if variance_1 > 2 * variance_2 or vice versa
    // if so send to welch's t-test

    Real dof = size_1 + size_2 - 2;
    Real pooled_std_dev = sqrt(((size_1-1)*variance_1 + (size_2-1)*variance_2) / dof);
    Real test_statistic = (mean_1-mean_2) / (pooled_std_dev*sqrt(1.0/static_cast<Real>(size_1) + 1.0/static_cast<Real>(size_2)));

    auto student = boost::math::students_t_distribution<Real, no_promote_policy>(dof);
    Real pvalue;
    if (test_statistic > 0) 
    {
        pvalue = 2*boost::math::cdf<Real>(student, -test_statistic);;
    }
    else 
    {
        pvalue = 2*boost::math::cdf<Real>(student, test_statistic);
    }

    return std::make_pair(test_statistic, pvalue);
}

template<typename ReturnType, typename ForwardIterator>
ReturnType two_sample_t_test_impl(ForwardIterator begin_1, ForwardIterator end_1, ForwardIterator begin_2, ForwardIterator end_2)
{
    using Real = typename std::tuple_element<0, ReturnType>::type;
    auto n1 = std::distance(begin_1, end_1);
    auto n2 = std::distance(begin_2, end_2);

    ReturnType temp_1 = mean_and_sample_variance(begin_1, end_1);
    Real mean_1 = std::get<0>(temp_1);
    Real variance_1 = std::get<1>(temp_1);

    ReturnType temp_2 = mean_and_sample_variance(begin_2, end_2);
    Real mean_2 = std::get<0>(temp_2);
    Real variance_2 = std::get<1>(temp_2);

    return two_sample_t_test_impl<ReturnType>(mean_1, variance_1, n1, mean_2, variance_2, n2);
}

template<typename ReturnType, typename ForwardIterator>
ReturnType paired_samples_t_test_impl(ForwardIterator begin_1, ForwardIterator end_1, ForwardIterator begin_2, ForwardIterator end_2)
{
    using Real = typename std::tuple_element<0, ReturnType>::type;
    using no_promote_policy = boost::math::policies::policy<boost::math::policies::promote_float<false>, boost::math::policies::promote_double<false>>;
    using std::sqrt;

    auto n1 = std::distance(begin_1, end_1);
    auto n2 = std::distance(begin_2, end_2);

    if(n1 != n2)
    {
        throw std::domain_error("Both sets must have the same number of values.");
    }
    
    std::vector<Real> delta(n1);
    ForwardIterator it_1 = begin_1;
    ForwardIterator it_2 = begin_2;
    std::size_t pos = 0;
    while(it_1 != end_1)
    {
        delta[pos++] = static_cast<Real>(*it_1++) - static_cast<Real>(*it_2++);
    }

    std::pair<Real, Real> temp = mean_and_sample_variance(delta.begin(), delta.end());
    Real delta_mean = std::get<0>(temp);
    Real delta_std_dev = sqrt(std::get<1>(temp));

    Real test_statistic = delta_mean/(delta_std_dev/sqrt(n1));

    auto student = boost::math::students_t_distribution<Real, no_promote_policy>(n1 - 1);
    Real pvalue;
    if (test_statistic > 0) 
    {
        pvalue = 2*boost::math::cdf<Real>(student, -test_statistic);;
    }
    else 
    {
        pvalue = 2*boost::math::cdf<Real>(student, test_statistic);
    }

    return std::make_pair(test_statistic, pvalue);
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

template<typename ForwardIterator, typename Real = typename std::iterator_traits<ForwardIterator>::value_type, 
         typename std::enable_if<std::is_integral<Real>::value, bool>::type = true>
inline auto two_sample_t_test(ForwardIterator begin_1, ForwardIterator end_1, ForwardIterator begin_2, ForwardIterator end_2) -> std::pair<double, double>
{
    return detail::two_sample_t_test_impl<std::pair<double, double>>(begin_1, end_1, begin_2, end_2);
}

template<typename ForwardIterator, typename Real = typename std::iterator_traits<ForwardIterator>::value_type, 
         typename std::enable_if<!std::is_integral<Real>::value, bool>::type = true>
inline auto two_sample_t_test(ForwardIterator begin_1, ForwardIterator end_1, ForwardIterator begin_2, ForwardIterator end_2) -> std::pair<Real, Real>
{
    return detail::two_sample_t_test_impl<std::pair<Real, Real>>(begin_1, end_1, begin_2, end_2);
}

template<typename Container, typename Real = typename Container::value_type, typename std::enable_if<std::is_integral<Real>::value, bool>::type = true>
inline auto two_sample_t_test(Container const & u, Container const & v) -> std::pair<double, double>
{
    return detail::two_sample_t_test_impl<std::pair<double, double>>(std::begin(u), std::end(u), std::begin(v), std::end(v));
}

template<typename Container, typename Real = typename Container::value_type, typename std::enable_if<!std::is_integral<Real>::value, bool>::type = true>
inline auto two_sample_t_test(Container const & u, Container const & v) -> std::pair<Real, Real>
{
    return detail::two_sample_t_test_impl<std::pair<Real, Real>>(std::begin(u), std::end(u), std::begin(v), std::end(v));
}

template<typename ForwardIterator, typename Real = typename std::iterator_traits<ForwardIterator>::value_type, 
         typename std::enable_if<std::is_integral<Real>::value, bool>::type = true>
inline auto paired_samples_t_test(ForwardIterator begin_1, ForwardIterator end_1, ForwardIterator begin_2, ForwardIterator end_2) -> std::pair<double, double>
{
    return detail::paired_samples_t_test_impl<std::pair<double, double>>(begin_1, end_1, begin_2, end_2);
}

template<typename ForwardIterator, typename Real = typename std::iterator_traits<ForwardIterator>::value_type, 
         typename std::enable_if<!std::is_integral<Real>::value, bool>::type = true>
inline auto paired_samples_t_test(ForwardIterator begin_1, ForwardIterator end_1, ForwardIterator begin_2, ForwardIterator end_2) -> std::pair<Real, Real>
{
    return detail::paired_samples_t_test_impl<std::pair<Real, Real>>(begin_1, end_1, begin_2, end_2);
}

template<typename Container, typename Real = typename Container::value_type, typename std::enable_if<std::is_integral<Real>::value, bool>::type = true>
inline auto paired_samples_t_test(Container const & u, Container const & v) -> std::pair<double, double>
{
    return detail::paired_samples_t_test_impl<std::pair<double, double>>(std::begin(u), std::end(u), std::begin(v), std::end(v));
}

template<typename Container, typename Real = typename Container::value_type, typename std::enable_if<!std::is_integral<Real>::value, bool>::type = true>
inline auto paired_samples_t_test(Container const & u, Container const & v) -> std::pair<Real, Real>
{
    return detail::paired_samples_t_test_impl<std::pair<Real, Real>>(std::begin(u), std::end(u), std::begin(v), std::end(v));
}

}}} // namespace boost::math::statistics
#endif
