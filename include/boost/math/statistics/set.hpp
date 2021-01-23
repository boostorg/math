//  (C) Copyright Matt Borland 2021.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_STATISTICS_SET
#define BOOST_MATH_STATISTICS_SET

#include <boost/math/statistics/univariate_statistics.hpp>
#include <tuple>
#include <iterator>
#include <utility>
#include <type_traits>
#include <list>
#include <execution>
#include <cmath>

namespace boost::math::statistics
{

enum class metrics : unsigned
{
    mean = 0b0000001,
    median = 0b0000010,
    mode = 0b0000100,
    stddev = 0b0001000,
    variance = 0b0010000,
    skewness = 0b0100000,
    first_four_moments = 0b1000000,
    all = 0b1111111
};

namespace detail
{
template<typename ExecutionPolicy, typename ForwardIterator, typename T>
class stats_base
{

protected:
    ExecutionPolicy exec_;
    ForwardIterator first_;
    ForwardIterator last_;
    T mean_ = T(0);
    T median_= T(0);
    std::list<T> mode_;
    T stddev_= T(0);
    T variance_= T(0);
    T skewness_= T(0);
    std::tuple<T, T, T, T> first_four_moments_ = std::make_tuple(T(0), T(0), T(0), T(0));

    void calc_all();

public:
    stats_base(ExecutionPolicy exec, ForwardIterator first, ForwardIterator last) noexcept : exec_ {exec}, first_ {first}, last_ {last} {}

    void calc(metrics);

    [[nodiscard]] inline T mean() const noexcept { return mean_; }
    [[nodiscard]] inline T median() const noexcept { return median_; }
    [[nodiscard]] inline std::list<T> mode() const noexcept { return mode_; }
    [[nodiscard]] inline T stddev() const noexcept { return stddev_; }
    [[nodiscard]] inline T variance() const noexcept { return variance_; }
    [[nodiscard]] inline T skewness() const noexcept { return skewness_; }
    [[nodiscard]] inline std::tuple<T, T, T, T> first_four_moments() const noexcept { return first_four_moments_; }
};

template<typename ExecutionPolicy, typename ForwardIterator, typename T>
void stats_base<ExecutionPolicy, ForwardIterator, T>::calc_all()
{
    using std::sqrt;
    
    if(mean_ == T(0))
    {
        mean_ = boost::math::statistics::mean(exec_, first_, last_);
    }

    if(median_ == T(0))
    {
        median_ = boost::math::statistics::median(exec_, first_, last_);
    }
    
    if(mode_.empty())
    {
        boost::math::statistics::mode(exec_, first_, last_, std::inserter(mode_, mode_.begin()));
    }
    
    if(variance_ == T(0))
    {
        variance_ = boost::math::statistics::variance(exec_, first_, last_);
    }
    
    if(stddev_ == T(0))
    {
        stddev_ = sqrt(variance_);
    }
    
    if(skewness_ == T(0))
    {
        skewness_ = boost::math::statistics::skewness(exec_, first_, last_);
    }

    if(first_four_moments_ == std::make_tuple(T(0), T(0), T(0), T(0)))
    {
        first_four_moments_ = boost::math::statistics::first_four_moments(exec_, first_, last_);
    }
}

template<typename ExecutionPolicy, typename ForwardIterator, typename T>
void stats_base<ExecutionPolicy, ForwardIterator, T>::calc(metrics metric)
{
    using std::sqrt;
    
    switch(metric)
    {
        case metrics::mean:
        {
            mean_ = boost::math::statistics::mean(exec_, first_, last_);
            break;
        }
        case metrics::median:
        {
            median_ = boost::math::statistics::median(exec_, first_, last_);
            break;
        }
        case metrics::mode:
        {
            boost::math::statistics::mode(exec_, first_, last_, std::inserter(mode_, mode_.begin()));
            break;
        }
        case metrics::stddev:
        {
            T temp = boost::math::statistics::variance(exec_, first_, last_);
            stddev_ = sqrt(temp);
            break;
        }
        case metrics::variance:
        {
            variance_ = boost::math::statistics::variance(exec_, first_, last_);
            break;
        }
        case metrics::skewness:
        {
            skewness_ = boost::math::statistics::skewness(exec_, first_, last_);
            break;
        }
        case metrics::first_four_moments:
        {
            first_four_moments_ = boost::math::statistics::first_four_moments(exec_, first_, last_);
            break;
        }
        case metrics::all:
            [[fallthrough]];
        default:
            calc_all();
    }
}

} // End namespace detail

// Allow multiple interfaces to base allowing an execution policy to be passed or not
template<typename ExecutionPolicy = decltype(std::execution::seq), typename ForwardIterator = void, 
         typename T = typename std::iterator_traits<ForwardIterator>::value_type>
class stats : public detail::stats_base<ExecutionPolicy, ForwardIterator, T>
{
public:
    stats(ExecutionPolicy exec, ForwardIterator first, ForwardIterator last) : 
        detail::stats_base<ExecutionPolicy, ForwardIterator, T>(exec, first, last) {}
       
    stats(ForwardIterator first, ForwardIterator last) : 
        detail::stats_base<decltype(std::execution::seq), ForwardIterator, T>(std::execution::seq, first, last) {}
};

} // namespace boost::math::statistics

#endif // BOOST_MATH_STATISTICS_SET
