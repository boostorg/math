// Copyright 2020 John Maddock and Matt Borland
//
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_TOOLS_STOPWATCH_HPP
#define BOOST_MATH_TOOLS_STOPWATCH_HPP

#include <chrono>
#include <vector>
#include <climits>

namespace boost::math::tools
{
template<typename T = double, typename Clock = std::chrono::high_resolution_clock>
class stopwatch
{
private:
    using Time = typename Clock::time_point;
    using duration = typename Clock::duration;

    std::vector<Time> laps_;
    Time m_start_;
    Time m_stop_;

public:
    stopwatch()
    {
        reset();
    }

    inline T elapsed() noexcept
    {
        return std::chrono::duration_cast<std::chrono::duration<T>>(Clock::now() - m_start_).count();
    }

    inline void reset()
    {
        m_start_ = Clock::now();
        m_stop_ = m_start_;
        laps_.clear();
        laps_.push_back(m_start_);
    }

    inline T lap()
    {
        const auto temp {Clock::now()};
        const auto delta {temp - laps_.back()};
        laps_.push_back(temp);
        return std::chrono::duration_cast<std::chrono::duration<T>>(delta).count();
    }

    inline T stop()
    {
        m_stop_ = Clock::now();
        laps_.emplace_back(m_stop_);
        return std::chrono::duration_cast<std::chrono::duration<T>>(m_stop_ - m_start_).count();
    }

    inline std::vector<T> laps() noexcept
    {
        return laps_;
    }
};
}

#endif //BOOST_MATH_TOOLS_STOPWATCH_HPP
