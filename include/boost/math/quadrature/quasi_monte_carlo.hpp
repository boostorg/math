/*
 *  (C) Copyright Nick Thompson 2018.
 *  Use, modification and distribution are subject to the
 *  Boost Software License, Version 1.0. (See accompanying file
 *  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
 */

#ifndef BOOST_MATH_QUADRATURE_QUASI_MONTE_CARLO_HPP
#define BOOST_MATH_QUADRATURE_QUASI_MONTE_CARLO_HPP
#include <algorithm>
#include <vector>
#include <atomic>
#include <functional>
#include <future>
#include <thread>
#include <initializer_list>
#include <utility>
#include <random>
#include <chrono>
#include <map>

namespace boost { namespace math { namespace quadrature {

template<class Real, class F, class Policy = boost::math::policies::policy<>>
class quasi_monte_carlo
{
public:
    quasi_monte_carlo(const F& f,
                      std::vector<std::pair<Real, Real>> const & bounds,
                      Real error_goal,
                      size_t threads = std::thread::hardware_concurrency()): m_f{f}, m_num_threads{threads}
    {
        using std::isfinite;
        size_t n = bounds.size();
        m_lbs.resize(n);
        m_dxs.resize(n);
        m_volume = 1;
        static const char* function = "boost::math::quadrature::naive_monte_carlo<%1%>";
        for (size_t i = 0; i < n; ++i)
        {
            if (!(isfinite(bounds[i].first) && isfinite(bounds[i].second)))
            {
                boost::math::policies::raise_domain_error(function, "The routine only support bounded domains. Rescaling infinite domains must be done by the user.\n", bounds[i].first, Policy());
                return;
            }
            if (bounds[i].second <= bounds[i].first)
            {
                boost::math::policies::raise_domain_error(function, "The upper bound is <= the lower bound.\n", bounds[i].second, Policy());
                return;
            }
            m_lbs[i] = bounds[i].first;
            m_dxs[i] = bounds[i].second - m_lbs[i];
            m_volume *= m_dxs[i];
        }
        m_error_goal = error_goal;
        m_start = std::chrono::system_clock::now();
        m_done = false;
        m_total_calls = 0;
        if (m_num_threads < 2)
        {
            m_num_threads = 2;
        }
        m_variance = std::numeric_limits<Real>::max();
    }

    std::future<Real> integrate()
    {
        // Set done to false in case we wish to restart:
        m_done = false;
        return std::async(std::launch::async,
                          &naive_monte_carlo::m_integrate, this);
    }

    void cancel()
    {
        m_done = true;
    }

    void update_target_error(Real new_target_error)
    {
        m_error_goal = new_target_error;
    }

    Real current_estimate() const
    {
        return m_avg.load()*m_volume;
    }

    size_t calls() const
    {
        return m_total_calls.load();
    }

private:

    Real m_integrate()
    {
        m_start = std::chrono::system_clock::now();
        std::vector<std::thread> threads(m_num_threads);
        for (size_t i = 0; i < threads.size(); ++i)
        {
            threads[i] = std::thread(&quasi_monte_carlo::m_thread_task, this, i);
        }
        do {
            std::this_thread::sleep_for(std::chrono::milliseconds(500));
            // Allow cancellation:
            if (m_done)
            {
                break;
            }
        } while (this->current_error_estimate() > m_error_goal);
        // Error bound met; signal the threads:
        m_done = true;
        // Wait for each one to finish:
        std::for_each(threads.begin(), threads.end(),
                      std::mem_fn(&std::thread::join));
        if (m_exception)
        {
            std::rethrow_exception(m_exception);
        }
        return m_avg.load()*m_volume;
    }

    void m_thread_task(size_t thread_index)
    {
        try
        {
            std::vector<Real> x(m_lbs.size());
            std::random_device rd;
            auto seed = rd();
            std::mt19937_64 gen(seed);
            Real inv_denom = (Real) 1/( (Real) gen.max() + (Real) 2);
            Real M1 = m_thread_averages[thread_index];
            Real S = m_thread_Ss[thread_index];
            // Kahan summation is required. See the implementation discussion.
            Real compensator = 0;
            size_t k = m_thread_calls[thread_index];
            while (!m_done)
            {
                int j = 0;
                int magic_calls_before_update = 2048;
                while (j++ < magic_calls_before_update)
                {
                    for (size_t i = 0; i < m_lbs.size(); ++i)
                    {
                        x[i] = m_lbs[i] + (gen()+1)*inv_denom*m_dxs[i];
                    }
                    Real f = m_f(x);
                    ++k;
                    Real term = (f - M1)/k;
                    Real y1 = term - compensator;
                    Real M2 = M1 + y1;
                    compensator = (M2 - M1) - y1;
                    S += (f - M1)*(f - M2);
                    M1 = M2;
                }
                m_thread_averages[thread_index] = M1;
                m_thread_Ss[thread_index] = S;
                m_thread_calls[thread_index] = k;
            }
        }
        catch (...)
        {
            // Signal the other threads that the computation is ruined:
            m_done = true;
            m_exception = std::current_exception();
        }
    }

    const F& m_f;
    size_t m_num_threads;
    std::atomic<Real> m_error_goal;
    std::atomic<bool> m_done;
    std::vector<Real> m_lbs;
    std::vector<Real> m_dxs;
    Real m_volume;
    std::atomic<size_t> m_total_calls;
    // I wanted these to be vectors rather than maps,
    // but you can't resize a vector of atomics.
    std::map<size_t, std::atomic<size_t>> m_thread_calls;
    std::atomic<Real> m_variance;
    std::map<size_t, std::atomic<Real>> m_thread_Ss;
    std::atomic<Real> m_avg;
    std::map<size_t, std::atomic<Real>> m_thread_averages;
    std::chrono::time_point<std::chrono::system_clock> m_start;
    std::exception_ptr m_exception;
};
}}}
#endif
