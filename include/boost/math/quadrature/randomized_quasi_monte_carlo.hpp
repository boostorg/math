/*
 *  (C) Copyright Nick Thompson 2018.
 *  Use, modification and distribution are subject to the
 *  Boost Software License, Version 1.0. (See accompanying file
 *  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
 */

#ifndef BOOST_MATH_QUADRATURE_RANDOMIZED_QUASI_MONTE_CARLO_HPP
#define BOOST_MATH_QUADRATURE_RANDOMIZED_QUASI_MONTE_CARLO_HPP
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
#include <boost/math/tools/halton_sequence.hpp>
#include <boost/math/tools/atanassov_sequence.hpp>
#include <boost/math/tools/cranley_patterson_rotation.hpp>

namespace boost { namespace math { namespace quadrature {

  namespace detail {
     enum class limit_classification {FINITE,
                                      LOWER_BOUND_INFINITE,
                                      UPPER_BOUND_INFINITE,
                                      DOUBLE_INFINITE};
   }

template<class Real, class F, class Policy = boost::math::policies::policy<>>
class randomized_quasi_monte_carlo
{
public:
    randomized_quasi_monte_carlo(const F& integrand,
                       std::vector<std::pair<Real, Real>> const & bounds,
                       Real error_goal,
                       size_t threads = std::thread::hardware_concurrency()): m_num_threads{threads}
     {
         using std::numeric_limits;
         // Each thread produces an "observation", and these are averaged to produce the error estimate.
         // So we need at least 2 threads:
         m_num_threads = std::max(m_num_threads, static_cast<size_t>(2));
         size_t n = bounds.size();
         m_lbs.resize(n);
         m_dxs.resize(n);
         m_limit_types.resize(n);
         m_volume = 1;
         static const char* function = "boost::math::quadrature::naive_monte_carlo<%1%>";
         for (size_t i = 0; i < n; ++i)
         {
             if (bounds[i].second <= bounds[i].first)
             {
                 boost::math::policies::raise_domain_error(function, "The upper bound is <= the lower bound.\n", bounds[i].second, Policy());
                 return;
             }
             if (bounds[i].first == -numeric_limits<Real>::infinity())
             {
                 if (bounds[i].second == numeric_limits<Real>::infinity())
                 {
                     m_limit_types[i] = detail::limit_classification::DOUBLE_INFINITE;
                 }
                 else
                 {
                     m_limit_types[i] = detail::limit_classification::LOWER_BOUND_INFINITE;
                     // Ok ok this is bad:
                     m_lbs[i] = bounds[i].second;
                     m_dxs[i] = numeric_limits<Real>::quiet_NaN();
                 }
             }
             else if (bounds[i].second == numeric_limits<Real>::infinity())
             {
                 m_limit_types[i] = detail::limit_classification::UPPER_BOUND_INFINITE;
                 m_lbs[i] = bounds[i].first;
                 m_dxs[i] = numeric_limits<Real>::quiet_NaN();
             }
             else
             {
                 m_limit_types[i] = detail::limit_classification::FINITE;
                 m_lbs[i] = bounds[i].first;
                 m_dxs[i] = bounds[i].second - m_lbs[i];
                 m_volume *= m_dxs[i];
             }
         }

         m_integrand = [this, &integrand](std::vector<Real> & x)->Real
         {
             Real coeff = m_volume;
             for (size_t i = 0; i < x.size(); ++i)
             {
                 // Variable transformation are listed at:
                 // https://en.wikipedia.org/wiki/Numerical_integration
                 if (m_limit_types[i] == detail::limit_classification::FINITE)
                 {
                     x[i] = m_lbs[i] + x[i]*m_dxs[i];
                 }
                 else if (m_limit_types[i] == detail::limit_classification::UPPER_BOUND_INFINITE)
                 {
                     Real t = x[i];
                     Real z = 1/(1-t);
                     coeff *= (z*z);
                     x[i] = m_lbs[i] + t*z;
                 }
                 else if (m_limit_types[i] == detail::limit_classification::LOWER_BOUND_INFINITE)
                 {
                     Real t = x[i];
                     Real z = 1/t;
                     coeff *= (z*z);
                     x[i] = m_lbs[i] + (t-1)*z;
                 }
                 else
                 {
                     Real t = 2*x[i] - 1;
                     Real tsq = t*t;
                     Real z = 1/(1-t);
                     z /= (1+t);
                     coeff *= 2*(1+tsq)*z*z;
                     x[i] = t*z;
                 }
             }
             return coeff*integrand(x);
         };

         m_avg = 0;
         m_error_goal = error_goal;
         m_start = std::chrono::system_clock::now();
         m_done = false;
         m_total_calls = 0;
         m_variance = numeric_limits<Real>::max();
    }

    std::future<Real> integrate()
    {
        // Set done to false in case we wish to restart:
        m_done = false;
        return std::async(std::launch::async,
                          &randomized_quasi_monte_carlo::m_integrate, this);
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
        return m_avg.load();
    }

    Real current_error_estimate() const
    {
        using std::sqrt;
        return sqrt(m_variance.load());
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
            threads[i] = std::thread(&randomized_quasi_monte_carlo::m_thread_task, this, i);
        }
        do {
            std::this_thread::sleep_for(std::chrono::milliseconds(500));
            size_t total_calls = 0;
            for (size_t i = 0; i < m_num_threads; ++i)
            {
                total_calls += m_thread_calls[i];
            }
            Real variance = 0;
            Real avg = 0;
            for (size_t i = 0; i < m_num_threads; ++i)
            {
                avg += m_thread_averages[i];
            }
            avg /= m_num_threads;
            for (size_t i = 0; i < m_num_threads; ++i)
            {
                variance += (m_thread_averages[i] - avg)*(m_thread_averages[i] - avg);
            }
            variance /= (m_num_threads - 1);
            m_avg = avg;
            m_variance = variance;
            m_total_calls = total_calls;
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
        size_t total_calls = 0;
        for (size_t i = 0; i < m_num_threads; ++i)
        {
           total_calls += m_thread_calls[i];
        }
        Real variance = 0;
        Real avg = 0;
        for (size_t i = 0; i < m_num_threads; ++i)
        {
           avg += m_thread_averages[i];
        }
        avg /= m_num_threads;
        m_avg = avg;

        for (size_t i = 0; i < m_num_threads; ++i)
        {
            variance += (m_thread_averages[i] - avg)*(m_thread_averages[i] - avg);
        }
        variance /= (m_num_threads -1);
        m_variance = variance;
        m_total_calls = total_calls;

        return m_avg.load();
    }

    void m_thread_task(size_t thread_index)
    {
        try
        {
            std::vector<Real> x(m_lbs.size());
            // Ok this must be made generic. But the leaped halton sequence is a reasonable sequence!
            atanassov_sequence<Real> atanassov(m_lbs.size(), true);
            // Randomize using Cranley-Patterson:
            cranley_patterson_rotation<Real> cranley(m_lbs.size());

            Real M1 = m_thread_averages[thread_index];
            size_t k = m_thread_calls[thread_index];
            // Kahan summation is required. See the implementation discussion.
            Real compensator = 0;
            do
            {
                int j = 0;
                int magic_calls_before_update = 2048;
                while (j++ < magic_calls_before_update)
                {
                    // Fill up with quasi-random numbers:
                    atanassov(x.begin(), x.end());
                    // Apply Cranley-Patterson rotation:
                    cranley(x.begin(), x.end());
                    Real f = m_integrand(x);
                    ++k;
                    Real term = (f - M1)/k;
                    Real y1 = term - compensator;
                    Real M2 = M1 + y1;
                    compensator = (M2 - M1) - y1;
                    M1 = M2;
                }
                m_thread_averages[thread_index] = M1;
                m_thread_calls[thread_index] = k;
            } while (!m_done);
        }
        catch (...)
        {
            // Signal the other threads that the computation is ruined:
            m_done = true;
            m_exception = std::current_exception();
        }
    }

    std::function<Real(std::vector<Real> &)> m_integrand;
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
    std::atomic<Real> m_avg;
    std::map<size_t, std::atomic<Real>> m_thread_averages;
    std::chrono::time_point<std::chrono::system_clock> m_start;
    std::exception_ptr m_exception;
    std::vector<detail::limit_classification> m_limit_types;
};
}}}
#endif
