//  (C) Copyright Matt Borland 2026.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
//  Multithreaded driver: the range is cut into chunks of about 200 * sqrt(stop) that are
//  handed out dynamically through an atomic counter; each worker owns a segment_sieve and
//  shares the read-only sieving primes. Output chunks are committed in order through a
//  bounded queue so the caller sees ascending primes without padding.

#ifndef BOOST_MATH_SF_DETAIL_PRIME_SIEVE_PARALLEL_HPP
#define BOOST_MATH_SF_DETAIL_PRIME_SIEVE_PARALLEL_HPP

#include <boost/math/tools/config.hpp>
#include <boost/math/special_functions/detail/prime_sieve/execution.hpp>

#ifdef BOOST_MATH_PRIME_SIEVE_HAS_THREADS

#include <boost/math/tools/bit.hpp>
#include <boost/math/special_functions/detail/prime_sieve/driver.hpp>

#ifndef BOOST_MATH_BUILD_MODULE
#include <cstdint>
#include <cstddef>
#include <cmath>
#include <vector>
#include <thread>
#include <future>
#include <atomic>
#include <mutex>
#include <condition_variable>
#include <exception>
#include <algorithm>
#include <utility>
#endif

namespace boost::math::detail::prime_sieve {

inline constexpr std::uint64_t min_chunk_width {10000000};

inline std::vector<std::uint32_t> sieving_primes_upto_parallel(std::uint64_t n, const prime_sieve_options& options);

struct chunk_plan
{
    std::uint64_t start {0};
    std::uint64_t stop {0};
    std::uint64_t chunk {0};   // width in integers, multiple of 30
    std::uint64_t iters {1};
    unsigned threads {1};
};

inline unsigned worker_count(const prime_sieve_options& options) noexcept
{
    if (options.max_threads != 0)
    {
        return options.max_threads;
    }
    const unsigned hw {std::thread::hardware_concurrency()};
    return hw == 0 ? 2u : hw;
}

// primesieve's balancing: chunks of at most 200 * sqrt(stop) so the O(pi(sqrt stop)) setup
// per chunk stays below one percent, an iteration count that is a multiple of the thread
// count, and a floor of 1e7 integers per chunk.
inline chunk_plan plan_chunks(std::uint64_t start, std::uint64_t stop, const prime_sieve_options& options, bool extracting) noexcept
{
    chunk_plan plan {};
    plan.start = start;
    plan.stop = stop;
    const std::uint64_t dist {stop - start + 1};
    const unsigned hw {worker_count(options)};
    const std::uint64_t root {boost::math::tools::isqrt(stop)};
    const std::uint64_t threshold {(std::max)(root / 5, min_chunk_width)};
    std::uint64_t threads {dist / threshold};
    threads = (std::max)(threads, std::uint64_t(1));
    threads = (std::min)(threads, static_cast<std::uint64_t>(hw));
    if (threads <= 1)
    {
        plan.chunk = dist;
        plan.iters = 1;
        plan.threads = 1;
        return plan;
    }

    std::uint64_t chunk {(std::min)(200 * root, dist / threads)};
    if (extracting)
    {
        // Keep one chunk's worth of primes near chunk_primes entries
        const double per_chunk {static_cast<double>(options.chunk_primes) * std::log(static_cast<double>(stop))};
        if (per_chunk < static_cast<double>(chunk))
        {
            chunk = static_cast<std::uint64_t>(per_chunk);
        }
    }
    chunk = (std::max)(chunk, min_chunk_width);
    std::uint64_t iters {dist / chunk};
    iters = (iters / threads) * threads;
    iters = (std::max)(iters, threads);
    chunk = (dist - 1) / iters + 1;
    chunk = (std::max)(chunk, min_chunk_width);
    chunk += 30 - chunk % 30;
    iters = (dist - 1) / chunk + 1;

    plan.chunk = chunk;
    plan.iters = iters;
    plan.threads = static_cast<unsigned>((std::min)(threads, iters));
    return plan;
}

// Bounds of chunk i, inclusive.
inline void chunk_bounds(const chunk_plan& plan, std::uint64_t i, std::uint64_t& lo, std::uint64_t& hi) noexcept
{
    lo = plan.start + i * plan.chunk;
    const std::uint64_t remaining {plan.stop - lo};
    hi = remaining < plan.chunk ? plan.stop : lo + plan.chunk - 1;
}

// Counts the primes in [start, stop] (start >= 7) on several threads.
inline std::uint64_t parallel_count(std::uint64_t start, std::uint64_t stop, const prime_sieve_options& options)
{
    if (stop < start)
    {
        return 0;
    }
    const chunk_plan plan {plan_chunks(start, stop, options, false)};
    const sieve_geometry g {make_geometry(start, stop, options)};
    const std::vector<std::uint32_t> primes {sieving_primes_upto_parallel(g.sqrt_stop, options)};
    if (plan.threads <= 1)
    {
        // Too narrow to split, but the sieving primes above were still generated in parallel
        segment_sieve engine {g, primes};
        count_sink sink {};
        engine.run(sink);
        return sink.count;
    }

    alignas(64) std::atomic<std::uint64_t> next {0};

    auto worker = [&]() -> std::uint64_t
    {
        segment_sieve engine {g, primes};
        count_sink sink {};
        std::uint64_t i {};
        while ((i = next.fetch_add(1, std::memory_order_relaxed)) < plan.iters)
        {
            std::uint64_t lo {};
            std::uint64_t hi {};
            chunk_bounds(plan, i, lo, hi);
            engine.reset_range(lo, hi);
            engine.run(sink);
        }
        return sink.count;
    };

    std::vector<std::future<std::uint64_t>> futures;
    futures.reserve(plan.threads);
    for (unsigned t {0}; t < plan.threads; ++t)
    {
        try
        {
            futures.push_back(std::async(std::launch::async, worker));
        }
        catch (const std::system_error&)
        {
            break;
        }
    }
    // Any chunks left over (thread creation failed) are processed here
    std::uint64_t total {worker()};
    for (auto& f : futures)
    {
        total += f.get();
    }
    return total;
}

// In-order hand-off of per-chunk prime buffers from workers to the consuming thread.
class ordered_chunk_queue
{
public:
    explicit ordered_chunk_queue(std::size_t capacity) : slots_(capacity), ready_(capacity, false)
    {
    }

    // Returns the empty buffer for chunk i once the chunk that previously used the slot is drained.
    std::vector<std::uint64_t>& acquire(std::uint64_t i)
    {
        std::unique_lock<std::mutex> lock {mutex_};
        free_.wait(lock, [&] { return failed_ || i < next_ + slots_.size(); });
        return slots_[static_cast<std::size_t>(i % slots_.size())];
    }

    void publish(std::uint64_t i)
    {
        {
            std::lock_guard<std::mutex> lock {mutex_};
            ready_[static_cast<std::size_t>(i % slots_.size())] = true;
        }
        ready_cv_.notify_all();
    }

    void fail(std::exception_ptr e)
    {
        {
            std::lock_guard<std::mutex> lock {mutex_};
            if (!failed_)
            {
                failed_ = true;
                error_ = e;
            }
        }
        ready_cv_.notify_all();
        free_.notify_all();
    }

    bool failed() const
    {
        std::lock_guard<std::mutex> lock {mutex_};
        return failed_;
    }

    // Blocks until chunk i is published; returns nullptr when a worker failed.
    std::vector<std::uint64_t>* wait_ready(std::uint64_t i)
    {
        std::unique_lock<std::mutex> lock {mutex_};
        ready_cv_.wait(lock, [&] { return failed_ || ready_[static_cast<std::size_t>(i % slots_.size())]; });
        if (failed_)
        {
            return nullptr;
        }
        return &slots_[static_cast<std::size_t>(i % slots_.size())];
    }

    void release(std::uint64_t i)
    {
        {
            std::lock_guard<std::mutex> lock {mutex_};
            const std::size_t s {static_cast<std::size_t>(i % slots_.size())};
            ready_[s] = false;
            slots_[s].clear();
            next_ = i + 1;
        }
        free_.notify_all();
    }

    void rethrow()
    {
        std::exception_ptr e {};
        {
            std::lock_guard<std::mutex> lock {mutex_};
            e = error_;
        }
        if (e)
        {
            std::rethrow_exception(e);
        }
    }

private:
    mutable std::mutex mutex_;
    std::condition_variable ready_cv_;
    std::condition_variable free_;
    std::vector<std::vector<std::uint64_t>> slots_;
    std::vector<bool> ready_;
    std::uint64_t next_ {0};
    bool failed_ {false};
    std::exception_ptr error_ {};
};

// Sieves [start, stop] (start >= 7) on several threads and passes the primes in ascending
// order to consume(const std::uint64_t*, std::size_t) on the calling thread.
template <class Consumer>
void parallel_range(std::uint64_t start, std::uint64_t stop, const prime_sieve_options& options, Consumer& consume)
{
    if (stop < start)
    {
        return;
    }
    const chunk_plan plan {plan_chunks(start, stop, options, true)};
    const sieve_geometry g {make_geometry(start, stop, options)};
    const std::vector<std::uint32_t> primes {sieving_primes_upto_parallel(g.sqrt_stop, options)};
    if (plan.threads <= 1)
    {
        segment_sieve engine {g, primes};
        extract_sink<Consumer> sink {consume};
        engine.run(sink);
        return;
    }

    alignas(64) std::atomic<std::uint64_t> next {0};
    ordered_chunk_queue queue {2 * static_cast<std::size_t>(plan.threads)};

    auto worker = [&]()
    {
        try
        {
            segment_sieve engine {g, primes};
            std::uint64_t i {};
            while ((i = next.fetch_add(1, std::memory_order_relaxed)) < plan.iters)
            {
                std::uint64_t lo {};
                std::uint64_t hi {};
                chunk_bounds(plan, i, lo, hi);
                std::vector<std::uint64_t>& buffer {queue.acquire(i)};
                if (queue.failed())
                {
                    return;
                }
                buffer.reserve(static_cast<std::size_t>(prime_count_upper_bound(lo, hi)));
                append_u64 appender {buffer};
                extract_sink<append_u64> sink {appender};
                engine.reset_range(lo, hi);
                engine.run(sink);
                queue.publish(i);
            }
        }
        catch (...)
        {
            queue.fail(std::current_exception());
        }
    };

    std::vector<std::future<void>> futures;
    futures.reserve(plan.threads);
    for (unsigned t {0}; t < plan.threads; ++t)
    {
        try
        {
            futures.push_back(std::async(std::launch::async, worker));
        }
        catch (const std::system_error&)
        {
            break;
        }
    }
    if (futures.empty())
    {
        // No worker could be started: fall back to the sequential path
        extract_sink<Consumer> sink {consume};
        run_u64(start, stop, options, sink);
        return;
    }

    try
    {
        for (std::uint64_t i {0}; i < plan.iters; ++i)
        {
            std::vector<std::uint64_t>* buffer {queue.wait_ready(i)};
            if (buffer == nullptr)
            {
                break;
            }
            consume(buffer->data(), buffer->size());
            queue.release(i);
        }
    }
    catch (...)
    {
        // The consumer threw: release the workers blocked on the full queue before waiting
        // on the futures, whose destructors would otherwise block forever.
        queue.fail(std::current_exception());
        for (auto& f : futures)
        {
            f.wait();
        }
        throw;
    }
    for (auto& f : futures)
    {
        f.wait();
    }
    queue.rethrow();
}

// Sieving primes in [167, n] generated on several threads when there are many of them.
inline std::vector<std::uint32_t> sieving_primes_upto_parallel(std::uint64_t n, const prime_sieve_options& options)
{
    if (n < 100000000)
    {
        return sieving_primes_upto(n);
    }
    std::vector<std::uint32_t> out;
    out.reserve(static_cast<std::size_t>(prime_count_upper_bound(n)));
    append_u32 appender {out};
    parallel_range(167, n, options, appender);
    return out;
}

} // namespace boost::math::detail::prime_sieve

#endif // BOOST_MATH_PRIME_SIEVE_HAS_THREADS
#endif // BOOST_MATH_SF_DETAIL_PRIME_SIEVE_PARALLEL_HPP
