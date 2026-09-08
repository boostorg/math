//  (C) Copyright Matt Borland 2026.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
//  The segmented sieve engine: iterates segments of sieve_bytes bytes over [start, stop],
//  pre-sieves, feeds sieving primes lazily (a prime joins when the segment reaches p * p),
//  runs the size classes, masks the range ends and hands each finished segment to a sink.

#ifndef BOOST_MATH_SF_DETAIL_PRIME_SIEVE_ENGINE_HPP
#define BOOST_MATH_SF_DETAIL_PRIME_SIEVE_ENGINE_HPP

#include <boost/math/tools/config.hpp>
#include <boost/math/special_functions/detail/prime_sieve/layout.hpp>
#include <boost/math/special_functions/detail/prime_sieve/options.hpp>
#include <boost/math/special_functions/detail/prime_sieve/presieve.hpp>
#include <boost/math/special_functions/detail/prime_sieve/erat_small.hpp>
#include <boost/math/special_functions/detail/prime_sieve/erat_medium.hpp>
#include <boost/math/special_functions/detail/prime_sieve/erat_big.hpp>

#ifndef BOOST_MATH_BUILD_MODULE
#include <cstdint>
#include <cstddef>
#include <cstring>
#include <vector>
#include <algorithm>
#endif

namespace boost::math::detail::prime_sieve {

class segment_sieve
{
public:
    // sieving_primes must hold the primes in [167, sqrt(stop)] in ascending order.
    segment_sieve(const sieve_geometry& geometry, const std::vector<std::uint32_t>& sieving_primes)
        : geo_(geometry), primes_(&sieving_primes)
    {
        words_.resize(geo_.sieve_bytes / 8);
        reset_range(geo_.start, geo_.stop);
    }

    // Prepares a new sub-range with the same geometry; the range must lie within the
    // geometry's [start, stop] so that sieving primes up to sqrt(stop) suffice.
    void reset_range(std::uint64_t start, std::uint64_t stop)
    {
        start_ = start < 7 ? 7 : start;
        stop_ = stop;
        next_prime_ = 0;
        small_.init(stop_, geo_.l1_bytes);
        medium_.init(stop_);
        big_.init(stop_, geo_.sqrt_stop, geo_.sieve_bytes, geo_.has_big);
    }

    template <class Sink>
    void run(Sink& sink)
    {
        if (stop_ < start_)
        {
            return;
        }
        // The byte holding start_ covers offsets 7..36 above segment_low
        std::uint64_t segment_low {30u * ((start_ - 7) / 30)};
        std::uint8_t* sieve {reinterpret_cast<std::uint8_t*>(words_.data())};
        const std::size_t sieve_bytes {geo_.sieve_bytes};

        while (true)
        {
            const std::uint64_t remaining {stop_ - segment_low};
            if (remaining < 7)
            {
                break;
            }
            const std::uint64_t needed {(remaining - 7) / 30 + 1};
            const bool last {needed <= sieve_bytes};
            const std::size_t bytes {last ? static_cast<std::size_t>(needed) : sieve_bytes};
            const std::uint64_t span {30u * static_cast<std::uint64_t>(bytes) + 1};
            const std::uint64_t segment_high {span <= remaining ? segment_low + span : stop_};

            presieve_segment(sieve, bytes, segment_low);
            feed(segment_low, segment_high);

            small_.cross_off(sieve, bytes);
            medium_.cross_off(sieve, bytes);
            if (geo_.has_big)
            {
                big_.cross_off(sieve);
            }

            if (segment_low < start_)
            {
                sieve[0] &= keep_low[static_cast<std::size_t>(start_ - segment_low)];
            }
            if (last)
            {
                const std::uint64_t e {stop_ - segment_low - 30u * static_cast<std::uint64_t>(bytes - 1)};
                sieve[bytes - 1] &= keep_high[e > 36 ? 36 : static_cast<std::size_t>(e)];
            }

            const std::size_t n_words {(bytes + 7) / 8};
            if (bytes < n_words * 8)
            {
                std::memset(sieve + bytes, 0, n_words * 8 - bytes);
            }
            sink.segment(words_.data(), n_words, segment_low);

            if (last)
            {
                break;
            }
            segment_low += 30u * static_cast<std::uint64_t>(sieve_bytes);
        }
        sink.flush();
    }

    const sieve_geometry& geometry() const noexcept
    {
        return geo_;
    }

private:
    // Adds the sieving primes whose square lies at or below segment_high.
    void feed(std::uint64_t segment_low, std::uint64_t segment_high)
    {
        const std::vector<std::uint32_t>& primes {*primes_};
        while (next_prime_ < primes.size())
        {
            const std::uint64_t p {primes[next_prime_]};
            if (p * p > segment_high)
            {
                break;
            }
            if (p <= geo_.max_small)
            {
                small_.add(p, segment_low);
            }
            else if (p <= geo_.max_medium)
            {
                medium_.add(p, segment_low);
            }
            else
            {
                big_.add(p, segment_low);
            }
            ++next_prime_;
        }
    }

    sieve_geometry geo_;
    const std::vector<std::uint32_t>* primes_;
    std::vector<std::uint64_t> words_;
    erat_small small_;
    erat_medium medium_;
    erat_big big_;
    std::uint64_t start_ {7};
    std::uint64_t stop_ {7};
    std::size_t next_prime_ {0};
};

} // namespace boost::math::detail::prime_sieve

#endif // BOOST_MATH_SF_DETAIL_PRIME_SIEVE_ENGINE_HPP
