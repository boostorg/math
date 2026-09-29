//  (C) Copyright Matt Borland 2026.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
//  Big sieving primes (at most one multiple per segment): Oliveira e Silva's bucket sieve.
//  Each prime is filed in the list of the segment holding its next multiple, so a segment
//  only touches the primes that actually hit it. Lists are chains of fixed 8 KiB buckets
//  from a memory pool; a write pointer sitting on an 8 KiB boundary means "bucket full".
//  Uses the mod 210 wheel.

#ifndef BOOST_MATH_SF_DETAIL_PRIME_SIEVE_ERAT_BIG_HPP
#define BOOST_MATH_SF_DETAIL_PRIME_SIEVE_ERAT_BIG_HPP

#include <boost/math/tools/config.hpp>
#include <boost/math/tools/bit.hpp>
#include <boost/math/special_functions/detail/prime_sieve/layout.hpp>

#ifndef BOOST_MATH_BUILD_MODULE
#include <cstdint>
#include <cstddef>
#include <cmath>
#include <memory>
#include <vector>
#include <utility>
#endif

namespace boost::math::detail::prime_sieve {

inline constexpr std::size_t bucket_bytes {8192};

// Entries looked ahead in a bucket; reading past a bucket's end stays inside the 8 KiB block
// (or the following one) and is harmless.
inline constexpr std::size_t prefetch_distance {16};

inline void prefetch_byte(const std::uint8_t* p) noexcept
{
#if defined(__GNUC__) || defined(__clang__)
    __builtin_prefetch(p, 1, 0);
#else
    (void)p;
#endif
}

// A bucket is exactly bucket_bytes long and aligned to bucket_bytes, so the bucket owning
// a write pointer is found by masking, and a pointer that reached the next boundary is full.
struct alignas(bucket_bytes) bucket
{
    static constexpr std::size_t capacity {(bucket_bytes - 2 * sizeof(void*)) / sizeof(sieving_prime)};

    sieving_prime* end;      // one past the last used slot, valid while the bucket is processed
    bucket* previous;        // older bucket of the same list
    sieving_prime primes[capacity];

    static bool is_full(const sieving_prime* write) noexcept
    {
        // also true for nullptr, which makes an empty list allocate its first bucket
        return (reinterpret_cast<std::uintptr_t>(write) & (bucket_bytes - 1)) == 0;
    }

    static bucket* owner(const sieving_prime* write) noexcept
    {
        return reinterpret_cast<bucket*>((reinterpret_cast<std::uintptr_t>(write) - 1) & ~static_cast<std::uintptr_t>(bucket_bytes - 1));
    }
};

static_assert(sizeof(bucket) == bucket_bytes, "bucket must be exactly one block");

// Hands out aligned buckets from large blocks; freed buckets go to a free list.
class bucket_pool
{
public:
    bucket* get()
    {
        if (free_ == nullptr)
        {
            grow();
        }
        bucket* b {free_};
        free_ = b->previous;
        b->previous = nullptr;
        return b;
    }

    void put(bucket* b) noexcept
    {
        b->previous = free_;
        free_ = b;
    }

    // Returns every bucket to the free list, keeping the blocks
    void release_all() noexcept
    {
        free_ = nullptr;
        for (const block& blk : blocks_)
        {
            for (std::size_t i {0}; i < blk.count; ++i)
            {
                put(blk.first + i);
            }
        }
    }

private:
    struct block
    {
        std::unique_ptr<std::uint8_t[]> storage;
        bucket* first;
        std::size_t count;
    };

    void grow()
    {
        // 512 KiB, then doubling up to 16 MiB per block
        const std::size_t count {blocks_.empty() ? std::size_t(64) : (std::min)(blocks_.back().count * 2, std::size_t(2048))};
        block blk {};
        blk.storage.reset(new std::uint8_t[(count + 2) * bucket_bytes]);   // one bucket of padding for alignment, one for the look-ahead reads
        const std::uintptr_t raw {reinterpret_cast<std::uintptr_t>(blk.storage.get())};
        const std::uintptr_t aligned {(raw + bucket_bytes - 1) & ~static_cast<std::uintptr_t>(bucket_bytes - 1)};
        blk.first = reinterpret_cast<bucket*>(aligned);
        blk.count = count;
        for (std::size_t i {0}; i < count; ++i)
        {
            put(blk.first + i);
        }
        blocks_.push_back(std::move(blk));
    }

    std::vector<block> blocks_;
    bucket* free_ {nullptr};
};

class erat_big
{
    // Dusart's bound on pi(x), used only for sizing
    static std::uint64_t prime_count_upper_bound_local(std::uint64_t x) noexcept
    {
        if (x < 60184)
        {
            return x / 8 + 64;
        }
        const double d {static_cast<double>(x)};
        return static_cast<std::uint64_t>(d / (std::log(d) - 1.1)) + 1;
    }

public:
    void init(std::uint64_t stop, std::uint64_t max_prime, std::size_t sieve_bytes, bool enabled)
    {
        stop_ = stop;
        enabled_ = enabled;
        staged_ = 0;
        pool_.release_all();
        lists_.clear();
        if (!enabled)
        {
            return;
        }
        sieve_bytes_ = sieve_bytes;
        log2_sieve_bytes_ = static_cast<unsigned>(63 - boost::math::tools::countl_zero(sieve_bytes));
        mask_ = sieve_bytes - 1;
        // Enough lists for the largest possible jump: 10 wheel steps of the largest prime
        const std::size_t lists {static_cast<std::size_t>(((max_prime / 30 + 1) * 10 + sieve_bytes) / sieve_bytes) + 2};
        lists_.assign(lists, nullptr);
        head_ = 0;
        staging_.reserve(static_cast<std::size_t>((std::min)(prime_count_upper_bound_local(max_prime), std::uint64_t(1u << 20))));
    }

    // Primes are staged unconditionally and the staging index advances by the accept flag:
    // for short windows at large magnitudes about half of the primes have no multiple in
    // range, and a branch on that outcome mispredicts (about 3 ns per prime measured).
    void add(std::uint64_t prime, std::uint64_t segment_low)
    {
        if (staged_ == staging_.size())
        {
            staging_.resize(staged_ + 65536);
        }
        staged_prime& st {staging_[staged_]};
        const bool hit {first_multiple(wheel210, prime, segment_low, stop_, st.multiple_index, st.wheel_index)};
        st.prime30 = static_cast<std::uint32_t>(prime / 30);
        staged_ += hit ? 1u : 0u;
    }

    void cross_off(std::uint8_t* sieve)
    {
        flush_staging();
        // A prime whose next multiple is still inside this segment is re-filed into the
        // current list, so keep draining it until nothing is left.
        while (lists_[head_] != nullptr)
        {
            sieving_prime* write {lists_[head_]};
            lists_[head_] = nullptr;
            bucket* b {bucket::owner(write)};
            b->end = write;
            while (b != nullptr)
            {
                for (const sieving_prime* sp {b->primes}; sp != b->end; ++sp)
                {
                    // the byte touched a few entries ahead is known: hide the cache miss
                    prefetch_byte(sieve + sp[prefetch_distance].multiple_index());
                    const auto& e = wheel210.wheel[sp->wheel_index()];
                    const std::size_t i {sp->multiple_index()};
                    sieve[i] &= e.unset_bit;
                    const std::uint64_t next {static_cast<std::uint64_t>(i) + static_cast<std::uint64_t>(sp->prime30) * e.next_multiple_factor + e.correct};
                    sieving_prime moved {*sp};
                    store(moved, next, e.next);
                }
                bucket* done {b};
                b = b->previous;
                if (b != nullptr)
                {
                    b->end = b->primes + bucket::capacity;
                }
                pool_.put(done);
            }
        }
        head_ = head_ + 1 == lists_.size() ? 0 : head_ + 1;
    }

private:
    struct staged_prime
    {
        std::uint64_t multiple_index;
        std::uint32_t prime30;
        unsigned wheel_index;
    };

    // Distributes the accepted staged primes into their lists (a predictable loop).
    void flush_staging()
    {
        for (std::size_t i {0}; i < staged_; ++i)
        {
            const staged_prime& st {staging_[i]};
            sieving_prime sp {};
            sp.prime30 = st.prime30;
            store(sp, st.multiple_index, st.wheel_index);
        }
        staged_ = 0;
    }

    // Files the prime in the list of the segment (relative to the current one) that holds
    // the byte index, keeping only the offset inside that segment.
    void store(sieving_prime& sp, std::uint64_t multiple_index, unsigned wheel_index)
    {
        const std::size_t segment {static_cast<std::size_t>(multiple_index >> log2_sieve_bytes_)};
        sp.set(static_cast<std::size_t>(multiple_index & mask_), wheel_index);
        // segment is below the list count by construction; the wrap compiles to a select
        const std::size_t slot {head_ + segment};
        const std::size_t wrapped {slot >= lists_.size() ? slot - lists_.size() : slot};
        sieving_prime*& write {lists_[wrapped]};
        if (bucket::is_full(write))
        {
            bucket* b {pool_.get()};
            b->previous = write == nullptr ? nullptr : bucket::owner(write);
            write = b->primes;
        }
        *write++ = sp;
    }

    std::vector<sieving_prime*> lists_;
    bucket_pool pool_;
    std::vector<staged_prime> staging_;
    std::size_t staged_ {0};
    std::size_t head_ {0};
    std::size_t sieve_bytes_ {0};
    std::size_t mask_ {0};
    unsigned log2_sieve_bytes_ {0};
    std::uint64_t stop_ {0};
    bool enabled_ {false};
};

} // namespace boost::math::detail::prime_sieve

#endif // BOOST_MATH_SF_DETAIL_PRIME_SIEVE_ERAT_BIG_HPP
