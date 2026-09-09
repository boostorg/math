//  (C) Copyright Matt Borland 2026.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
//  Small sieving primes (many multiples per segment). The segment is processed in L1-sized
//  chunks; for each prime one full wheel turn (8 multiples, exactly p bytes) is unrolled, and
//  the state is carried into the next chunk.

#ifndef BOOST_MATH_SF_DETAIL_PRIME_SIEVE_ERAT_SMALL_HPP
#define BOOST_MATH_SF_DETAIL_PRIME_SIEVE_ERAT_SMALL_HPP

#include <boost/math/tools/config.hpp>
#include <boost/math/special_functions/detail/prime_sieve/layout.hpp>

#ifndef BOOST_MATH_BUILD_MODULE
#include <cstdint>
#include <cstddef>
#include <array>
#include <vector>
#include <utility>
#include <algorithm>
#endif

namespace boost::math::detail::prime_sieve {

namespace small_detail {

// Cumulative byte offsets, bit masks and wheel states of the eight multiples of one wheel turn
// starting at State.
struct turn_chain
{
    std::array<unsigned, 8> factor {};   // multiply by p / 30
    std::array<unsigned, 8> correct {};  // add
    std::array<std::uint8_t, 8> mask {};
    std::array<unsigned, 8> state {};
};

template <unsigned State>
constexpr turn_chain make_chain() noexcept
{
    turn_chain c {};
    unsigned s {State};
    unsigned f {0};
    unsigned k {0};
    for (unsigned j {0}; j < 8; ++j)
    {
        c.factor[j] = f;
        c.correct[j] = k;
        c.mask[j] = wheel30.wheel[s].unset_bit;
        c.state[j] = s;
        f += wheel30.wheel[s].next_multiple_factor;
        k += wheel30.wheel[s].correct;
        s = wheel30.wheel[s].next;
    }
    return c;
}

// Full wheel turns of a prime in state State over sieve[i, n): eight multiples per iteration
// with one bounds check. Leaves i at the first multiple that does not have a whole turn of room.
template <unsigned State>
BOOST_MATH_FORCEINLINE void turn_loop(std::uint8_t* sieve, std::size_t n, std::size_t& i, std::size_t p) noexcept
{
    constexpr turn_chain c {make_chain<State>()};
    constexpr unsigned cls {State / 8};
    const std::size_t turn {30u * p + class_residue(cls)};
    const std::size_t o1 {c.factor[1] * p + c.correct[1]};
    const std::size_t o2 {c.factor[2] * p + c.correct[2]};
    const std::size_t o3 {c.factor[3] * p + c.correct[3]};
    const std::size_t o4 {c.factor[4] * p + c.correct[4]};
    const std::size_t o5 {c.factor[5] * p + c.correct[5]};
    const std::size_t o6 {c.factor[6] * p + c.correct[6]};
    const std::size_t o7 {c.factor[7] * p + c.correct[7]};
    if (n > o7)
    {
        const std::size_t limit {n - o7};
        std::size_t j {i};
        for (; j < limit; j += turn)
        {
            sieve[j] &= c.mask[0];
            sieve[j + o1] &= c.mask[1];
            sieve[j + o2] &= c.mask[2];
            sieve[j + o3] &= c.mask[3];
            sieve[j + o4] &= c.mask[4];
            sieve[j + o5] &= c.mask[5];
            sieve[j + o6] &= c.mask[6];
            sieve[j + o7] &= c.mask[7];
        }
        i = j;
    }
}

// One step of the chain of class C at phase K: exits with the state when the multiple leaves
// the chunk, otherwise clears the bit and advances to the next phase.
#define BOOST_MATH_PRIME_SIEVE_STEP_BODY(C, K)                                    \
    if (i >= n)                                                                                    \
    {                                                                                              \
        state = (C) * 8 + (K);                                                                     \
        goto done;                                                                                 \
    }                                                                                              \
    sieve[i] &= wheel30.wheel[(C) * 8 + (K)].unset_bit;                                            \
    i += p * wheel30.wheel[(C) * 8 + (K)].next_multiple_factor + wheel30.wheel[(C) * 8 + (K)].correct;

#define BOOST_MATH_PRIME_SIEVE_STEP(C, K)                                         \
    step_##C##_##K:                                                                \
    BOOST_MATH_PRIME_SIEVE_STEP_BODY(C, K)

// Entry points of class C: phase 0 goes straight to the unrolled loop, other phases take
// single steps down the chain until it wraps to phase 0.
#define BOOST_MATH_PRIME_SIEVE_CASES(C)                                                           \
    case (C) * 8 + 0: goto loop_##C;                                                               \
    case (C) * 8 + 1: goto step_##C##_1;                                                           \
    case (C) * 8 + 2: goto step_##C##_2;                                                           \
    case (C) * 8 + 3: goto step_##C##_3;                                                           \
    case (C) * 8 + 4: goto step_##C##_4;                                                           \
    case (C) * 8 + 5: goto step_##C##_5;                                                           \
    case (C) * 8 + 6: goto step_##C##_6;                                                           \
    case (C) * 8 + 7: goto step_##C##_7;

// The body of class C: whole turns from phase 0, then one chain of single steps that wraps.
// One unrolled loop per class (eight in total) keeps the function within the micro-op cache;
// a loop per state (64) measured 20 percent slower on x86.
#define BOOST_MATH_PRIME_SIEVE_CHAIN(C)                                                           \
    loop_##C:                                                                                      \
    turn_loop<(C) * 8>(sieve, n, i, p);                                                            \
    BOOST_MATH_PRIME_SIEVE_STEP_BODY(C, 0)                                                         \
    BOOST_MATH_PRIME_SIEVE_STEP(C, 1)                                                              \
    BOOST_MATH_PRIME_SIEVE_STEP(C, 2)                                                              \
    BOOST_MATH_PRIME_SIEVE_STEP(C, 3)                                                              \
    BOOST_MATH_PRIME_SIEVE_STEP(C, 4)                                                              \
    BOOST_MATH_PRIME_SIEVE_STEP(C, 5)                                                              \
    BOOST_MATH_PRIME_SIEVE_STEP(C, 6)                                                              \
    BOOST_MATH_PRIME_SIEVE_STEP(C, 7)                                                              \
    goto loop_##C;

// Crosses off every prime of [first, last) over the chunk sieve[0, n) and stores each prime's
// offset into the next chunk. One function for all 64 states (the switch is a jump table).
inline void cross_off_primes(std::uint8_t* sieve, std::size_t n, sieving_prime* first, sieving_prime* last) noexcept
{
    for (sieving_prime* sp {first}; sp != last; ++sp)
    {
        std::size_t i {sp->multiple_index()};
        unsigned state {sp->wheel_index()};
        if (i >= n)
        {
            sp->set(i - n, state);
            continue;
        }
        const std::size_t p {static_cast<std::size_t>(sp->prime30)};
        switch (state)
        {
            BOOST_MATH_PRIME_SIEVE_CASES(0)
            BOOST_MATH_PRIME_SIEVE_CASES(1)
            BOOST_MATH_PRIME_SIEVE_CASES(2)
            BOOST_MATH_PRIME_SIEVE_CASES(3)
            BOOST_MATH_PRIME_SIEVE_CASES(4)
            BOOST_MATH_PRIME_SIEVE_CASES(5)
            BOOST_MATH_PRIME_SIEVE_CASES(6)
            BOOST_MATH_PRIME_SIEVE_CASES(7)
            default: goto done;
        }
        BOOST_MATH_PRIME_SIEVE_CHAIN(0)
        BOOST_MATH_PRIME_SIEVE_CHAIN(1)
        BOOST_MATH_PRIME_SIEVE_CHAIN(2)
        BOOST_MATH_PRIME_SIEVE_CHAIN(3)
        BOOST_MATH_PRIME_SIEVE_CHAIN(4)
        BOOST_MATH_PRIME_SIEVE_CHAIN(5)
        BOOST_MATH_PRIME_SIEVE_CHAIN(6)
        BOOST_MATH_PRIME_SIEVE_CHAIN(7)
    done:
        sp->set(i - n, state);
    }
}

#undef BOOST_MATH_PRIME_SIEVE_CHAIN
#undef BOOST_MATH_PRIME_SIEVE_CASES
#undef BOOST_MATH_PRIME_SIEVE_STEP
#undef BOOST_MATH_PRIME_SIEVE_STEP_BODY

} // namespace small_detail

class erat_small
{
public:
    void init(std::uint64_t stop, std::size_t chunk_bytes)
    {
        stop_ = stop;
        chunk_bytes_ = chunk_bytes;
        primes_.clear();
    }

    void add(std::uint64_t prime, std::uint64_t segment_low)
    {
        std::uint64_t multiple_index {};
        unsigned wheel_index {};
        if (first_multiple(wheel30, prime, segment_low, stop_, multiple_index, wheel_index))
        {
            sieving_prime sp {};
            sp.prime30 = static_cast<std::uint32_t>(prime / 30);
            sp.set(static_cast<std::size_t>(multiple_index), wheel_index);
            primes_.push_back(sp);
        }
    }

    void cross_off(std::uint8_t* sieve, std::size_t bytes) noexcept
    {
        for (std::size_t chunk {0}; chunk < bytes; chunk += chunk_bytes_)
        {
            const std::size_t n {(std::min)(chunk_bytes_, bytes - chunk)};
            std::uint8_t* s {sieve + chunk};
            small_detail::cross_off_primes(s, n, primes_.data(), primes_.data() + primes_.size());
        }
    }

    bool empty() const noexcept
    {
        return primes_.empty();
    }

private:
    std::vector<sieving_prime> primes_;
    std::uint64_t stop_ {0};
    std::size_t chunk_bytes_ {32768};
};

} // namespace boost::math::detail::prime_sieve

#endif // BOOST_MATH_SF_DETAIL_PRIME_SIEVE_ERAT_SMALL_HPP
