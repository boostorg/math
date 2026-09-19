//  (C) Copyright Matt Borland 2026.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
//  Medium sieving primes (a handful of multiples per segment). Primes are kept in 64 lists
//  keyed by wheel state; each list is processed by one function specialized for that state
//  (Duff's device over the eight phases), so the branch predictor sees one pattern and the
//  per-prime overhead is a loop iteration rather than a call. After a segment each prime is
//  re-filed under its new state.

#ifndef BOOST_MATH_SF_DETAIL_PRIME_SIEVE_ERAT_MEDIUM_HPP
#define BOOST_MATH_SF_DETAIL_PRIME_SIEVE_ERAT_MEDIUM_HPP

#include <boost/math/tools/config.hpp>
#include <boost/math/special_functions/detail/prime_sieve/layout.hpp>

#ifndef BOOST_MATH_BUILD_MODULE
#include <cstdint>
#include <cstddef>
#include <array>
#include <vector>
#include <utility>
#endif

namespace boost::math::detail::prime_sieve {

namespace medium_detail {

using list_array = std::array<std::vector<sieving_prime>, 64>;

// Crosses off the multiples of every prime in [first, last) (all in wheel state State) over
// sieve[0, n) and files each prime, with its offset into the next segment, into next.
template <unsigned State>
void cross_off_list(std::uint8_t* sieve, std::size_t n, const sieving_prime* first, const sieving_prime* last, list_array& next) noexcept
{
    constexpr unsigned base {(State / 8) * 8};
    constexpr const auto& w = wheel30.wheel;

    for (const sieving_prime* sp {first}; sp != last; ++sp)
    {
        const std::size_t p {static_cast<std::size_t>(sp->prime30)};
        std::size_t i {sp->multiple_index()};
        unsigned state {State};
        if (i < n)
        {
            const std::size_t d0 {p * w[base + 0].next_multiple_factor + w[base + 0].correct};
            const std::size_t d1 {p * w[base + 1].next_multiple_factor + w[base + 1].correct};
            const std::size_t d2 {p * w[base + 2].next_multiple_factor + w[base + 2].correct};
            const std::size_t d3 {p * w[base + 3].next_multiple_factor + w[base + 3].correct};
            const std::size_t d4 {p * w[base + 4].next_multiple_factor + w[base + 4].correct};
            const std::size_t d5 {p * w[base + 5].next_multiple_factor + w[base + 5].correct};
            const std::size_t d6 {p * w[base + 6].next_multiple_factor + w[base + 6].correct};
            const std::size_t d7 {p * w[base + 7].next_multiple_factor + w[base + 7].correct};
            switch (State % 8)
            {
                for (;;)
                {
                case 0:
                    if (i >= n) { state = base + 0; break; }
                    sieve[i] &= w[base + 0].unset_bit;
                    i += d0;
                    [[fallthrough]];
                case 1:
                    if (i >= n) { state = base + 1; break; }
                    sieve[i] &= w[base + 1].unset_bit;
                    i += d1;
                    [[fallthrough]];
                case 2:
                    if (i >= n) { state = base + 2; break; }
                    sieve[i] &= w[base + 2].unset_bit;
                    i += d2;
                    [[fallthrough]];
                case 3:
                    if (i >= n) { state = base + 3; break; }
                    sieve[i] &= w[base + 3].unset_bit;
                    i += d3;
                    [[fallthrough]];
                case 4:
                    if (i >= n) { state = base + 4; break; }
                    sieve[i] &= w[base + 4].unset_bit;
                    i += d4;
                    [[fallthrough]];
                case 5:
                    if (i >= n) { state = base + 5; break; }
                    sieve[i] &= w[base + 5].unset_bit;
                    i += d5;
                    [[fallthrough]];
                case 6:
                    if (i >= n) { state = base + 6; break; }
                    sieve[i] &= w[base + 6].unset_bit;
                    i += d6;
                    [[fallthrough]];
                case 7:
                    if (i >= n) { state = base + 7; break; }
                    sieve[i] &= w[base + 7].unset_bit;
                    i += d7;
                }
            }
        }
        sieving_prime moved {*sp};
        moved.set(i - n, state);
        next[state].push_back(moved);
    }
}

using cross_off_fn = void (*)(std::uint8_t*, std::size_t, const sieving_prime*, const sieving_prime*, list_array&) noexcept;

// Indirection keeps the pack expansion simple enough for every front end (nvcc included)
template <unsigned State>
constexpr cross_off_fn list_function() noexcept
{
    return &cross_off_list<State>;
}

template <std::size_t... I>
constexpr std::array<cross_off_fn, 64> make_table(std::index_sequence<I...>) noexcept
{
    return {{list_function<static_cast<unsigned>(I)>()...}};
}

inline constexpr std::array<cross_off_fn, 64> table = make_table(std::make_index_sequence<64>{});

} // namespace medium_detail

class erat_medium
{
public:
    void init(std::uint64_t stop)
    {
        stop_ = stop;
        for (auto& l : lists_)
        {
            l.clear();
        }
        for (auto& l : next_)
        {
            l.clear();
        }
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
            lists_[wheel_index].push_back(sp);
        }
    }

    void cross_off(std::uint8_t* sieve, std::size_t bytes)
    {
        for (unsigned w {0}; w < 64; ++w)
        {
            const std::vector<sieving_prime>& list {lists_[w]};
            if (!list.empty())
            {
                medium_detail::table[w](sieve, bytes, list.data(), list.data() + list.size(), next_);
            }
        }
        for (auto& l : lists_)
        {
            l.clear();
        }
        lists_.swap(next_);
    }

private:
    medium_detail::list_array lists_ {};
    medium_detail::list_array next_ {};
    std::uint64_t stop_ {0};
};

} // namespace boost::math::detail::prime_sieve

#endif // BOOST_MATH_SF_DETAIL_PRIME_SIEVE_ERAT_MEDIUM_HPP
