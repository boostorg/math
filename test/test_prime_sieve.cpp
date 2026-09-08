//  (C) Copyright Matt Borland 2026.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <boost/math/special_functions/prime_sieve.hpp>
#include <boost/math/special_functions/detail/prime_sieve/primality.hpp>
#include <boost/math/special_functions/prime.hpp>
#include <boost/core/lightweight_test.hpp>
#include <cstdint>
#include <cstddef>
#include <vector>
#include <list>
#include <iterator>
#include <random>
#include <limits>

using namespace boost::math;
namespace ps = boost::math::detail::prime_sieve;

// Known values of the prime counting function
struct pi_entry
{
    std::uint64_t x;
    std::uint64_t pi;
};

constexpr pi_entry pi_table[] =
{
    {10, 4}, {100, 25}, {1000, 168}, {10000, 1229}, {100000, 9592}, {1000000, 78498},
    {10000000, 664579}, {100000000, 5761455}, {1000000000, 50847534}
};

void test_primality_helpers()
{
    // The two 64-bit tests agree with each other, with brute force, and with the prime table
    for (unsigned n {0}; n < 5000; ++n)
    {
        bool brute {n >= 2};
        for (unsigned f {2}; f * f <= n; ++f)
        {
            if (n % f == 0)
            {
                brute = false;
                break;
            }
        }
        BOOST_TEST_EQ(ps::is_prime_u64(n), brute);
        BOOST_TEST_EQ(ps::is_prime_u64_miller_rabin(n), brute);
    }
    for (unsigned i {0}; i < boost::math::max_prime; i += 97)
    {
        BOOST_TEST(ps::is_prime_u64(boost::math::prime(i)));
    }
    // Carmichael numbers and strong pseudoprimes to base 2
    for (const std::uint64_t c : {561ull, 1105ull, 1729ull, 2465ull, 2821ull, 2047ull, 3277ull, 4033ull, 4681ull, 8321ull, 3215031751ull})
    {
        BOOST_TEST(!ps::is_prime_u64(c));
    }
    BOOST_TEST(ps::is_prime_u64(18446744073709551557ull));
    BOOST_TEST(!ps::is_prime_u64(18446744073709551615ull));
    BOOST_TEST(!ps::is_prime_u64(18446744073709551559ull));
    std::mt19937_64 rng {12345};
    for (int i {0}; i < 20000; ++i)
    {
        const std::uint64_t n {rng()};
        BOOST_TEST_EQ(ps::is_prime_u64(n), ps::is_prime_u64_miller_rabin(n));
    }
}

template <class Integer>
void test_pi_table()
{
    for (const pi_entry& e : pi_table)
    {
        if (e.x > static_cast<std::uint64_t>((std::numeric_limits<Integer>::max)()))
        {
            continue;
        }
        const Integer x {static_cast<Integer>(e.x)};
        BOOST_TEST_EQ(prime_count(x), e.pi);

        std::vector<Integer> v;
        prime_reserve(x, v);
        BOOST_TEST(v.capacity() >= e.pi);
        prime_sieve(x, v);
        BOOST_TEST_EQ(v.size(), e.pi);
        for (std::size_t i {1}; i < v.size(); ++i)
        {
            if (!(v[i - 1] < v[i]))
            {
                BOOST_ERROR("output not ascending");
                break;
            }
        }
        if (!v.empty())
        {
            BOOST_TEST_EQ(v.front(), Integer(2));
        }

        std::vector<Integer> w;
        prime_sieve(x, std::back_inserter(w));
        BOOST_TEST(v == w);
    }
}

// Every value in [lo, hi) is checked individually against the primality test
void check_range(std::uint64_t lo, std::uint64_t hi, const prime_sieve_options& options = {})
{
    std::vector<std::uint64_t> v;
    prime_range(lo, hi, v, options);
    std::uint64_t expected {0};
    std::size_t k {0};
    bool ok {true};
    for (std::uint64_t n {lo}; n < hi; ++n)
    {
        if (ps::is_prime_u64(n))
        {
            ++expected;
            if (k >= v.size() || v[k] != n)
            {
                ok = false;
            }
            ++k;
        }
        if (n == (std::numeric_limits<std::uint64_t>::max)())
        {
            break;
        }
    }
    BOOST_TEST(ok);
    BOOST_TEST_EQ(v.size(), expected);
    BOOST_TEST_EQ(prime_count(lo, hi, options), expected);
}

void test_ranges()
{
    std::mt19937_64 rng {7};
    const std::uint64_t lows[] =
    {
        0, 1, 2, 3, 6, 7, 8, 29, 30, 31, 32, 1000, 65536, 1000000, 1000000000,
        4294967296ull - 50000, 1000000000000ull, 1000000000000000ull, 1000000000000000000ull,
        (std::numeric_limits<std::uint64_t>::max)() - 100000
    };
    for (const std::uint64_t lo : lows)
    {
        for (int rep {0}; rep < 2; ++rep)
        {
            const std::uint64_t width {rep == 0 ? 100000u : static_cast<std::uint64_t>(rng() % 100000)};
            const std::uint64_t hi {lo > (std::numeric_limits<std::uint64_t>::max)() - width ? (std::numeric_limits<std::uint64_t>::max)() : lo + width};
            check_range(lo, hi);
        }
    }
    // every residue of the lower bound modulo 30
    for (std::uint64_t r {0}; r < 30; ++r)
    {
        check_range(1000000 + r, 1000000 + r + 3000);
    }
    // both strategies must agree on a short interval at a large magnitude
    prime_sieve_options full {};
    full.range_strategy = prime_range_strategy::full_sieve;
    prime_sieve_options test {};
    test.range_strategy = prime_range_strategy::test_survivors;
    BOOST_TEST_EQ(prime_count(1000000000000000000ull, 1000000000001000000ull, full), 24280u);
    BOOST_TEST_EQ(prime_count(1000000000000000000ull, 1000000000001000000ull, test), 24280u);
    BOOST_TEST_EQ(prime_count(1000000000000000000ull, 1000000000001000000ull), 24280u);
    // forced geometries: a small segment engages the bucket sieve, a large one the medium class only
    prime_sieve_options small_segment {};
    small_segment.sieve_bytes = 16384;
    BOOST_TEST_EQ(prime_count(1000000000ull, small_segment), 50847534u);
    check_range(1000000000000ull, 1000000000000ull + 300000, small_segment);
    prime_sieve_options large_segment {};
    large_segment.sieve_bytes = 8u * 1024u * 1024u;
    BOOST_TEST_EQ(prime_count(100000000ull, large_segment), 5761455u);
    prime_sieve_options small_l1 {};
    small_l1.l1d_bytes = 16384;
    BOOST_TEST_EQ(prime_count(100000000ull, small_l1), 5761455u);
    // windows that need the bucket sieve across many segments
    BOOST_TEST_EQ(prime_count(1000000000000000ull, 1000000000000000ull + 100000000ull), 2893937u);
}

void test_edge_cases()
{
    std::vector<int> s;
    prime_range(0, 2, s);
    BOOST_TEST(s.empty());
    prime_range(0, 3, s);
    BOOST_TEST(s.size() == 1 && s[0] == 2);
    s.clear();
    prime_range(2, 8, s);
    BOOST_TEST_EQ(s.size(), 4u);
    s.clear();
    prime_range(7, 8, s);
    BOOST_TEST(s.size() == 1 && s[0] == 7);
    s.clear();
    prime_range(30, 31, s);
    BOOST_TEST(s.empty());
    prime_range(31, 32, s);
    BOOST_TEST(s.size() == 1 && s[0] == 31);
    s.clear();
    prime_range(100, 100, s);
    BOOST_TEST(s.empty());
    prime_range(100, 50, s);
    BOOST_TEST(s.empty());
    prime_range(-100, 10, s);
    BOOST_TEST_EQ(s.size(), 4u);
    s.clear();
    prime_sieve(-5, s);
    BOOST_TEST(s.empty());

    BOOST_TEST_EQ(prime_count(0), 0u);
    BOOST_TEST_EQ(prime_count(2), 0u);
    BOOST_TEST_EQ(prime_count(3), 1u);
    BOOST_TEST_EQ(prime_count(8), 4u);
    BOOST_TEST_EQ(prime_count(168), 39u);
    BOOST_TEST_EQ(prime_count(-7), 0u);
    BOOST_TEST_EQ(prime_count(5, 5), 0u);
    BOOST_TEST_EQ(prime_count(5, 6), 1u);
    BOOST_TEST_EQ(prime_count(6, 7), 0u);

    // the largest 64-bit prime and the top of the range
    std::vector<std::uint64_t> top;
    prime_range(static_cast<std::uint64_t>(18446744073709551557ull), (std::numeric_limits<std::uint64_t>::max)(), top);
    BOOST_TEST(top.size() == 1 && top[0] == 18446744073709551557ull);
    BOOST_TEST_EQ(prime_count(static_cast<std::uint64_t>(18446744073709551558ull), (std::numeric_limits<std::uint64_t>::max)()), 0u);

    // other output iterators
    std::list<long long> l;
    prime_range(10, 100, std::back_inserter(l));
    BOOST_TEST_EQ(l.size(), 21u);
    BOOST_TEST_EQ(l.front(), 11);
    BOOST_TEST_EQ(l.back(), 97);
    std::vector<unsigned> buffer(200);
    unsigned* end {prime_sieve(1000u, buffer.data())};
    BOOST_TEST_EQ(static_cast<std::size_t>(end - buffer.data()), 168u);
    BOOST_TEST_EQ(buffer[167], 997u);

    // narrower integer types
    BOOST_TEST_EQ(prime_count(static_cast<std::int16_t>(30000)), 3245u);
    std::vector<std::uint8_t> bytes;
    prime_sieve(static_cast<std::uint8_t>(255), bytes);
    BOOST_TEST_EQ(bytes.size(), 54u);
}

void test_policies()
{
#ifdef BOOST_MATH_PRIME_SIEVE_HAS_STD_EXECUTION
    for (const std::uint64_t n : {1000ull, 100000ull, 100000000ull, 1000000000ull})
    {
        std::vector<std::uint64_t> a;
        std::vector<std::uint64_t> b;
        prime_sieve(std::execution::seq, n, a);
        prime_sieve(std::execution::par, n, b);
        BOOST_TEST(a == b);
        BOOST_TEST_EQ(prime_count(std::execution::seq, n), a.size());
        BOOST_TEST_EQ(prime_count(std::execution::par, n), a.size());
        BOOST_TEST_EQ(prime_count(std::execution::par_unseq, n), a.size());
    }
    for (const std::uint64_t lo : {0ull, 1000000ull, 1000000000000ull, 1000000000000000000ull, (std::numeric_limits<std::uint64_t>::max)() - 400000000ull})
    {
        const std::uint64_t hi {lo > (std::numeric_limits<std::uint64_t>::max)() - 300000000ull ? (std::numeric_limits<std::uint64_t>::max)() : lo + 300000000ull};
        std::vector<std::uint64_t> a;
        std::vector<std::uint64_t> b;
        prime_range(std::execution::seq, lo, hi, std::back_inserter(a));
        prime_range(std::execution::par, lo, hi, b);
        BOOST_TEST(a == b);
        BOOST_TEST_EQ(prime_count(std::execution::par, lo, hi), a.size());
        std::list<std::uint64_t> l;
        prime_range(std::execution::par, lo, hi, std::back_inserter(l));
        BOOST_TEST_EQ(l.size(), a.size());
    }
    prime_sieve_options three {};
    three.max_threads = 3;
    BOOST_TEST_EQ(prime_count(std::execution::par, 1000000000ull, three), 50847534u);
    prime_sieve_options one {};
    one.max_threads = 1;
    BOOST_TEST_EQ(prime_count(std::execution::par, 100000000ull, one), 5761455u);
    std::vector<int> vi;
    prime_sieve(std::execution::par_unseq, 1000, vi);
    BOOST_TEST(vi.size() == 168 && vi.front() == 2 && vi.back() == 997);
#endif
    // the CUDA tag is always accepted; without nvcc it runs on the CPU
    BOOST_TEST_EQ(prime_count(execution::cuda, 1000000ull), 78498u);
    std::vector<std::uint64_t> vc;
    prime_sieve(execution::cuda, 100000ull, vc);
    BOOST_TEST_EQ(vc.size(), 9592u);
    prime_range(execution::cuda, 1000000000000ull, 1000000001000ull, vc);
    BOOST_TEST_EQ(vc.size(), 9592u + 37u);
}

int main()
{
    test_primality_helpers();
    test_pi_table<int>();
    test_pi_table<unsigned>();
    test_pi_table<long long>();
    test_pi_table<std::uint64_t>();
    test_pi_table<std::int32_t>();
    test_ranges();
    test_edge_cases();
    test_policies();
    return boost::report_errors();
}
