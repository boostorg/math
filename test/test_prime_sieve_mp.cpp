//  (C) Copyright Matt Borland 2026.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
//  Prime sieve with Boost.Multiprecision integers, including bounds beyond 2^64.
//  Define BOOST_MATH_TEST_GMP to also run with mpz_int.

#include <boost/math/special_functions/prime_sieve.hpp>
#include <boost/math/special_functions/detail/prime_sieve/primality.hpp>
#include <boost/multiprecision/cpp_int.hpp>
#include <boost/multiprecision/miller_rabin.hpp>
#ifdef BOOST_MATH_TEST_GMP
#include <boost/multiprecision/gmp.hpp>
#endif
#include <boost/core/lightweight_test.hpp>
#include <cstdint>
#include <vector>
#include <iterator>

using namespace boost::math;
namespace ps = boost::math::detail::prime_sieve;

template <class Integer>
void test_helpers()
{
    // Jacobi symbol against a brute-force Legendre symbol
    for (const int q : {3, 5, 7, 11, 13, 17, 19, 23, 101, 997})
    {
        for (int a {0}; a < q; ++a)
        {
            int legendre {0};
            for (int x {1}; x < q && legendre == 0; ++x)
            {
                if ((x * x) % q == a)
                {
                    legendre = 1;
                }
            }
            if (a != 0 && legendre == 0)
            {
                legendre = -1;
            }
            BOOST_TEST_EQ(ps::jacobi(Integer(a), Integer(q)), legendre);
        }
    }
    BOOST_TEST_EQ(ps::jacobi_signed(-1, Integer(13)), 1);
    BOOST_TEST_EQ(ps::jacobi_signed(-1, Integer(7)), -1);

    // powm against Boost.Multiprecision
    for (int i {1}; i < 50; ++i)
    {
        const Integer b {Integer(1234567) * i + 89};
        const Integer e {Integer(987654321) * i};
        const Integer m {Integer(1000000007) * i + 2};
        BOOST_TEST(ps::powm(b, e, m) == boost::multiprecision::powm(b, e, m));
    }
    // integer roots and perfect powers
    for (int i {1}; i < 50; ++i)
    {
        const Integer n {Integer(123456789) * Integer(987654321) * i};
        const Integer r {ps::isqrt(n)};
        BOOST_TEST(r * r <= n && (r + 1) * (r + 1) > n);
        const Integer c {ps::iroot(n, 3)};
        BOOST_TEST(c * c * c <= n && (c + 1) * (c + 1) * (c + 1) > n);
    }
    BOOST_TEST(ps::is_perfect_power(Integer(1194649), 3));      // 1093^2
    BOOST_TEST(ps::is_perfect_power(Integer(1194649), 2));
    BOOST_TEST(!ps::is_perfect_power(Integer(1194651), 7));
    BOOST_TEST(ps::is_perfect_power(Integer(Integer(3511) * 3511 * 3511), 3));

    // Baillie-PSW: pseudoprimes, Carmichael numbers, strong Lucas pseudoprimes, small primes
    for (const unsigned c : {561u, 1105u, 1729u, 2465u, 2821u, 2047u, 3277u, 4033u, 4681u, 8321u, 5459u, 5777u, 10877u, 16109u, 18971u})
    {
        BOOST_TEST(!ps::is_probable_prime_bpsw(Integer(c)));
    }
    BOOST_TEST(!ps::is_probable_prime_bpsw(Integer(3215031751ull)));
    for (unsigned n {0}; n < 3000; ++n)
    {
        BOOST_TEST_EQ(ps::is_probable_prime_bpsw(Integer(n)), ps::is_prime_u64(n));
    }
    // the pseudosquares themselves are classified like Miller-Rabin does
    for (std::size_t i {0}; i < 49; ++i)
    {
        const Integer L {ps::pseudosquare_value<Integer>(i)};
        BOOST_TEST_EQ(ps::is_probable_prime_bpsw(L), boost::multiprecision::miller_rabin_test(L, 25));
    }
    // Wieferich squares must fail the pseudosquares test through the perfect power check
    const Integer w1 {1194649};
    BOOST_TEST(!ps::pseudosquares_prime_test(w1, ps::pseudosquare_index_for(Integer(w1 / 1000)), 1000));
    const Integer w2 {Integer(3511) * 3511};
    BOOST_TEST(!ps::pseudosquares_prime_test(w2, ps::pseudosquare_index_for(Integer(w2 / 1000)), 1000));
}

template <class Integer>
void test_sieve()
{
    // results equal the builtin ones
    BOOST_TEST_EQ(prime_count(Integer(1000000)), 78498u);
    std::vector<Integer> s;
    prime_sieve(Integer(1000), s);
    BOOST_TEST(s.size() == 168 && s.back() == 997);
    BOOST_TEST_EQ(prime_count(Integer("1000000000000"), Integer("1000000100000")), prime_count(1000000000000ull, 1000000100000ull));
    BOOST_TEST_EQ(prime_count(Integer("1000000000000000000"), Integer("1000000000001000000")), 24280u);

    const Integer two64 {Integer(1) << 64};

    // primes around 2^64 through an output iterator and through a vector
    const char* expected[] =
    {
        "18446744073709551557", "18446744073709551629", "18446744073709551653", "18446744073709551667",
        "18446744073709551697", "18446744073709551709", "18446744073709551757", "18446744073709551923",
        "18446744073709551947", "18446744073709552009"
    };
    std::vector<Integer> v;
    prime_range(Integer(two64 - 60), Integer(two64 + 400), std::back_inserter(v));
    BOOST_TEST_EQ(v.size(), 10u);
    for (std::size_t i {0}; i < v.size() && i < 10; ++i)
    {
        BOOST_TEST(v[i] == Integer(expected[i]));
    }
    std::vector<Integer> w;
    prime_range(Integer(two64 - 60), Integer(two64 + 400), w);
    BOOST_TEST(v == w);
    BOOST_TEST_EQ(prime_count(Integer(two64 - 60), Integer(two64 + 400)), 10u);

    // a range straddling 2^64 is the sum of its two halves
    std::vector<Integer> straddle;
    prime_range(Integer(two64 - 5000), Integer(two64 + 5000), straddle);
    const std::uint64_t below {prime_count(18446744073709551615ull - 4999, 18446744073709551615ull)};
    const std::uint64_t above {prime_count(two64, Integer(two64 + 5000))};
    BOOST_TEST_EQ(straddle.size(), below + above);

    // reference counts in windows beyond 2^64 (pseudosquares path and BPSW path)
    struct reference
    {
        const char* lower;
        std::uint64_t width;
        std::uint64_t count;
    };
    const reference references[] =
    {
        {"18446744073709551616", 200000, 4335},
        {"100000000000000000000", 200000, 4294},
        {"10000000000000000000000", 100000, 1979},
        {"100000000000000000000000", 50000, 942},
        {"1000000000000000000000000", 50000, 919}
    };
    prime_sieve_options probable {};
    probable.probable_prime_only = true;
    for (const reference& r : references)
    {
        const Integer lo {r.lower};
        const Integer hi {lo + r.width};
        BOOST_TEST_EQ(prime_count(lo, hi), r.count);
        BOOST_TEST_EQ(prime_count(lo, hi, probable), r.count);
        std::vector<Integer> primes;
        prime_range(lo, hi, primes);
        BOOST_TEST_EQ(primes.size(), r.count);
        bool ok {true};
        for (std::size_t i {0}; i < primes.size(); ++i)
        {
            if (!boost::multiprecision::miller_rabin_test(primes[i], 25) || (i > 0 && !(primes[i - 1] < primes[i])))
            {
                ok = false;
            }
        }
        BOOST_TEST(ok);
    }
}

int main()
{
    test_helpers<boost::multiprecision::cpp_int>();
    test_sieve<boost::multiprecision::cpp_int>();
#ifdef BOOST_MATH_TEST_GMP
    test_helpers<boost::multiprecision::mpz_int>();
    test_sieve<boost::multiprecision::mpz_int>();
#endif
    return boost::report_errors();
}
