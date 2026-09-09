//  (C) Copyright Matt Borland 2026.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
//  Primality tests used by the prime sieve:
//    * deterministic Miller-Rabin for 64-bit integers
//    * generic modular arithmetic, Jacobi symbol, integer roots
//    * Baillie-PSW probable prime test
//    * the Lukes-Patterson-Williams pseudosquares test (deterministic, table bounded)
//  The generic templates require an Integer type that can hold m * m for the modulus m
//  in use (Boost.Multiprecision integers qualify). Builtin 64-bit values use the
//  dedicated *_u64 functions.

#ifndef BOOST_MATH_SF_DETAIL_PRIME_SIEVE_PRIMALITY_HPP
#define BOOST_MATH_SF_DETAIL_PRIME_SIEVE_PRIMALITY_HPP

#include <boost/math/tools/config.hpp>

#ifndef BOOST_MATH_HAS_NVRTC

#include <boost/math/tools/bit.hpp>
#ifndef BOOST_MATH_BUILD_MODULE
#include <cstdint>
#include <cstddef>
#include <utility>
#include <type_traits>
#if defined(_MSC_VER) && !defined(__clang__) && defined(_M_X64)
#  include <intrin.h>
#endif
#endif

namespace boost::math::detail::prime_sieve {

// The 71 primes not exceeding 353: bases for the pseudosquares test and trial division.
inline constexpr std::uint16_t small_primes_to_353[71] =
{
      2,   3,   5,   7,  11,  13,  17,  19,  23,  29,  31,  37,  41,  43,  47,  53,  59,  61,  67,  71,
     73,  79,  83,  89,  97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173,
    179, 181, 191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281,
    283, 293, 307, 311, 313, 317, 331, 337, 347, 349, 353
};

// (a + b) mod m for a, b < m without overflow.
inline std::uint64_t addmod_u64(std::uint64_t a, std::uint64_t b, std::uint64_t m) noexcept
{
    return a >= m - b ? a - (m - b) : a + b;
}

// (a * b) mod m for any 64-bit operands.
inline std::uint64_t mulmod_u64(std::uint64_t a, std::uint64_t b, std::uint64_t m) noexcept
{
#if defined(BOOST_MATH_HAS_INT128)
    return static_cast<std::uint64_t>((static_cast<unsigned __int128>(a) * b) % m);
#elif defined(_MSC_VER) && !defined(__clang__) && defined(_M_X64)
    std::uint64_t high {};
    const std::uint64_t low {_umul128(a, b, &high)};
    if (high >= m)
    {
        high %= m;
    }
    std::uint64_t remainder {};
    _udiv128(high, low, m, &remainder);
    return remainder;
#else
    a %= m;
    b %= m;
    std::uint64_t result {0};
    while (b != 0)
    {
        if (b & 1u)
        {
            result = addmod_u64(result, a, m);
        }
        a = addmod_u64(a, a, m);
        b >>= 1;
    }
    return result;
#endif
}

// base^exp mod m.
inline std::uint64_t powm_u64(std::uint64_t base, std::uint64_t exp, std::uint64_t m) noexcept
{
    std::uint64_t result {1 % m};
    base %= m;
    while (exp != 0)
    {
        if (exp & 1u)
        {
            result = mulmod_u64(result, base, m);
        }
        base = mulmod_u64(base, base, m);
        exp >>= 1;
    }
    return result;
}

// Deterministic Miller-Rabin for all 64-bit values with the seven bases of Sinclair.
// Slower than is_prime_u64; kept as an independent oracle for testing.
inline bool is_prime_u64_miller_rabin(std::uint64_t n) noexcept
{
    if (n < 2)
    {
        return false;
    }
    for (std::size_t i {0}; i < 12; ++i)
    {
        const std::uint64_t p {small_primes_to_353[i]};
        if (n % p == 0)
        {
            return n == p;
        }
    }
    if (n < 41u * 41u)
    {
        return true;
    }

    constexpr std::uint64_t bases[7] = {2u, 325u, 9375u, 28178u, 450775u, 9780504u, 1795265022u};
    const int r {boost::math::tools::countr_zero(n - 1)};
    const std::uint64_t d {(n - 1) >> r};

    for (const std::uint64_t base : bases)
    {
        const std::uint64_t a {base % n};
        if (a == 0)
        {
            continue;
        }
        std::uint64_t x {powm_u64(a, d, n)};
        if (x == 1 || x == n - 1)
        {
            continue;
        }
        bool composite {true};
        for (int i {1}; i < r; ++i)
        {
            x = mulmod_u64(x, x, n);
            if (x == n - 1)
            {
                composite = false;
                break;
            }
        }
        if (composite)
        {
            return false;
        }
    }
    return true;
}

// Jacobi symbol (a / n) for odd n > 0.
inline int jacobi_u64(std::uint64_t a, std::uint64_t n) noexcept
{
    a %= n;
    int result {1};
    while (a != 0)
    {
        const int twos {boost::math::tools::countr_zero(a)};
        a >>= twos;
        if ((twos & 1) != 0)
        {
            const std::uint64_t r {n % 8};
            if (r == 3 || r == 5)
            {
                result = -result;
            }
        }
        const std::uint64_t t {a};
        a = n;
        n = t;
        if (a % 4 == 3 && n % 4 == 3)
        {
            result = -result;
        }
        a %= n;
    }
    return n == 1 ? result : 0;
}

// Montgomery arithmetic modulo an odd 64-bit n. Where no 128-bit product is available
// the "Montgomery" form degenerates to plain residues with a slow multiply.
class montgomery_u64
{
public:
    explicit montgomery_u64(std::uint64_t n) noexcept : n_(n)
    {
#if defined(BOOST_MATH_HAS_INT128) || (defined(_MSC_VER) && !defined(__clang__) && defined(_M_X64))
        // Newton iteration for n^-1 mod 2^64: n is odd so n * n == 1 (mod 8) seeds 3 bits
        std::uint64_t inv {n};
        for (int i {0}; i < 5; ++i)
        {
            inv *= 2u - n * inv;
        }
        n_inv_neg_ = 0u - inv;
        one_ = (0u - n) % n;                 // 2^64 mod n
        r2_ = mulmod_u64(one_, one_, n);     // 2^128 mod n
#else
        one_ = 1;
        r2_ = 1;
#endif
    }

    std::uint64_t modulus() const noexcept
    {
        return n_;
    }

    std::uint64_t one() const noexcept
    {
        return one_;
    }

    std::uint64_t to_form(std::uint64_t a) const noexcept
    {
        return mul(a % n_, r2_);
    }

    std::uint64_t add(std::uint64_t a, std::uint64_t b) const noexcept
    {
        return addmod_u64(a, b, n_);
    }

    std::uint64_t sub(std::uint64_t a, std::uint64_t b) const noexcept
    {
        return a >= b ? a - b : a + (n_ - b);
    }

    // a / 2 modulo n (n odd)
    std::uint64_t half(std::uint64_t a) const noexcept
    {
        return (a & 1u) ? (a >> 1) + (n_ >> 1) + 1 : a >> 1;
    }

    std::uint64_t mul(std::uint64_t a, std::uint64_t b) const noexcept
    {
#if defined(BOOST_MATH_HAS_INT128)
        const unsigned __int128 t {static_cast<unsigned __int128>(a) * b};
        const std::uint64_t m {static_cast<std::uint64_t>(t) * n_inv_neg_};
        const unsigned __int128 mn {static_cast<unsigned __int128>(m) * n_};
        const std::uint64_t t_hi {static_cast<std::uint64_t>(t >> 64)};
        const std::uint64_t mn_hi {static_cast<std::uint64_t>(mn >> 64)};
        return reduce(t_hi, mn_hi, static_cast<std::uint64_t>(t) != 0);
#elif defined(_MSC_VER) && !defined(__clang__) && defined(_M_X64)
        std::uint64_t t_hi {};
        const std::uint64_t t_lo {_umul128(a, b, &t_hi)};
        const std::uint64_t m {t_lo * n_inv_neg_};
        std::uint64_t mn_hi {};
        _umul128(m, n_, &mn_hi);
        return reduce(t_hi, mn_hi, t_lo != 0);
#else
        return mulmod_u64(a, b, n_);
#endif
    }

    std::uint64_t pow(std::uint64_t base, std::uint64_t exp) const noexcept
    {
        std::uint64_t result {one_};
        while (exp != 0)
        {
            if (exp & 1u)
            {
                result = mul(result, base);
            }
            base = mul(base, base);
            exp >>= 1;
        }
        return result;
    }

private:
    // The low halves of t and m * n sum to a multiple of 2^64, so only the carry matters.
    std::uint64_t reduce(std::uint64_t t_hi, std::uint64_t mn_hi, bool carry) const noexcept
    {
        const std::uint64_t r {t_hi + mn_hi + (carry ? 1u : 0u)};
        if (r < t_hi || r >= n_)
        {
            return r - n_;
        }
        return r;
    }

    std::uint64_t n_;
    std::uint64_t n_inv_neg_ {0};
    std::uint64_t one_;
    std::uint64_t r2_;
};

// Deterministic primality for all 64-bit values: trial division, then Baillie-PSW (a strong
// base-2 test and a strong Lucas test with Selfridge's parameters), which has been verified
// to have no counterexample below 2^64.
inline bool is_prime_u64(std::uint64_t n) noexcept
{
    if (n < 2)
    {
        return false;
    }
    for (std::size_t i {0}; i < 16; ++i)
    {
        const std::uint64_t p {small_primes_to_353[i]};
        if (n % p == 0)
        {
            return n == p;
        }
    }
    if (n < 59u * 59u)
    {
        return true;
    }

    const montgomery_u64 mont {n};
    const std::uint64_t one {mont.one()};
    const std::uint64_t minus_one {n - one};

    // strong Fermat test to base 2
    const int r {boost::math::tools::countr_zero(n - 1)};
    const std::uint64_t d {(n - 1) >> r};
    std::uint64_t x {mont.pow(mont.to_form(2), d)};
    if (x != one && x != minus_one)
    {
        bool witness {true};
        for (int i {1}; i < r; ++i)
        {
            x = mont.mul(x, x);
            if (x == minus_one)
            {
                witness = false;
                break;
            }
        }
        if (witness)
        {
            return false;
        }
    }

    // perfect squares would make the search for D below loop forever
    const std::uint64_t root {boost::math::tools::isqrt(n)};
    if (root * root == n)
    {
        return false;
    }

    // Selfridge method A: D = 5, -7, 9, -11, ... with (D / n) == -1; P = 1, Q = (1 - D) / 4
    long d_small {5};
    while (true)
    {
        const std::uint64_t magnitude {static_cast<std::uint64_t>(d_small < 0 ? -d_small : d_small)};
        int j {jacobi_u64(magnitude % n, n)};
        if (d_small < 0 && n % 4 == 3)
        {
            j = -j;
        }
        if (j == -1)
        {
            break;
        }
        if (j == 0 && magnitude % n != 0)
        {
            return false;
        }
        d_small = d_small > 0 ? -(d_small + 2) : -(d_small - 2);
    }
    const long q_small {(1 - d_small) / 4};
    const std::uint64_t d_res {d_small < 0 ? n - (static_cast<std::uint64_t>(-d_small) % n) : static_cast<std::uint64_t>(d_small) % n};
    const std::uint64_t q_res {q_small < 0 ? n - (static_cast<std::uint64_t>(-q_small) % n) : static_cast<std::uint64_t>(q_small) % n};
    const std::uint64_t big_d {mont.to_form(d_res)};
    const std::uint64_t big_q {mont.to_form(q_res)};

    const int s {boost::math::tools::countr_zero(n + 1)};
    const std::uint64_t k {(n + 1) >> s};
    std::uint64_t u {one};
    std::uint64_t v {one};
    std::uint64_t qk {big_q};
    for (int i {63 - boost::math::tools::countl_zero(k)}; i-- > 0;)
    {
        u = mont.mul(u, v);
        v = mont.sub(mont.mul(v, v), mont.add(qk, qk));
        qk = mont.mul(qk, qk);
        if (((k >> i) & 1u) != 0)
        {
            const std::uint64_t u_next {mont.half(mont.add(u, v))};
            const std::uint64_t v_next {mont.half(mont.add(mont.mul(big_d, u), v))};
            u = u_next;
            v = v_next;
            qk = mont.mul(qk, big_q);
        }
    }
    if (u == 0 || v == 0)
    {
        return true;
    }
    for (int i {1}; i < s; ++i)
    {
        v = mont.sub(mont.mul(v, v), mont.add(qk, qk));
        if (v == 0)
        {
            return true;
        }
        qk = mont.mul(qk, qk);
    }
    return false;
}

// Reduces x into [0, m) for signed Integer types.
template <class Integer>
inline Integer mod_positive(Integer x, const Integer& m)
{
    x %= m;
    if (x < 0)
    {
        x += m;
    }
    return x;
}

// base^exp mod m by square and multiply.
template <class Integer>
Integer powm(Integer base, Integer exp, const Integer& m)
{
    Integer result {1};
    base = mod_positive(base, m);
    while (exp > 0)
    {
        if (exp % 2 == 1)
        {
            result = (result * base) % m;
        }
        exp /= 2;
        if (exp > 0)
        {
            base = (base * base) % m;
        }
    }
    return result;
}

// Jacobi symbol (a / n) for a >= 0 and odd n > 0.
template <class Integer>
int jacobi(Integer a, Integer n)
{
    a %= n;
    int result {1};
    while (a != 0)
    {
        while (a % 2 == 0)
        {
            a /= 2;
            const Integer r {n % 8};
            if (r == 3 || r == 5)
            {
                result = -result;
            }
        }
        std::swap(a, n);
        if (a % 4 == 3 && n % 4 == 3)
        {
            result = -result;
        }
        a %= n;
    }
    return n == 1 ? result : 0;
}

// Jacobi symbol for a small signed numerator.
template <class Integer>
int jacobi_signed(long a, const Integer& n)
{
    if (a < 0)
    {
        const int sign {n % 4 == 3 ? -1 : 1};
        return sign * jacobi(Integer(-a), n);
    }
    return jacobi(Integer(a), n);
}

// Number of significant bits in n >= 0.
template <class Integer>
unsigned bit_length(Integer n)
{
    unsigned bits {0};
    while (n > 0)
    {
        n /= 2;
        ++bits;
    }
    return bits;
}

// Largest r with r * r <= n.
template <class Integer>
Integer isqrt(const Integer& n)
{
    if (n < 2)
    {
        return n;
    }
    Integer x {1};
    x <<= (bit_length(n) + 1) / 2;
    while (true)
    {
        const Integer y {(x + n / x) / 2};
        if (y >= x)
        {
            return x;
        }
        x = y;
    }
}

// Largest r with r^k <= n, for k >= 2.
template <class Integer>
Integer iroot(const Integer& n, unsigned k)
{
    if (k == 2)
    {
        return isqrt(n);
    }
    Integer low {1};
    Integer high {1};
    high <<= (bit_length(n) / k + 1);
    while (low < high)
    {
        const Integer mid {(low + high + 1) / 2};
        Integer power {1};
        bool too_big {false};
        for (unsigned i {0}; i < k; ++i)
        {
            power *= mid;
            if (power > n)
            {
                too_big = true;
                break;
            }
        }
        if (too_big)
        {
            high = mid - 1;
        }
        else
        {
            low = mid;
        }
    }
    return low;
}

// True when n == r^k for some r and some prime k <= max_k.
template <class Integer>
bool is_perfect_power(const Integer& n, unsigned max_k)
{
    for (std::size_t i {0}; i < 71 && small_primes_to_353[i] <= max_k; ++i)
    {
        const unsigned k {small_primes_to_353[i]};
        const Integer r {iroot(n, k)};
        Integer power {1};
        for (unsigned j {0}; j < k; ++j)
        {
            power *= r;
        }
        if (power == n)
        {
            return true;
        }
    }
    return false;
}

// Strong Fermat probable prime test to base 2 for odd n > 2.
template <class Integer>
bool strong_fermat_base2(const Integer& n)
{
    Integer d {n - 1};
    unsigned r {0};
    while (d % 2 == 0)
    {
        d /= 2;
        ++r;
    }
    Integer x {powm(Integer(2), d, n)};
    if (x == 1 || x == n - 1)
    {
        return true;
    }
    for (unsigned i {1}; i < r; ++i)
    {
        x = (x * x) % n;
        if (x == n - 1)
        {
            return true;
        }
    }
    return false;
}

// Strong Lucas probable prime test with Selfridge's method A parameters for odd n > 2.
// n must not be a perfect square (the D search would not terminate).
template <class Integer>
bool strong_lucas_selfridge(const Integer& n)
{
    // Find D = 5, -7, 9, -11, ... with (D / n) == -1
    long d_small {5};
    while (true)
    {
        const int j {jacobi_signed(d_small, n)};
        if (j == -1)
        {
            break;
        }
        if (j == 0)
        {
            // gcd(|D|, n) > 1: composite unless n is |D| itself
            return n == Integer(d_small < 0 ? -d_small : d_small);
        }
        d_small = d_small > 0 ? -(d_small + 2) : -(d_small - 2);
    }

    const Integer big_d {static_cast<Integer>(d_small)};
    const Integer big_q {mod_positive(Integer((1 - d_small) / 4), n)};

    Integer d {n + 1};
    unsigned s {0};
    while (d % 2 == 0)
    {
        d /= 2;
        ++s;
    }

    // Binary ladder for U_d, V_d and Q^d with P == 1
    Integer u {1};
    Integer v {1};
    Integer qk {big_q};
    const unsigned bits {bit_length(d)};
    for (unsigned i {bits - 1}; i-- > 0;)
    {
        // k -> 2k
        u = (u * v) % n;
        v = mod_positive(Integer(v * v - 2 * qk), n);
        qk = (qk * qk) % n;
        if (((d >> i) & 1) == 1)
        {
            // k -> k + 1
            Integer u_next {u + v};
            if (u_next % 2 != 0)
            {
                u_next += n;
            }
            u_next /= 2;
            Integer v_next {big_d * u + v};
            v_next = mod_positive(v_next, n);
            if (v_next % 2 != 0)
            {
                v_next += n;
            }
            v_next /= 2;
            u = u_next % n;
            v = v_next % n;
            qk = (qk * big_q) % n;
        }
    }

    if (u == 0 || v == 0)
    {
        return true;
    }
    for (unsigned i {1}; i < s; ++i)
    {
        v = mod_positive(Integer(v * v - 2 * qk), n);
        if (v == 0)
        {
            return true;
        }
        qk = (qk * qk) % n;
    }
    return false;
}

// Baillie-PSW probable prime test: no counterexample is known, deterministic below 2^64.
template <class Integer>
bool is_probable_prime_bpsw(const Integer& n)
{
    if (n < 2)
    {
        return false;
    }
    for (std::size_t i {0}; i < 71; ++i)
    {
        const Integer p {small_primes_to_353[i]};
        if (n % p == 0)
        {
            return n == p;
        }
    }
    if (n < Integer(359u * 359u))
    {
        return true;
    }
    if (!strong_fermat_base2(n))
    {
        return false;
    }
    const Integer r {isqrt(n)};
    if (r * r == n)
    {
        return false;
    }
    return strong_lucas_selfridge(n);
}

// Pseudosquares L_p: the least non-square m == 1 (mod 8) that is a quadratic residue modulo
// every odd prime q <= p. Each row stores the smallest prime attaining a distinct value, so
// choosing the first row with L_p > n / s uses the fewest bases. Values from Hugh Williams,
// "Edouard Lucas and Primality Testing", table 16.3.1, extended by Wooding (2003).
struct pseudosquare_entry
{
    std::uint16_t p;
    std::uint64_t high;
    std::uint64_t low;
};

inline constexpr pseudosquare_entry pseudosquares[49] =
{
    {  3, 0, 73ULL },                     {  5, 0, 241ULL },                    {  7, 0, 1009ULL },
    { 11, 0, 2641ULL },                   { 13, 0, 8089ULL },                   { 17, 0, 18001ULL },
    { 19, 0, 53881ULL },                  { 23, 0, 87481ULL },                  { 29, 0, 117049ULL },
    { 31, 0, 515761ULL },                 { 37, 0, 1083289ULL },                { 41, 0, 3206641ULL },
    { 43, 0, 3818929ULL },                { 47, 0, 9257329ULL },                { 53, 0, 22000801ULL },
    { 59, 0, 48473881ULL },               { 67, 0, 175244281ULL },              { 71, 0, 427733329ULL },
    { 79, 0, 898716289ULL },              { 83, 0, 2805544681ULL },             { 101, 0, 10310263441ULL },
    { 103, 0, 23616331489ULL },           { 107, 0, 85157610409ULL },           { 113, 0, 196265095009ULL },
    { 131, 0, 2871842842801ULL },         { 149, 0, 26250887023729ULL },        { 157, 0, 112434732901969ULL },
    { 173, 0, 178936222537081ULL },       { 181, 0, 696161110209049ULL },       { 193, 0, 2854909648103881ULL },
    { 197, 0, 6450045516630769ULL },      { 211, 0, 11641399247947921ULL },     { 227, 0, 190621428905186449ULL },
    { 229, 0, 196640148121928601ULL },    { 233, 0, 712624335095093521ULL },    { 239, 0, 1773855791877850321ULL },
    { 241, 0, 2327687064124474441ULL },   { 251, 0, 6384991873059836689ULL },   { 257, 0, 8019204661305419761ULL },
    { 263, 0, 10198100582046287689ULL },
    { 277, 0x3ULL, 0xc956f827e0524359ULL },       // 69848288320900186969
    { 281, 0xbULL, 0x539315b3b1268d59ULL },       // 208936365799044975961
    { 283, 0x1cULL, 0xec87d86ca60b50a1ULL },      // 533552663339828203681
    { 293, 0x32ULL, 0xc6d3496f20db3d81ULL },      // 936664079266714697089
    { 311, 0x74ULL, 0x210967a12ba94be1ULL },      // 2142202860370269916129
    { 331, 0x2e3ULL, 0xec11ddc09fd65c51ULL },     // 13649154491558298803281
    { 337, 0x753ULL, 0x641c14b397c27bf1ULL },     // 34594858801670127778801
    { 347, 0x1511ULL, 0x85fdf38d1fc9ce21ULL },    // 99492945930479213334049
    { 353, 0x3e8bULL, 0xaba417e222ca5091ULL }     // 295363187400900310880401
};

// Builds the pseudosquare value of a table row in the Integer type.
template <class Integer>
Integer pseudosquare_value(std::size_t index)
{
    Integer value {pseudosquares[index].high};
    value <<= 32;
    value <<= 32;
    value += Integer(pseudosquares[index].low);
    return value;
}

// Smallest table index whose pseudosquare exceeds n_over_s, or 49 when the table is exhausted.
template <class Integer>
std::size_t pseudosquare_index_for(const Integer& n_over_s)
{
    for (std::size_t i {0}; i < 49; ++i)
    {
        if (pseudosquare_value<Integer>(i) > n_over_s)
        {
            return i;
        }
    }
    return 49;
}

// Lukes-Patterson-Williams test. Preconditions: n odd, n has no prime factor <= s,
// s >= 2, and pseudosquares[index] satisfies L_p > n / s. Returns true iff n is prime.
// Conditions (Sorenson, "The pseudosquares prime sieve", theorem 2.1):
//   (1) q^((n-1)/2) == +-1 (mod n) for every prime q <= p
//   (2) 2^((n-1)/2) == -1 (mod n) when n == 5 (mod 8)
//   (3) some q gives -1 when n == 1 (mod 8)
// which prove n is a prime or a prime power; the perfect power check removes the latter.
template <class Integer>
bool pseudosquares_prime_test(const Integer& n, std::size_t index, std::uint64_t s)
{
    const unsigned p {pseudosquares[index].p};
    const Integer exponent {(n - 1) / 2};
    const Integer n_mod_8 {n % 8};
    bool saw_minus_one {false};

    for (std::size_t i {0}; i < 71 && small_primes_to_353[i] <= p; ++i)
    {
        const Integer q {small_primes_to_353[i]};
        const Integer e {powm(q, exponent, n)};
        if (e == n - 1)
        {
            saw_minus_one = true;
        }
        else if (e != 1)
        {
            return false;
        }
        else if (i == 0 && n_mod_8 == 5)
        {
            return false;
        }
    }
    if (n_mod_8 == 1 && !saw_minus_one)
    {
        return false;
    }

    // Any prime power q^k with q > s has k <= log(n) / log(s)
    const unsigned log2_s {s < 2 ? 1u : static_cast<unsigned>(63 - boost::math::tools::countl_zero(s))};
    const unsigned max_k {bit_length(n) / log2_s};
    if (max_k >= 2 && is_perfect_power(n, max_k))
    {
        return false;
    }
    return true;
}

} // namespace boost::math::detail::prime_sieve

#endif // BOOST_MATH_HAS_NVRTC
#endif // BOOST_MATH_SF_DETAIL_PRIME_SIEVE_PRIMALITY_HPP
