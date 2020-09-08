// Copyright 2020 Matt Borland and Jonathan Sorenson
//
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_SPECIAL_FUNCTIONS_INTERVAL_SIEVE_HPP
#define BOOST_MATH_SPECIAL_FUNCTIONS_INTERVAL_SIEVE_HPP

#include <boost/dynamic_bitset.hpp>
#include <boost/multiprecision/cpp_int.hpp>
#include <boost/math/special_functions/prime_wheel.hpp>
#include <cmath>
#include <memory>
#include <array>
#include <cstdint>
#include <bitset>

namespace boost::math::detail
{
// boost::multiprecision and POD
template<class Integer>
inline double get_double(const Integer &x) noexcept
{
    return static_cast<double>(x);
}

template<class Integer, class PrimeContainer, class Container>
class IntervalSieve
{  
    
#ifdef BOOST_HAS_INT128    // Defined in GCC 4.6+, clang, intel. MSVC does not define. 
using int_128t = unsigned __int128; // One machine word smaller than the boost equivalent
#else
using int_128t = boost::multiprecision::uint128_t;
#endif

private:
    // Table of pseudo-sqares (https://mathworld.wolfram.com/Pseudosquare.html)
    // This table is from page 421, table 16.3.1, Hugh Williams' book
    // Last 8 entries added from Wooding's MS thesis, 2003, pp. 92-93
    struct pssentry
    {
        static constexpr std::size_t len {49};
        static constexpr std::array<std::int_fast16_t, len> prime
        {
            3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 67, 71, 79, 83, 101, 103, 107, 113, 131, 149, 157,
            173, 181, 193, 197, 211, 227, 229, 233, 239, 241, 251, 257, 263, 277, 281, 283, 293, 311, 331, 337, 347, 353
        };

        static constexpr std::array<int_128t, len> ps
        {
            73, 241, 1'009, 2'641, 8'089, 18'001, 53'881, 87'481, 117'049, 515'761, 1'083'289, 3'206'641, 3'818'929,
            9'257'329, 22'000'801, 48'473'881, 175'244'281, 427'733'329, 898'716'289, 2'805'544'681, 10'310'263'441,
            23'616'331'489, 85'157'610'409, 196'265'095'009, 2'871'842'842'801, 26'250'887'023'729, 112'434'732'901'969,
            178'936'222'537'081, 696'161'110'209'049, 2'854'909'648'103'881, 6'450'045'516'630'769, 11'641'399'247'947'921,
            190'621'428'905'186'449, 196'640'148'121'928'601, 712'624'335'095'093'521, 1'773'855'791'877'850'321,
            2'327'687'064'124'474'441, 6'384'991'873'059'836'689, 8'019'204'661'305'419'761, 10'198'100'582'046'287'689u,

            (static_cast<int_128t>(0x3uLL) << 64) | 0xc956f827e0524359uLL,      // 69'848'288'320'900'186'969
            (static_cast<int_128t>(0xbuLL) << 64) | 0x539315b3b1268d59uLL,      // 208'936'365'799'044'975'961
            (static_cast<int_128t>(0x1cuLL) << 64) | 0xec87d86ca60b50a1uLL,     // 533'552'663'339'828'203'681
            (static_cast<int_128t>(0x32uLL) << 64) | 0xc6d3496f20db3d81uLL,     // 936'664'079'266'714'697'089
            (static_cast<int_128t>(0x74uLL) << 64) | 0x210967a12ba94be1uLL,     // 2'142'202'860'370'269'916'129
            (static_cast<int_128t>(0x2e3uLL) << 64) | 0xec11ddc09fd65c51uLL,    // 13'649'154'491'558'298'803'281
            (static_cast<int_128t>(0x753uLL) << 64) | 0x641c14b397c27bf1uLL,    // 34'594'858'801'670'127'778'801
            (static_cast<int_128t>(0x1511uLL) << 64) | 0x85fdf38d1fc9ce21uLL,   // 99'492'945'930'479'213'334'049
            (static_cast<int_128t>(0x3e8buLL) << 64) | 0xaba417e222ca5091uLL    // 295'363'187'400'900'310'880'401
        };
    };

    static constexpr pssentry pss_{};
    static constexpr boost::math::detail::MOD210Wheel<Integer> w_{};
    std::size_t tdlimit_;

    Integer delta_;
    Integer left_;
    Integer right_;

    // https://www.researchgate.net/publication/220803585_Performance_of_C_bit-vector_implementations
    boost::dynamic_bitset<> b_;
    
    const PrimeContainer& primes_;
    std::int_fast64_t plimit_;

    void Settdlimit() noexcept;
    void SeiveLength(const Integer d) noexcept;
    void Sieve() noexcept;
    bool Psstest(const std::size_t pos) noexcept;
    void Psstestall() noexcept;
    void WriteOutput(Container &resultant_primes) noexcept;
    
public:
    IntervalSieve(const Integer left, const Integer right, const PrimeContainer &primes, Container &resultant_primes) noexcept;
    void NewRange(const Integer left, const Integer right, Container &resultant_primes) noexcept;
};

template<class Integer, class PrimeContainer, class Container>
void IntervalSieve<Integer, PrimeContainer, Container>::Settdlimit() noexcept
{
    const double dr = get_double(right_);
    const double delta = get_double(delta_);
    const double tdest = delta * std::log(dr);

    // Small cases
    if(tdest * tdest >= dr)
    {
        tdlimit_ = static_cast<std::size_t>(std::sqrt(dr));
        plimit_ = 0;
        return;
    }

    // First guess
    if(tdest <= 1ul<<30)
    {
        tdlimit_ = static_cast<std::size_t>(tdest);
    }

    else
    {
        tdlimit_ = 1ul<<30;
    }

    // Find the corresponding prime
    std::size_t i;
    for(i = pss_.len - 1; i > 0; --i)
    {
        if(static_cast<double>(pss_.ps[i]) * tdlimit_ < dr)
        {
            break;
        }
    }
    plimit_ = pss_.prime[i];

    double tdlimit_guess = 1 + std::fmod(dr, static_cast<double>(pss_.ps[i]));
    if(tdlimit_guess * tdlimit_guess >= dr)
    {
        tdlimit_ = static_cast<std::size_t>(std::sqrt(dr));
        plimit_ = 0;
    }
}

template<class Integer, class PrimeContainer, class Container>
void IntervalSieve<Integer, PrimeContainer, Container>::SeiveLength(const Integer d) noexcept
{
    Integer r {left_ % d};
    Integer start {0};

    if(r != 0)
    {
        start = d - r;
    }

    for(Integer i {start}; i >= 0 && i < b_.size(); i += d)
    {
        b_[static_cast<std::size_t>(i)] = 0;
    }
}

template<class Integer, class PrimeContainer, class Container>
void IntervalSieve<Integer, PrimeContainer, Container>::Sieve() noexcept
{
    std::int_fast64_t primes_range {};
    if(plimit_ <= 10)
    {
        primes_range = 10;
    }

    else
    {
        primes_range = plimit_;
    }

    // Sieve with pre-computed primes and then use the wheel for the remainder    
    std::size_t i {};
    for(; primes_[i] < primes_range; ++i)
    {
        SeiveLength(primes_[i]);
    }

    for(Integer j = w_.Next(primes_[--i]); j <= tdlimit_; j = w_.Next(j))
    {
        SeiveLength(j);
    }
}

template<class Integer, class PrimeContainer, class Container>
void IntervalSieve<Integer, PrimeContainer, Container>::WriteOutput(Container &resultant_primes) noexcept
{
    for(std::size_t i {}; i < b_.size(); ++i)
    {
        if(b_[i])
        {
            resultant_primes.emplace_back(left_ + i);
        }
    }
}

// Performs the pseduosqaure prime test on n = left + pos
// return 1 if prime or prime power, 0 otherwise
// Begins with a base-2 test
template<class Integer, class PrimeContainer, class Container>
bool IntervalSieve<Integer, PrimeContainer, Container>::Psstest(const std::size_t pos) noexcept
{
    const Integer n {static_cast<Integer>(left_ + pos)};
    const Integer exponent {(n - 1) / 2};
    const std::int_fast64_t nmod8 = static_cast<std::int_fast64_t>(n % 8);

    std::int_fast64_t negative_one_count {0};

    for(std::size_t i {}; i < primes_.size(); ++i)
    {
        Integer temp = primes_[i];
        temp = static_cast<Integer>(std::pow(get_double(temp), get_double(n)));

        if(temp == 1)
        {
            if(i == 0 && nmod8 == 5)
            {
                return false;
            }
        }

        else
        {
            ++temp;
            if(temp == n)
            {
                if(i > 0)
                {
                    ++negative_one_count;
                }
            }
            else
            {
                return false;
            }
        }
    }

    return (nmod8 != 1 || negative_one_count > 0);
}

template<class Integer, class PrimeContainer, class Container>
void IntervalSieve<Integer, PrimeContainer, Container>::Psstestall() noexcept
{
    for(std::size_t i {}; i < b_.size(); ++i)
    {
        if(b_[i])
        {
            if(!Psstest(i))
            {
                b_[i] = 0;
            }
        }
    }
}

template<class Integer, class PrimeContainer, class Container>
IntervalSieve<Integer, PrimeContainer, Container>::IntervalSieve(const Integer left, const Integer right, const PrimeContainer &primes, Container &resultant_primes) noexcept : 
    left_ {left}, right_ {right}, primes_ {primes}
{
    delta_ = right_ - left_;
    b_.resize(static_cast<std::size_t>(delta_), 1);
    Settdlimit();
    Sieve();
    
    if(plimit_ != 0 )
    {
        Psstestall();
    }
    
    WriteOutput(resultant_primes);
}

template<class Integer, class PrimeContainer, class Container>
void IntervalSieve<Integer, PrimeContainer, Container>::NewRange(const Integer left, const Integer right, Container &resultant_primes) noexcept
{
    left_ = left;
    right_ = right;
    delta_ = right_ - left_;

    b_.resize(static_cast<std::size_t>(delta_));
    b_.set();
    Settdlimit();
    Sieve();
    
    if(plimit_ != 0)
    {
        Psstestall();
    }
    
    WriteOutput(resultant_primes);
}
}

#endif // BOOST_MATH_SPECIAL_FUNCTIONS_INTERVAL_SIEVE_HPP
