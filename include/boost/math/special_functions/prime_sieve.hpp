// Copyright 2020 Matt Borland
//
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_SPECIAL_FUNCTIONS_PRIME_SIEVE_HPP
#define BOOST_MATH_SPECIAL_FUNCTIONS_PRIME_SIEVE_HPP

#include <boost/math/special_functions/prime.hpp>
#include <boost/assert.hpp>
#include <vector>
#include <iterator>
#include <cmath>
#include <thread>
#include <memory>
#include <future>
#include <deque>
#include <numeric>
#include <iostream>
#include <algorithm>

#ifdef _MSVC_LANG
#if _MSVC_LANG >= 201703 // _MSVC_LANG == __cplusplus: https://devblogs.microsoft.com/cppblog/msvc-now-correctly-reports-__cplusplus/
#include <execution>
#endif
#else
#if __cplusplus >= 201703
#include <execution>
#endif
#endif

namespace boost { namespace math { namespace detail
{
// https://mathworld.wolfram.com/SieveofEratosthenes.html
// https://www.cs.utexas.edu/users/misra/scannedPdf.dir/linearSieve.pdf
template<class Z, class Container>
void linear_sieve(Z upper_bound, Container &c)
{
    size_t least_divisors_size{static_cast<size_t>(upper_bound + 1)};
    std::unique_ptr<Z[]> least_divisors{new Z[least_divisors_size]{0}};

    for (Z i{2}; i <= upper_bound; ++i)
    {
        if (least_divisors[i] == 0)
        {
            least_divisors[i] = i;
            c.emplace_back(i);
        }

        for (size_t j{}; j < least_divisors_size; ++j)
        {
            if (j >= c.size())
            {
                break;
            }

            else if (c[j] > least_divisors[i])
            {
                break;
            }

            else if (i * c[j] > upper_bound)
            {
                break;
            }

            else
            {
                least_divisors[i * c[j]] = c[j];
            }
        }
    }
}

//https://core.ac.uk/download/pdf/62440589.pdf
template<class Z, class PrimeContainer, class Container>
void mask_sieve(Z lower_bound, Z upper_bound, const PrimeContainer& primes, Container &c)
{
    Z limit {static_cast<Z>(std::floor(std::sqrt(static_cast<double>(upper_bound)))) + 1};

    size_t primes_size {};
    auto it{primes.begin()};
    while(it != primes.end() && *it < limit)
    {
        ++primes_size;
        ++it;
    }

    const size_t n {static_cast<size_t>(upper_bound - lower_bound + 1)};
    std::unique_ptr<bool[]> is_prime {new bool[n]};
    memset(is_prime.get(), true, sizeof(*is_prime.get()) * (n));

    #if __cplusplus >= 201703 || _MSVC_LANG >= 201703
    // Enable use of thread pool, not SIMD compatible
    std::for_each(std::execution::par, primes.begin(), it, [&is_prime, lower_bound, upper_bound](auto prime)
    {
        for(Z j {std::max(prime * prime, (lower_bound + prime - 1) / prime * prime)}; j <= upper_bound; j += prime)
        {
            is_prime[j - lower_bound] = false;
        }
    });
    #else

    for(size_t i{}; i < primes_size; ++i)
    {
        Z current_prime{primes[i]};
        for(Z j {std::max(current_prime * current_prime, (lower_bound + current_prime - 1) / current_prime * current_prime)};
            j <= upper_bound; j += current_prime)
        {
            is_prime[j - lower_bound] = false;
        }
    }
    #endif


    if(lower_bound == 1)
    {
        is_prime[0] = false;
    }

    for(Z i{lower_bound}; i <= upper_bound; ++i)
    {
        if(is_prime[i - lower_bound])
        {
            c.emplace_back(i);
        }
    }
}

// https://minds.wisconsin.edu/bitstream/handle/1793/59248/TR909.pdf?sequence=1
// Implementing S
template<class Z>
class SetS
{
private:
    std::vector<Z> srec_;
    
public:
    SetS() = default;

    constexpr explicit SetS(const Z limit)
    {
        srec_.reserve(limit);
    }

    constexpr Z next(const Z x) const noexcept
    {
        return *std::upper_bound(srec_.cbegin(), srec_.cend(), x);
    }

    constexpr auto next_it(const Z x) const noexcept
    {
        return std::upper_bound(srec_.cbegin(), srec_.cend(), x);
    }

    constexpr Z prev(const Z x) const noexcept
    {
        return *--std::lower_bound(srec_.cbegin(), srec_.cend(), x);
    }

    constexpr auto prev_it(const Z x) const noexcept
    {
        return --std::lower_bound(srec_.cbegin(), srec_.cend(), x);
    }
    
    constexpr Z max() const noexcept
    {
        return srec_.back();
    }

    void remove(const size_t current_index, const Z x) noexcept
    { 
        auto index {std::lower_bound(srec_.cbegin() + current_index + 1, srec_.cend(), x)};

        if(index != srec_.cend() && *index == x)
        {
            srec_.erase(index); 
        } 
    }

    void insert(const Z x) noexcept
    {
        srec_.emplace_back(x);
    }

    constexpr Z operator[] (const Z index) const noexcept
    {
        return srec_[index];
    }

    constexpr size_t size() const noexcept
    {
        return srec_.size();
    }
};

template<class Z>
constexpr bool is_prime(const Z n)
{
    if(n <= 1)
    {
        return false;
    }

    for(Z factor{2}; factor * factor <= n; ++factor)
    {
        if(n % factor == 0)
        {
            return false;
        }
    }

    return true;
}

//https://minds.wisconsin.edu/bitstream/handle/1793/59248/TR909.pdf?sequence=1
// 7 - A segmented Wheel Sieve [Pritchard '87]
template<class Z, class Container>
void sub_linear_wheel_sieve(Z upper_bound, Container &resultant_primes)
{
    Z limit{static_cast<Z>(std::floor(std::sqrt(static_cast<double>(upper_bound)))) + 1};
    resultant_primes.reserve(upper_bound / std::log(upper_bound));

    // Step 1 - Compute the wheel
    Z Mk {2};
    Z k {3};
    resultant_primes.emplace_back(static_cast<Z>(2));

    if(upper_bound == 2)
    {
        return;
    }

    while(Mk * k <= limit)
    {
        Mk *= k;
        resultant_primes.emplace_back(k);
        k += 2;
    }  

    // Initialze wheel wk   
    std::unique_ptr<Z[]> wk{new Z[Mk]{0}};

    for(Z i{1}; i < Mk; ++i)
    {
        if(std::gcd(i, Mk) == 1)
        {
            wk[i] = 1;
        }
    }

    // Part 3 of init wheel
    wk[Mk - 1] = static_cast<Z>(2);

    for(Z x{Mk - 2}; x > 0; --x)
    {
        if(wk[x] == 1)
        {
            Z i{x + 1};
            while(wk[i] == 0)
            {
                ++i;
            }
            wk[x] = i - x;
        }
    }

    // Step 2 - Init set S to the kth wheel extended to n
    Z d {1};
    SetS S(upper_bound);

    while(d < upper_bound)
    {
        S.insert(d);
        d += wk[d % Mk];
    }

    // Step 3 - Run the linear algorithm starting with p := next(S, 1), which is p_k+1
    // next(S, 1) = S[1]
    // 4 - A linear Algorithm
    Z p {S[1]};
    size_t p_index {1};

    while(p * p <= upper_bound)
    {
        //Remove Multiples
        Z f {p};
        size_t f_index {p_index};

        while(p * f <= upper_bound)
        {
            f = S[++f_index];
        }

        // Loop down through the values of f
        while(f >= p)
        {
            S.remove(f_index, p * f);
            f = S[--f_index];
        }

        p = S[++p_index];
    }

    // Step 4 - Write S - {1} and the first k primes to resultant_primes
    for(size_t i{1}; i < S.size(); ++i)
    {   
        resultant_primes.emplace_back(S[i]);
    }
}

// https://minds.wisconsin.edu/bitstream/handle/1793/59248/TR909.pdf?sequence=1
// 8 - A segmented Wheel Sieve [Pritchard '83]
template<class Z, class Container>
void linear_segmented_wheel_sieve(Z lower_bound, Z upper_bound, Container &resultant_primes)
{
    const Z limit{static_cast<Z>(std::floor(std::sqrt(static_cast<double>(upper_bound)))) + 1};
    const Z interval {upper_bound - lower_bound};

    // Solves small cases for benchmark, but needs better remedy
    if(lower_bound == 2 && upper_bound <= 10)
    {
        //boost::math::detail::sub_linear_wheel_sieve(upper_bound, resultant_primes);
        boost::math::detail::linear_sieve(upper_bound, resultant_primes);
        return;
    }

    else if(lower_bound == 2 && upper_bound > 10)
    {
        //boost::math::detail::sub_linear_wheel_sieve(static_cast<Z>(10), resultant_primes);
        boost::math::detail::linear_sieve(static_cast<Z>(10), resultant_primes);
        boost::math::detail::linear_segmented_wheel_sieve(static_cast<Z>(11), upper_bound, resultant_primes);
        return;
    }

    // Pre-processing
    // 1
    std::vector<Z> primes;
    boost::math::detail::linear_sieve(limit, primes); 

    Z Mk {2};
    Z k {3};
    Z k_index {1};

    while(Mk * k <= limit)
    {
        Mk *= k;
        ++k_index;
        k = primes[k_index];
    }
    
    // Initialze wheel wk
    std::vector<Z> wk;
    wk.emplace_back(static_cast<Z>(0));
    for(Z i{1}; i < Mk; ++i)
    {
        // If the number is not prime
        if(std::gcd(i, Mk) != 1)
        {
            wk.emplace_back(static_cast<Z>(0));
        }

        else
        {
            wk.emplace_back(static_cast<Z>(1));
        }
    }

    // Part 3 of init wheel
    wk.back() = static_cast<Z>(2);
    for(Z x{Mk - 2}; x > 0; --x)
    {
        if(wk[x] == 0)
        {
            continue;
        }
        
        else
        {
            Z i{x + 1};
            while(wk[i] == 0)
            {
                ++i;
            }
            wk[x] = i - x;
        }
    }

    // Pre-processing step 2
    // Done as part of step 1 for performance improvement

    // Pre-processing step 3
    std::vector<Z> factor;
    for(size_t i{static_cast<size_t>(k_index)}; i < primes.size(); ++i)
    {
        factor.emplace_back(primes[i]);
    }

    // Sieving the interval
    // Step 1
    std::unique_ptr<bool[]> mark {new bool [static_cast<size_t>(interval + 1)]};

    for(Z x {lower_bound}, i {}; x <= upper_bound; ++x, ++i)
    {
        if(std::gcd(Mk, x) == 1)
        {
            mark[i] = 1;
        }
    }

    // Step 2
    for(Z p {0}; p < static_cast<Z>(factor.size()); ++p)
    {
        Z f {factor[p]};
        Z current_prime{factor[p]};
        while(f * current_prime <= upper_bound)
        {
            if(f * current_prime > lower_bound)
            {
                mark[f * current_prime - lower_bound] = 0;
            }

            f += wk[f % Mk];
        }
        factor[p] = f;
    }

    for(Z i {}; i < interval; ++i)
    {
        if(mark[i] == 1)
        {
            resultant_primes.emplace_back(i + lower_bound);
        }
    }
}

template<class Z, class Container>
void mask_sieve(Z lower_bound, Z upper_bound, Container &c)
{
    Z limit{static_cast<Z>(std::floor(std::sqrt(static_cast<double>(upper_bound)))) + 1};
    std::vector<Z> primes {};
    primes.reserve(limit / std::log(limit));

    boost::math::detail::linear_sieve(limit, primes);

    boost::math::detail::mask_sieve(lower_bound, upper_bound, primes, c);
}

template<class Z, class PrimesContainer, class Container>
void segmented_sieve(Z lower_bound, Z upper_bound, const PrimesContainer &primes, Container &c)
{
    const Z L1_SIZE {32648};
    const Z interval {L1_SIZE * 4};
    Z current_lower_bound{lower_bound};
    Z current_upper_bound{current_lower_bound + interval};

    if(current_upper_bound > upper_bound)
    {
        current_upper_bound = upper_bound;
    }

    size_t ranges {static_cast<size_t>((upper_bound - lower_bound) / interval)};

    std::vector<std::vector<Z>> prime_vectors(ranges + 1);
    std::vector<std::future<void>> future_manager(ranges);

    Z primes_in_range {static_cast<Z>(current_upper_bound / std::log(static_cast<double>(current_upper_bound)) -
                                      current_lower_bound / std::log(static_cast<double>(current_lower_bound)))};

    for(size_t i {}; i < ranges; ++i)
    {
        prime_vectors[i].reserve(primes_in_range);

        future_manager.emplace_back(std::async([current_lower_bound, current_upper_bound, &primes, &prime_vectors, i]{
            boost::math::detail::mask_sieve(current_lower_bound, current_upper_bound, primes, prime_vectors[i]);
        }));

        current_lower_bound = current_upper_bound + 1;
        current_upper_bound += interval;
    }

    prime_vectors[ranges].reserve(primes_in_range);
    future_manager.emplace_back(std::async([current_lower_bound, upper_bound, &primes, &prime_vectors]{
        boost::math::detail::mask_sieve(current_lower_bound, upper_bound, primes, prime_vectors.back());
    }));

    for(auto &&future : future_manager)
    {
        if(future.valid())
        {
            future.get();
        }
    }

    for(auto &v : prime_vectors)
    {
        c.insert(c.end(), v.begin(), v.end());
    }
}

template<class Z, class Container>
void segmented_sieve(Z lower_bound, Z upper_bound, Container &c)
{
    Z limit{static_cast<Z>(std::floor(std::sqrt(static_cast<double>(upper_bound)))) + 1};
    std::vector<Z> primes {};
    primes.reserve(limit / std::log(limit));

    // Prepare for max value so you do not have to calculate this again
    if(limit < 4096)
    {
        boost::math::detail::linear_sieve(limit, primes);
    }

    else
    {
        boost::math::detail::mask_sieve(static_cast<Z>(2), limit, primes);
    }

    boost::math::detail::segmented_sieve(lower_bound, upper_bound, primes, c);
}
} // End namespace detail

#if __cplusplus >= 201703 || _MSVC_LANG >= 201703
template<class ExecutionPolicy, class Z, class OutputIterator>
auto prime_sieve(ExecutionPolicy&& policy, Z upper_bound, OutputIterator output) -> decltype(output)
{
    static_assert(std::is_integral<Z>::value, "No primes for floating point types");
    BOOST_ASSERT_MSG(upper_bound + 1 < std::numeric_limits<Z>::max(), "Type Overflow");

    --upper_bound; // Not inclusive, but several methods in boost::math::detail need to be

    std::vector<Z> primes {};
    primes.reserve(upper_bound / std::log(static_cast<double>(upper_bound)));

    // Range for when the linear sieve is no longer faster than threading
    if (upper_bound < 4096)
    {
        boost::math::detail::linear_sieve(upper_bound, primes);
    }

    else if (std::is_same_v<decltype(policy), std::execution::sequenced_policy>)
    {
        boost::math::detail::segmented_sieve(static_cast<Z>(2), upper_bound, primes);
    }

    else
    {
        unsigned processor_count {std::thread::hardware_concurrency()};

        // May return 0 when unable to detect
        if(processor_count == 0)
        {
            processor_count = 2;
        }

        std::vector<Z> small_primes {};
        small_primes.reserve(1000);

        // Threshold for when 2 thread performance begins to be non-linear, or when the system can only support two threads
        if(upper_bound < 16777216 || processor_count == 2)
        {
            // Split into two vectors and merge after joined to avoid data races
            std::thread t1([&small_primes] {
                boost::math::detail::linear_sieve(static_cast<Z>(8192), small_primes);
            });
            std::thread t2([upper_bound, &primes] {
                boost::math::detail::segmented_sieve(static_cast<Z>(8192), upper_bound, primes);
            });

            t1.join();
            t2.join();
            primes.insert(primes.begin(), small_primes.begin(), small_primes.end());
        }

        //If sufficiently large upper bound spawn as many threads as the system has processors for
        else
        {
            //Pre-generate all of the primes so that each thread does not have to
            Z limit{static_cast<Z>(std::floor(std::sqrt(static_cast<double>(upper_bound)))) + 1};
            std::vector<Z> pre_generated_primes {};
            pre_generated_primes.reserve(limit / std::log(limit));

            if(limit < 4096)
            {
                boost::math::detail::linear_sieve(limit, pre_generated_primes);
            }

            else
            {
                std::thread t1([&small_primes] {
                    boost::math::detail::linear_sieve(static_cast<Z>(4096), small_primes);
                });
                std::thread t2([limit, &pre_generated_primes] {
                    boost::math::detail::segmented_sieve(static_cast<Z>(4096), limit, pre_generated_primes);
                });

                t1.join();
                t2.join();
                pre_generated_primes.insert(pre_generated_primes.begin(), small_primes.begin(), small_primes.end());
            }

            std::vector<std::thread> thread_manager {};
            std::vector<std::vector<Z>> prime_vectors(processor_count);
            const Z range_per_thread = upper_bound / (processor_count);
            Z current_lower_bound {limit};
            Z current_upper_bound {current_lower_bound + range_per_thread};
            Z primes_in_range {static_cast<Z>(current_upper_bound / std::log(static_cast<double>(current_upper_bound)) -
                               current_lower_bound / std::log(static_cast<double>(current_lower_bound)))};

            for(size_t i{}; i < processor_count - 1; ++i)
            {
                prime_vectors[i].reserve(primes_in_range);

                std::thread t([current_lower_bound, current_upper_bound, &prime_vectors, i, &pre_generated_primes] {
                    boost::math::detail::segmented_sieve(current_lower_bound, current_upper_bound, pre_generated_primes,
                                                         prime_vectors[i]);
                });

                thread_manager.push_back(std::move(t));

                current_lower_bound = current_upper_bound;
                current_upper_bound += range_per_thread;
                primes_in_range = static_cast<Z>(current_upper_bound / std::log(static_cast<double>(current_upper_bound)) -
                                                 current_lower_bound / std::log(static_cast<double>(current_lower_bound)));
            }

            prime_vectors.back().reserve(primes_in_range);

            std::thread t([current_lower_bound, upper_bound, &prime_vectors, &pre_generated_primes] {
                boost::math::detail::segmented_sieve(current_lower_bound, upper_bound, pre_generated_primes, prime_vectors.back());
            });
            thread_manager.push_back(std::move(t));

            for(auto &thread : thread_manager)
            {
                thread.join();
            }

            primes.insert(primes.begin(), pre_generated_primes.begin(), pre_generated_primes.end());
            for(auto &v : prime_vectors)
            {
                primes.insert(primes.end(), v.begin(), v.end());
            }
        }
    }

    return std::move(primes.begin(), primes.end(), output);
}

template<class ExecutionPolicy, class Z, class OutputIterator>
auto prime_range(ExecutionPolicy&& policy, Z lower_bound, Z upper_bound, OutputIterator output) -> decltype(output)
{
    --upper_bound; // Not inclusive, but several methods in boost::math::detail need to be
    std::vector<Z> primes {};
    primes.reserve(upper_bound / std::log(static_cast<double>(upper_bound)));

    boost::math::prime_sieve(policy, upper_bound, std::back_inserter(primes));

    auto it{primes.begin()};
    while(*it < lower_bound && it != primes.end())
    {
        ++it;
    }

    return std::move(it, primes.end(), output);
}
#endif //__cplusplus >= 201703

template<class Z, class OutputIterator>
auto prime_sieve(Z upper_bound, OutputIterator output) -> decltype(output)
{
    static_assert(std::is_integral<Z>::value, "No primes for floating point types");
    BOOST_ASSERT_MSG(upper_bound + 1 < std::numeric_limits<Z>::max(), "Type Overflow");

    --upper_bound; // Not inclusive, but several methods in boost::math::detail need to be
    std::vector<Z> primes{};
    primes.reserve(upper_bound / std::log(upper_bound));

    if (upper_bound <= 4096)
    {
        boost::math::detail::linear_sieve(upper_bound, primes);
    }

    else
    {
        boost::math::detail::segmented_sieve(static_cast<Z>(2), upper_bound, primes);
    }

    return std::move(primes.begin(), primes.end(), output);
}

template<class Z, class OutputIterator>
auto prime_range(Z lower_bound, Z upper_bound, OutputIterator output) -> decltype(output)
{
    --upper_bound; // Not inclusive, but several methods in boost::math::detail need to be
    std::vector<Z> primes{};
    primes.reserve(upper_bound / std::log(static_cast<double>(upper_bound)));

    boost::math::prime_sieve(upper_bound, std::back_inserter(primes));

    auto it{primes.begin()};
    while(*it < lower_bound && it != primes.end())
    {
        ++it;
    }

    return std::move(it, primes.end(), output);
}
}}

#endif //BOOST_MATH_SPECIAL_FUNCTIONS_PRIME_SIEVE_HPP
