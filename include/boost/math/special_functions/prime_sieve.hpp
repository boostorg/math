//  (C) Copyright Matt Borland 2026.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
//  Segmented sieve of Eratosthenes with wheel factorization.
//
//    prime_sieve(policy, upper, out)          primes in [2, upper)
//    prime_range(policy, lower, upper, out)   primes in [lower, upper)
//    prime_count(policy, upper)               number of primes in [2, upper)
//    prime_count(policy, lower, upper)        number of primes in [lower, upper)
//
//  Integer may be any builtin integer or a Boost.Multiprecision integer. The policy is a
//  std::execution policy (seq / unseq sequential, par / par_unseq threaded) or
//  boost::math::execution::cuda. Requires C++17.

#ifndef BOOST_MATH_SPECIAL_FUNCTIONS_PRIME_SIEVE_HPP
#define BOOST_MATH_SPECIAL_FUNCTIONS_PRIME_SIEVE_HPP

#include <boost/math/tools/config.hpp>

#ifndef BOOST_MATH_HAS_NVRTC

#include <boost/math/special_functions/detail/prime_sieve/execution.hpp>
#include <boost/math/special_functions/detail/prime_sieve/options.hpp>
#include <boost/math/special_functions/detail/prime_sieve/integer_traits.hpp>
#include <boost/math/special_functions/detail/prime_sieve/sinks.hpp>
#include <boost/math/special_functions/detail/prime_sieve/driver.hpp>
#include <boost/math/special_functions/detail/prime_sieve/parallel.hpp>
#include <boost/math/special_functions/detail/prime_sieve/big_range.hpp>
#include <boost/math/special_functions/detail/prime_sieve/cuda.hpp>

#ifndef BOOST_MATH_BUILD_MODULE
#include <cstdint>
#include <cstddef>
#include <vector>
#include <type_traits>
#include <utility>
#include <iterator>
#include <limits>
#endif

namespace boost::math {

namespace detail::prime_sieve {

template <class T>
struct is_std_vector : std::false_type
{
};

template <class T, class Alloc>
struct is_std_vector<std::vector<T, Alloc>> : std::true_type
{
};

// Runs the 64-bit engine on [first, last] with the backend selected by Mode, feeding
// consume(const std::uint64_t*, std::size_t) in ascending order.
template <exec_mode Mode, class Consumer>
void run_mode(std::uint64_t first, std::uint64_t last, const prime_sieve_options& options, Consumer& consume)
{
    if (prefer_test_path(first, last, options))
    {
        test_range_u64(first, last, consume);
        return;
    }
#ifdef BOOST_MATH_HAS_CUDA_PRIME_SIEVE
    if constexpr (Mode == exec_mode::cuda)
    {
        cuda_range(first, last, options, consume);
        return;
    }
#endif
#ifdef BOOST_MATH_PRIME_SIEVE_HAS_THREADS
    if constexpr (Mode != exec_mode::sequential)
    {
        parallel_range(first, last, options, consume);
        return;
    }
#endif
    extract_sink<Consumer> sink {consume};
    run_u64(first, last, options, sink);
}

template <exec_mode Mode>
std::uint64_t count_mode(std::uint64_t first, std::uint64_t last, const prime_sieve_options& options)
{
    if (prefer_test_path(first, last, options))
    {
        struct counter
        {
            std::uint64_t count {0};
            void operator()(const std::uint64_t*, std::size_t n) noexcept
            {
                count += n;
            }
        } consume {};
        test_range_u64(first, last, consume);
        return consume.count;
    }
#ifdef BOOST_MATH_HAS_CUDA_PRIME_SIEVE
    if constexpr (Mode == exec_mode::cuda)
    {
        return cuda_count(first, last, options);
    }
#endif
#ifdef BOOST_MATH_PRIME_SIEVE_HAS_THREADS
    if constexpr (Mode != exec_mode::sequential)
    {
        return parallel_count(first, last, options);
    }
#endif
    return count_u64(first, last, options);
}

// Emits 2, 3, 5 as needed and narrows [lower, upper) to the engine's [first, last] (first >= 7).
// Returns false when nothing beyond the small primes remains.
template <class Integer, class Small>
bool prepare_bounds(Integer& lower, const Integer& upper, Small&& small, Integer& start, Integer& stop)
{
    lower = clamp_non_negative(lower);
    if (lower < Integer(2))
    {
        lower = Integer(2);
    }
    if (!(lower < upper))
    {
        return false;
    }
    for (const unsigned q : {2u, 3u, 5u})
    {
        if (lower <= Integer(q) && Integer(q) < upper)
        {
            small(q);
        }
    }
    start = lower < Integer(7) ? Integer(7) : lower;
    stop = upper - Integer(1);
    return !(stop < start);
}

// Appends the primes in [lower, upper) to a vector, converting from 64-bit values in batches.
template <exec_mode Mode, class Integer, class T, class Alloc>
void range_to_vector(Integer lower, Integer upper, std::vector<T, Alloc>& out, const prime_sieve_options& options)
{
    Integer start {};
    Integer stop {};
    if (!prepare_bounds(lower, upper, [&](unsigned q) { out.push_back(static_cast<T>(q)); }, start, stop))
    {
        return;
    }
    struct appender
    {
        std::vector<T, Alloc>& v;
        void operator()(const std::uint64_t* primes, std::size_t n)
        {
            v.insert(v.end(), primes, primes + n);
        }
    } consume {out};
    if (fits_u64(stop))
    {
        run_mode<Mode>(to_u64(start), to_u64(stop), options, consume);
        return;
    }
    if constexpr (!std::is_integral<Integer>::value)
    {
        if (fits_u64(start))
        {
            // Head below 2^64 through the engine (2^64 - 1 is composite), tail by the window method
            run_mode<Mode>(to_u64(start), (std::numeric_limits<std::uint64_t>::max)(), options, consume);
            start = Integer((std::numeric_limits<std::uint64_t>::max)()) + Integer(1);
        }
        big_range_impl(start, stop, options, [&](const Integer& p) { out.push_back(static_cast<T>(p)); });
    }
}

// Writes the primes in [lower, upper) to out and returns the advanced iterator.
template <exec_mode Mode, class Integer, class OutputIterator>
OutputIterator range_dispatch(Integer lower, Integer upper, OutputIterator out, const prime_sieve_options& options)
{
    Integer start {};
    Integer stop {};
    if (!prepare_bounds(lower, upper, [&](unsigned q) { *out = Integer(q); ++out; }, start, stop))
    {
        return out;
    }
    output_converter<Integer, OutputIterator> consumer {out};
    if (fits_u64(stop))
    {
        run_mode<Mode>(to_u64(start), to_u64(stop), options, consumer);
        return consumer.out;
    }
    if constexpr (std::is_integral<Integer>::value)
    {
        return consumer.out;
    }
    else
    {
        if (fits_u64(start))
        {
            run_mode<Mode>(to_u64(start), (std::numeric_limits<std::uint64_t>::max)(), options, consumer);
            start = Integer((std::numeric_limits<std::uint64_t>::max)()) + Integer(1);
        }
        return big_range(start, stop, consumer.out, options);
    }
}

// Counts the primes in [lower, upper).
template <exec_mode Mode, class Integer>
std::uint64_t count_dispatch(Integer lower, Integer upper, const prime_sieve_options& options)
{
    std::uint64_t count {0};
    Integer start {};
    Integer stop {};
    if (!prepare_bounds(lower, upper, [&](unsigned) { ++count; }, start, stop))
    {
        return count;
    }
    if (fits_u64(stop))
    {
        return count + count_mode<Mode>(to_u64(start), to_u64(stop), options);
    }
    if constexpr (std::is_integral<Integer>::value)
    {
        return count;
    }
    else
    {
        if (fits_u64(start))
        {
            count += count_mode<Mode>(to_u64(start), (std::numeric_limits<std::uint64_t>::max)(), options);
            start = Integer((std::numeric_limits<std::uint64_t>::max)()) + Integer(1);
        }
        return count + big_count(start, stop, options);
    }
}

} // namespace detail::prime_sieve

BOOST_MATH_EXPORT template <class ExecutionPolicy, class Integer, class OutputIterator,
          std::enable_if_t<detail::prime_sieve::is_execution_policy_v<std::decay_t<ExecutionPolicy>> &&
                           !detail::prime_sieve::is_std_vector<std::decay_t<OutputIterator>>::value, bool> = true>
inline OutputIterator prime_sieve(ExecutionPolicy&&, Integer upper_bound, OutputIterator out, const prime_sieve_options& options = {})
{
    constexpr auto mode = detail::prime_sieve::mode_of<std::decay_t<ExecutionPolicy>>();
    return detail::prime_sieve::range_dispatch<mode>(Integer(0), upper_bound, out, options);
}

BOOST_MATH_EXPORT template <class Integer, class OutputIterator,
          std::enable_if_t<detail::prime_sieve::is_integer_like_v<Integer> &&
                           !detail::prime_sieve::is_std_vector<std::decay_t<OutputIterator>>::value, bool> = true>
inline OutputIterator prime_sieve(Integer upper_bound, OutputIterator out, const prime_sieve_options& options = {})
{
    return detail::prime_sieve::range_dispatch<detail::prime_sieve::exec_mode::sequential>(Integer(0), upper_bound, out, options);
}

// Vector overloads append in bulk, which is noticeably faster than a back_insert_iterator.
BOOST_MATH_EXPORT template <class ExecutionPolicy, class Integer, class T, class Alloc,
          std::enable_if_t<detail::prime_sieve::is_execution_policy_v<std::decay_t<ExecutionPolicy>>, bool> = true>
inline void prime_sieve(ExecutionPolicy&&, Integer upper_bound, std::vector<T, Alloc>& out, const prime_sieve_options& options = {})
{
    constexpr auto mode = detail::prime_sieve::mode_of<std::decay_t<ExecutionPolicy>>();
    detail::prime_sieve::range_to_vector<mode>(Integer(0), upper_bound, out, options);
}

BOOST_MATH_EXPORT template <class Integer, class T, class Alloc,
          std::enable_if_t<detail::prime_sieve::is_integer_like_v<Integer>, bool> = true>
inline void prime_sieve(Integer upper_bound, std::vector<T, Alloc>& out, const prime_sieve_options& options = {})
{
    detail::prime_sieve::range_to_vector<detail::prime_sieve::exec_mode::sequential>(Integer(0), upper_bound, out, options);
}

// The two bounds may have different integer types (e.g. a std::uint64_t and a literal);
// they are converted to their common type first.
BOOST_MATH_EXPORT template <class ExecutionPolicy, class Lower, class Upper, class OutputIterator,
          std::enable_if_t<detail::prime_sieve::is_execution_policy_v<std::decay_t<ExecutionPolicy>> &&
                           detail::prime_sieve::is_integer_like_v<Lower> && detail::prime_sieve::is_integer_like_v<Upper> &&
                           !detail::prime_sieve::is_std_vector<std::decay_t<OutputIterator>>::value, bool> = true>
inline OutputIterator prime_range(ExecutionPolicy&&, Lower lower_bound, Upper upper_bound, OutputIterator out, const prime_sieve_options& options = {})
{
    using Integer = detail::prime_sieve::common_integer_t<Lower, Upper>;
    constexpr auto mode = detail::prime_sieve::mode_of<std::decay_t<ExecutionPolicy>>();
    return detail::prime_sieve::range_dispatch<mode>(detail::prime_sieve::to_common<Integer>(lower_bound), detail::prime_sieve::to_common<Integer>(upper_bound), out, options);
}

BOOST_MATH_EXPORT template <class Lower, class Upper, class OutputIterator,
          std::enable_if_t<detail::prime_sieve::is_integer_like_v<Lower> && detail::prime_sieve::is_integer_like_v<Upper> &&
                           !detail::prime_sieve::is_std_vector<std::decay_t<OutputIterator>>::value, bool> = true>
inline OutputIterator prime_range(Lower lower_bound, Upper upper_bound, OutputIterator out, const prime_sieve_options& options = {})
{
    using Integer = detail::prime_sieve::common_integer_t<Lower, Upper>;
    return detail::prime_sieve::range_dispatch<detail::prime_sieve::exec_mode::sequential>(detail::prime_sieve::to_common<Integer>(lower_bound), detail::prime_sieve::to_common<Integer>(upper_bound), out, options);
}

BOOST_MATH_EXPORT template <class ExecutionPolicy, class Lower, class Upper, class T, class Alloc,
          std::enable_if_t<detail::prime_sieve::is_execution_policy_v<std::decay_t<ExecutionPolicy>> &&
                           detail::prime_sieve::is_integer_like_v<Lower> && detail::prime_sieve::is_integer_like_v<Upper>, bool> = true>
inline void prime_range(ExecutionPolicy&&, Lower lower_bound, Upper upper_bound, std::vector<T, Alloc>& out, const prime_sieve_options& options = {})
{
    using Integer = detail::prime_sieve::common_integer_t<Lower, Upper>;
    constexpr auto mode = detail::prime_sieve::mode_of<std::decay_t<ExecutionPolicy>>();
    detail::prime_sieve::range_to_vector<mode>(detail::prime_sieve::to_common<Integer>(lower_bound), detail::prime_sieve::to_common<Integer>(upper_bound), out, options);
}

BOOST_MATH_EXPORT template <class Lower, class Upper, class T, class Alloc,
          std::enable_if_t<detail::prime_sieve::is_integer_like_v<Lower> && detail::prime_sieve::is_integer_like_v<Upper>, bool> = true>
inline void prime_range(Lower lower_bound, Upper upper_bound, std::vector<T, Alloc>& out, const prime_sieve_options& options = {})
{
    using Integer = detail::prime_sieve::common_integer_t<Lower, Upper>;
    detail::prime_sieve::range_to_vector<detail::prime_sieve::exec_mode::sequential>(detail::prime_sieve::to_common<Integer>(lower_bound), detail::prime_sieve::to_common<Integer>(upper_bound), out, options);
}

BOOST_MATH_EXPORT template <class ExecutionPolicy, class Integer,
          std::enable_if_t<detail::prime_sieve::is_execution_policy_v<std::decay_t<ExecutionPolicy>>, bool> = true>
inline std::uint64_t prime_count(ExecutionPolicy&&, Integer upper_bound, const prime_sieve_options& options = {})
{
    constexpr auto mode = detail::prime_sieve::mode_of<std::decay_t<ExecutionPolicy>>();
    return detail::prime_sieve::count_dispatch<mode>(Integer(0), upper_bound, options);
}

BOOST_MATH_EXPORT template <class ExecutionPolicy, class Lower, class Upper,
          std::enable_if_t<detail::prime_sieve::is_execution_policy_v<std::decay_t<ExecutionPolicy>> &&
                           detail::prime_sieve::is_integer_like_v<Lower> && detail::prime_sieve::is_integer_like_v<Upper>, bool> = true>
inline std::uint64_t prime_count(ExecutionPolicy&&, Lower lower_bound, Upper upper_bound, const prime_sieve_options& options = {})
{
    using Integer = detail::prime_sieve::common_integer_t<Lower, Upper>;
    constexpr auto mode = detail::prime_sieve::mode_of<std::decay_t<ExecutionPolicy>>();
    return detail::prime_sieve::count_dispatch<mode>(detail::prime_sieve::to_common<Integer>(lower_bound), detail::prime_sieve::to_common<Integer>(upper_bound), options);
}

BOOST_MATH_EXPORT template <class Integer,
          std::enable_if_t<detail::prime_sieve::is_integer_like_v<Integer>, bool> = true>
inline std::uint64_t prime_count(Integer upper_bound, const prime_sieve_options& options = {})
{
    return detail::prime_sieve::count_dispatch<detail::prime_sieve::exec_mode::sequential>(Integer(0), upper_bound, options);
}

BOOST_MATH_EXPORT template <class Lower, class Upper,
          std::enable_if_t<detail::prime_sieve::is_integer_like_v<Lower> && detail::prime_sieve::is_integer_like_v<Upper>, bool> = true>
inline std::uint64_t prime_count(Lower lower_bound, Upper upper_bound, const prime_sieve_options& options = {})
{
    using Integer = detail::prime_sieve::common_integer_t<Lower, Upper>;
    return detail::prime_sieve::count_dispatch<detail::prime_sieve::exec_mode::sequential>(detail::prime_sieve::to_common<Integer>(lower_bound), detail::prime_sieve::to_common<Integer>(upper_bound), options);
}

// Reserves room for every prime below upper_bound (Dusart's bound on the prime counting function).
BOOST_MATH_EXPORT template <class Integer, class T, class Alloc>
inline void prime_reserve(Integer upper_bound, std::vector<T, Alloc>& primes)
{
    upper_bound = detail::prime_sieve::clamp_non_negative(upper_bound);
    if (detail::prime_sieve::fits_u64(upper_bound))
    {
        primes.reserve(static_cast<std::size_t>(detail::prime_sieve::prime_count_upper_bound(detail::prime_sieve::to_u64(upper_bound))));
    }
}

} // namespace boost::math

#endif // BOOST_MATH_HAS_NVRTC
#endif // BOOST_MATH_SPECIAL_FUNCTIONS_PRIME_SIEVE_HPP
