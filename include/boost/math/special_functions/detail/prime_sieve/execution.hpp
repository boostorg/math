//  (C) Copyright Matt Borland 2026.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
//  Execution policy plumbing: std::execution policies plus a CUDA tag.

#ifndef BOOST_MATH_SF_DETAIL_PRIME_SIEVE_EXECUTION_HPP
#define BOOST_MATH_SF_DETAIL_PRIME_SIEVE_EXECUTION_HPP

#include <boost/math/tools/config.hpp>

// Threads are available unless the library configuration disabled them
#if defined(BOOST_MATH_HAS_THREADS) && !defined(BOOST_MATH_DISABLE_THREADS)
#  define BOOST_MATH_PRIME_SIEVE_HAS_THREADS
#endif

// std::execution policies. Only the policy types are used, never the parallel algorithms,
// so no TBB is needed. Besides the library-wide check, accept libstdc++ 9 and later, MSVC,
// and libc++ built with -fexperimental-library (which omits the feature-test macro).
#if defined(BOOST_MATH_EXEC_COMPATIBLE)
#  define BOOST_MATH_PRIME_SIEVE_HAS_STD_EXECUTION
#elif defined(BOOST_MATH_PRIME_SIEVE_HAS_THREADS) && defined(__has_include)
#  if __has_include(<execution>)
#    if (defined(_GLIBCXX_RELEASE) && (_GLIBCXX_RELEASE >= 9)) || (defined(_MSC_VER) && (_MSC_VER >= 1914)) || (defined(_LIBCPP_VERSION) && defined(_LIBCPP_ENABLE_EXPERIMENTAL))
#      define BOOST_MATH_PRIME_SIEVE_HAS_STD_EXECUTION
#    endif
#  endif
#endif

#ifndef BOOST_MATH_BUILD_MODULE
#include <type_traits>
#ifdef BOOST_MATH_PRIME_SIEVE_HAS_STD_EXECUTION
#include <execution>
#endif
#endif

namespace boost::math {

namespace execution {

// Requests the CUDA backend; falls back to the CPU when the translation unit is not
// compiled with nvcc and BOOST_MATH_ENABLE_CUDA (test BOOST_MATH_HAS_CUDA_PRIME_SIEVE).
BOOST_MATH_EXPORT struct cuda_policy
{
};

BOOST_MATH_EXPORT inline constexpr cuda_policy cuda {};

} // namespace execution

namespace detail::prime_sieve {

enum class exec_mode
{
    sequential,
    parallel,
    cuda
};

template <class T>
inline constexpr bool is_cuda_policy_v = std::is_same<T, execution::cuda_policy>::value;

#ifdef BOOST_MATH_PRIME_SIEVE_HAS_STD_EXECUTION
template <class T>
inline constexpr bool is_std_policy_v = std::is_execution_policy<T>::value;
#else
template <class T>
inline constexpr bool is_std_policy_v = false;
#endif

template <class T>
inline constexpr bool is_execution_policy_v = is_cuda_policy_v<T> || is_std_policy_v<T>;

template <class T>
constexpr exec_mode mode_of() noexcept
{
    if constexpr (is_cuda_policy_v<T>)
    {
        return exec_mode::cuda;
    }
#ifdef BOOST_MATH_PRIME_SIEVE_HAS_STD_EXECUTION
    else if constexpr (std::is_same<T, std::execution::sequenced_policy>::value)
    {
        return exec_mode::sequential;
    }
#if defined(__cpp_lib_execution) && (__cpp_lib_execution >= 201902L)
    else if constexpr (std::is_same<T, std::execution::unsequenced_policy>::value)
    {
        return exec_mode::sequential;
    }
#endif
    else
    {
        return exec_mode::parallel;
    }
#else
    else
    {
        return exec_mode::sequential;
    }
#endif
}

} // namespace detail::prime_sieve
} // namespace boost::math

#endif // BOOST_MATH_SF_DETAIL_PRIME_SIEVE_EXECUTION_HPP
