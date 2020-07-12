// Copyright 2020, Madhur Chauhan

// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_SPECIAL_FIBO_HPP
#define BOOST_MATH_SPECIAL_FIBO_HPP

#include <boost/math/constants/constants.hpp>
#include <boost/math/policies/error_handling.hpp>
#include <cmath>
#include <limits>

#ifdef _MSC_VER
#pragma once
#endif

namespace boost {
namespace math {
namespace detail {
// this should be constexpr in future
const double fib_bits_phi = std::log2(boost::math::constants::phi<double>()),
             fib_bits_deno = std::log2(5.0) / 2.0;
} // namespace detail
template <typename T>
inline T fibonacci_unchecked(unsigned long long n) {
    // This function is called by the rest and computes the actual nth fibonacci number
    // First few fibonacci numbers: 0 (0th), 1 (1st), 1 (2nd), 2 (3rd), 3 (4th), 5, 8, 13, 21, 34, 55 (10th), 89
    if (n <= 2) return n == 0 ? 0 : 1;
    using ull = unsigned long long;
    ull mask = 1;
    for (int ct = 1; ct != std::numeric_limits<ull>::digits && (mask << 1) <= n; ++ct, mask <<= 1)
        ;
    T a = 1, b = 1;
    for (mask >>= 1; mask; mask >>= 1) {
        T t1 = a * a;
        a = 2 * a * b - t1, b = b * b + t1;
        if (mask & n) t1 = b, b += a, a = t1; // equivalent to: swap(a,b), b += a;
    }
    return a;
}

template <typename T, class Policy>
T inline fibonacci(unsigned long long n, const Policy &pol) {
    // check for overflow using approximation to binet's formula: F_n ~ phi^n / sqrt(5)
    if (n > 20 && n * detail::fib_bits_phi - detail::fib_bits_deno > std::numeric_limits<T>::digits)
        return policies::raise_overflow_error<T>("boost::math::fibonacci<%1%>(unsigned long long)", "Possible overflow detected.", pol);
    return fibonacci_unchecked<T>(n);
}

template <typename T>
T inline fibonacci(unsigned long long n) {
    return fibonacci<T>(n, policies::policy<>());
}
} // namespace math
} // namespace boost

#endif