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
const double log_2 = std::log(2.0),
             fib_bits_phi = std::log(boost::math::constants::phi<double>()) / log_2,
             fib_bits_deno = std::log(5.0) / 2.0 / log_2;
} // namespace detail
template <typename T>
inline T fibonacci_unchecked(unsigned long long n) {
    // This function is called by the rest and computes the actual nth fibonacci number
    // First few fibonacci numbers: 0 (0th), 1 (1st), 1 (2nd), 2 (3rd), 3 (4th), 5, 8, 13, 21, 34, 55 (10th), 89
    if (n <= 2) return n == 0 ? 0 : 1;
    unsigned long long mask = 1;
    for (int ct = 1; ct != std::numeric_limits<unsigned long long>::digits && (mask << 1) <= n; ++ct, mask <<= 1)
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

// generator for next fibonacci number (see examples/reciprocal_fibonacci_constant.hpp)
template <typename T>
class fibonacci_next {
  public:
    // return next fibonacci number
    T operator()() {
        T ret = a;
        a = b, b += ret; // could've simply: swap(a, b), b += a;
        return ret;
    }

    // after set(nth), subsequent calls to the generator returns consecutive
    // fibonacci numbers starting with the nth fibonacci number
    void set(unsigned long long nth) {
        n = nth;
        a = fibonacci<T>(n);
        b = fibonacci<T>(n + 1);
    }

  private:
    unsigned long long n = 0;
    T a = 0, b = 1;
};

} // namespace math
} // namespace boost

#endif
