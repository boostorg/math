// Copyright 2020, Madhur Chauhan

// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_SPECIAL_FIBO_HPP
#define BOOST_MATH_SPECIAL_FIBO_HPP

#include "boost/math/constants/constants.hpp"
#include <boost/assert.hpp>
#include <boost/math/policies/error_handling.hpp>
#include <cmath>
#include <limits>

#ifdef _MSC_VER
#pragma once
#endif

namespace boost {
namespace math {

template <typename T, class Policy>
T fibonacci(unsigned long long n, const Policy &pol) {
    using ull = unsigned long long;
    BOOST_ASSERT(n >= 0);
    if (n <= 3) return (n + 1) >> 1;
    if (std::ceil(n * std::log2(boost::math::constants::phi<long double>()) - std::log2(5) / 2.0) > std::numeric_limits<T>::digits)
        policies::raise_overflow_error<T>("boost::math::fibonacci<%1%>(unsigned long long)", "Result is too large to represent.", pol);
    ull mask = 1;
    for (int ct = 1; ct != std::numeric_limits<ull>::digits && (mask << 1) <= n; ++ct, mask <<= 1)
        ;
    T a = 1, b = 1;
    for (mask >>= 1; mask; mask >>= 1) {
        T t1 = a * a;
        a = 2 * a * b - t1, b = b * b + t1;
        if (mask & n) t1 = b, b += a, a = t1;
    }
    return a;
}
template <typename T>
T fibonacci(unsigned long long n) {
    return fibonacci<T>(n, policies::policy<>());
}
} // namespace math
} // namespace boost

#endif