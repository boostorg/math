// Copyright 2020 Matt Borland
//
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <vector>
#include <type_traits>
#include <cmath>

#ifndef BOOST_MATH_SPECIAL_FUNCTIONS_PRIME_APPROXIMATION_HPP
#define BOOST_MATH_SPECIAL_FUNCTIONS_PRIME_APPROXIMATION_HPP

namespace boost::math
{
// Constexpr does not work with multiprecision types
template<typename Integer, std::enable_if_t<std::is_integral_v<Integer>, bool> = true>
constexpr Integer prime_approximation(const Integer upper_bound) noexcept
{
    constexpr auto c = 30 * ::log(113) / 113; // Magic numbers from wikipedia
    return static_cast<Integer>(::floor(c * static_cast<double>(upper_bound) / ::log(static_cast<double>(upper_bound))));
}

template<typename Integer, std::enable_if_t<std::is_integral_v<Integer>, bool> = true>
constexpr Integer prime_approximation(const Integer lower_bound, const Integer upper_bound) noexcept
{
    return prime_approximation(upper_bound) - prime_approximation(lower_bound);
}

template<typename Integer, std::enable_if_t<!std::is_integral_v<Integer>, bool> = true>
Integer prime_approximation(const Integer upper_bound) noexcept
{
    const auto c = 30 * std::log(113) / 113;
    return static_cast<Integer>(std::floor(c * static_cast<long double>(upper_bound) / std::log(static_cast<long double>(upper_bound))));
}

template<typename Integer, std::enable_if_t<!std::is_integral_v<Integer>, bool> = true>
Integer prime_approximation(const Integer lower_bound, const Integer upper_bound) noexcept
{
    return prime_approximation(upper_bound) - prime_approximation(lower_bound);
}

template<typename Integer>
inline void prime_reserve(Integer upper_bound, std::vector<Integer>& prime_container)
{
    prime_container.reserve(static_cast<std::size_t>(prime_approximation(upper_bound)));
}
} 

#endif // BOOST_MATH_SPECIAL_FUNCTIONS_PRIME_APPROXIMATION_HPP
