// Copyright 2020 Matt Borland
//
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_SPECIAL_FUNCTIONS_DETAIL_TRIAL_DIVISION_HPP
#define BOOST_MATH_SPECIAL_FUNCTIONS_DETAIL_TRIAL_DIVISION_HPP

namespace boost{ namespace math{ namespace detail{
template<typename T>
bool is_prime(const T n) noexcept
{
    if (n <= 1)
    {
        return false;
    }

    for (T factor{2}; factor * factor <= n; ++factor)
    {
        if (n % factor == 0)
        {
            return false;
        }
    }

    return true;
}
}}}
#endif // BOOST_MATH_SPECIAL_FUNCTIONS_DETAIL_TRIAL_DIVISION_HPP
