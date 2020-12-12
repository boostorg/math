// Copyright 2020 Matt Borland
//
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_SPECIAL_FUNCTIONS_DETAIL_FAST_MOD_HPP
#define BOOST_MATH_SPECIAL_FUNCTIONS_DETAIL_FAST_MOD_HPP

// https://stackoverflow.com/questions/14997165/fastest-way-to-get-a-positive-modulo-in-c-c#14997413

namespace boost::math::detail
{
template<typename T>
inline unsigned modulo(const T value, const unsigned m)
{
    int mod = static_cast<int>(value % static_cast<T>(m));
    if(mod < 0)
    {
        mod += m;
    }

    return mod;
}

// Value must be positive and m must be a power of two
template<typename T>
inline unsigned modulo_power_of_two(const T value, const unsigned m)
{
    return value & static_cast<T>(m-1);
}
}
#endif // BOOST_MATH_SPECIAL_FUNCTIONS_DETAIL_FAST_MOD_HPP
