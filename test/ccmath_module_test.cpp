//  (C) Copyright Matt Borland 2021.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

import boost_math_ccmath;

int main()
{
    constexpr double test = boost::math::ccmath::abs(-2.0);
    static_assert(test == 2.0);
}
