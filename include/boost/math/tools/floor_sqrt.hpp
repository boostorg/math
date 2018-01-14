/*
 *  (C) Copyright Nick Thompson 2017.
 *  Use, modification and distribution are subject to the
 *  Boost Software License, Version 1.0. (See accompanying file
 *  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
 *
 *  The integer floor_sqrt doesn't lose precision like a cast does.
 *  Based on Algorithm 5.9 of "The Joy of Factoring".
 */


#ifndef BOOST_MATH_TOOLS_FLOOR_SQRT_HPP
#define BOOST_MATH_TOOLS_FLOOR_SQRT_HPP
#include <limits>

namespace boost { namespace math {

template<class Z>
Z floor_sqrt(Z N)
{
    static_assert(std::numeric_limits<Z>::is_integer,
                  "The floor_sqrt function is for taking square roots of integers.\n");

     Z x = N;
     Z y = x/2 + (x&1);
     while (y < x) {
        x = y;
        y = (x + N / x)/2;
     }
     return x;
}
}}

#endif
