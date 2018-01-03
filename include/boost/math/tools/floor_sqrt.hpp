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
//#include <boost/integer/integer_log2.hpp>

namespace boost { namespace math {

template<class Z>
Z floor_sqrt(Z N)
{
    using std::pow;
    static_assert(std::numeric_limits<Z>::is_integer,
                  "The floor_sqrt function is for taking square roots of integers.\n");


     Z x = N;
     Z y = x/2 + (x&1);
     while (y < x) {
        x = y;
        y = (x + N / x)/2;
     }
     return x;

    // Starting with x = 2^ceil(log_2(N)/2) gets rid of a few iterations. But it's delicate.
    // We need a log2 that returns a real, but that's precisely the problem we're trying to avoid.
    /*Z l = boost::integer_log2(N);
    Z exponent = l/2;
    if (pow(2, l) < N)
    {
        exponent += 1;
    }
    Z x = pow(2, exponent);
    Z y = (x + N/x)/2;
    while (y < x)
    {
        x = y;
        y = (x + N/x)/2;
    }
    return x;*/
}

}}

#endif
