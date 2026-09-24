//  (C) Copyright Matt Borland 2026.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// Basic sanity check that the header <boost/math/special_functions/prime_sieve.hpp>
// #includes all the files that it needs to.

#include <boost/math/special_functions/prime_sieve.hpp>
#include <vector>

void compile_and_link_test()
{
    std::vector<unsigned> v;
    boost::math::prime_sieve(100u, v);
    boost::math::prime_range(10u, 100u, v);
    (void)boost::math::prime_count(100u);
    (void)boost::math::prime_count(10u, 100u);
    (void)boost::math::prime_count(boost::math::execution::cuda, 100u);
}
