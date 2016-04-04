//  (C) Copyright Jeremy William Murphy 2016.

//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_COMMON_FACTOR_RT_HPP
#define BOOST_MATH_COMMON_FACTOR_RT_HPP

#include <boost/assert.hpp>
#include <boost/core/enable_if.hpp>
#include <boost/type_traits.hpp>
#include <boost/mpl/and.hpp>

#include <algorithm>

namespace boost {
namespace math{
namespace detail
{
    template <typename T>
    bool odd(T const &x)
    {
        return static_cast<bool>(x & 0x1);
    }
    
    
    template <typename T>
    bool even(T const &x)
    {
        return !odd(x);
    }
    
    
    /** Stein gcd (aka 'binary gcd')
     * 
     * From Mathematics to Generic Programming, Alexander Stepanov, Daniel Rose
     */
    template <typename SteinDomain>
    SteinDomain Stein_gcd(SteinDomain m, SteinDomain n)
    {
        BOOST_ASSERT(m >= 0);
        BOOST_ASSERT(n >= 0);
        if (m == SteinDomain(0))
            return n;
        if (n == SteinDomain(0))
            return m;
        // m > 0 && n > 0
        int d_m = 0;
        while (even(m))
        {
            m >>= 1;
            ++d_m;
        }
        int d_n = 0;
        while (even(n))
        {
            n >>= 1;
            ++d_n;
        }
        // odd(m) && odd(n)
        while (m != n)
        {
            if (n > m)
                swap(n, m);
            m -= n;
            do
                m >>= 1;
            while (even(m));
        }
        // m == n
        return m << std::min(d_m, d_n);
    }

    
    /** Euclidean algorithm
     * 
     * From Mathematics to Generic Programming, Alexander Stepanov, Daniel Rose
     * 
     */
    template <typename EuclideanDomain>
    EuclideanDomain Euclid_gcd(EuclideanDomain a, EuclideanDomain b)
    {
        while (b != EuclideanDomain(0))
        {
            a = a % b;
            std::swap(a, b);
        }
        return a;
    }
}


template <typename T>
BOOST_DEDUCED_TYPENAME enable_if_c<mpl::and_< has_right_shift_assign<T>, has_less<T> >::value, T>::type
optimal_gcd(T const &a, T const &b)
{
    return detail::Stein_gcd(a, b);
}


template <typename T>
BOOST_DEDUCED_TYPENAME disable_if_c<mpl::and_< has_right_shift_assign<T>, has_less<T> >::value, T>::type
optimal_gcd(T const &a, T const &b)
{
    return detail::Euclid_gcd(a, b);
}


template <typename Integer>
BOOST_DEDUCED_TYPENAME enable_if_c<is_signed<Integer>::value, Integer>::type
gcd(Integer const &a, Integer const &b)
{
    return optimal_gcd(abs(a), abs(b));
}


template <typename Integer>
BOOST_DEDUCED_TYPENAME disable_if_c<is_signed<Integer>::value, Integer>::type
gcd(Integer const &a, Integer const &b)
{
    return optimal_gcd(a, b);
}


}}
#endif  // BOOST_MATH_COMMON_FACTOR_RT_HPP
