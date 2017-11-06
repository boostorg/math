//  (C) Copyright Jeremy William Murphy 2016.

//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_TOOLS_FFT_HPP
#define BOOST_MATH_TOOLS_FFT_HPP

#include <boost/assert.hpp>
#include <boost/config.hpp>
#include <boost/core/enable_if.hpp>
#include <boost/iterator/iterator_adaptor.hpp>
#include <boost/lambda/lambda.hpp>

#include <algorithm>
#include <cmath>
#include <complex>
#include <iterator>
#include <vector>

namespace boost { namespace math { namespace tools {

/**
 * Based on code from Henry Warren Jr, Hacker's Delight, 2nd ed., 2013, p 138
 * @param x Reversed integer to increment.
 * @param m Increment.
 */
template <typename N>
void increment_reversed_integer(N &x, N m)
{
    x ^= m;
    while (x < m)
    {
        m >>= 1;
        x ^= m;
    }
}


/** 
 * Iterative FFT based on Cormen et al, Introduction to algorithms, 3rd ed., 2009
 * @tparam  I   InputIterator
 * @tparam  J   Mutable RandomAccessIterator
 * @tparam  R1  Unary Transformation f(T) -> T
 * @param   first   First element of data
 * @param   last    One after last element of data
 */
template <typename I, typename J, typename R1>
void iterative_fft_radix2(I first, I last, J A, R1 sign)
{
    typedef BOOST_DEDUCED_TYPENAME std::iterator_traits<J>::value_type output_type;
    typedef BOOST_DEDUCED_TYPENAME std::iterator_traits<I>::difference_type N;
    // Copy and permute the input.
    // TODO: Presumably can be done more elegantly with permutation_iterator.
    N const n = std::distance(first, last);
    for (N i_rev = 0; first != last; first++)
    {
        A[i_rev] = *first;
        increment_reversed_integer(i_rev, n / 2);
    }

    output_type const i(0, 1);
    N const log2_n = log2(n); // TODO: Err... C99?
    for (N s = 1; s <= log2_n; s++)
    {
        N const m = N(1) << s;
        output_type const w_m = std::exp(sign(2.0 * M_PI * i) / static_cast<double>(m));
        for (N k = 0; k < n; k += m)
        {
            output_type w = 1;
            N const m_2 = m / 2;
            for (N j = 0; j < m_2; j++)
            {
                output_type const t = w * A[k + j + m_2];
                // TODO: Assignment into arrays probably not good for performance.
                A[k + j + m_2] = A[k + j] - t;
                A[k + j] += t;
                w *= w_m;
            }
        }
    }
}


template <typename I, typename J>
void iterative_fft_radix2_forward(I first, I last, J A)
{
    typedef BOOST_DEDUCED_TYPENAME std::iterator_traits<J>::value_type complex_t;
    iterative_fft_radix2(first, last, A, std::negate<complex_t>());
}


template <typename I, typename J>
void iterative_fft_radix2_inverse(I first, I last, J A)
{
    using namespace boost::lambda;
    typedef BOOST_DEDUCED_TYPENAME std::iterator_traits<I>::value_type complex_t;

    iterative_fft_radix2(first, last, A, detail::identity<complex_t>());
    std::transform(A, A + (last - first), A, _1 / static_cast<double>(last - first));
}

}}}

#endif
