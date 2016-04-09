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

#include <cmath>
#include <complex>
#include <iterator>
#include <vector>

namespace boost { namespace math { namespace tools {
    
namespace detail
{
    template <typename Iterator>
    class stride_iterator : public iterator_adaptor<stride_iterator<Iterator>, Iterator>
    {
        typedef iterator_adaptor<stride_iterator<Iterator>, Iterator> super_t;
        typedef BOOST_DEDUCED_TYPENAME std::iterator_traits<Iterator>::difference_type base_distance;
        friend class iterator_core_access;
        
        base_distance stride;
        
    public:
        stride_iterator() {}
        
        stride_iterator(Iterator x, base_distance stride) : super_t(x), stride(stride) {}
        
        template <typename OtherIterator>
        stride_iterator(stride_iterator<OtherIterator> const &x, typename enable_if_convertible<OtherIterator, Iterator>::type* = 0) : super_t(x.base())
        {}
        
        template <typename OtherIterator>
        stride_iterator(stride_iterator<OtherIterator> const &x, base_distance stride, typename enable_if_convertible<OtherIterator, Iterator>::type* = 0) : super_t(x.base()), stride(x.stride * stride)
        {}

        // TODO: These *should* be private but compilation fails unexpectedly.
        template <class OtherDerived, class OtherIterator, class V, class C, class R, class D>
        BOOST_DEDUCED_TYPENAME super_t::difference_type
        distance_to(iterator_adaptor<OtherDerived, OtherIterator, V, C, R, D> const& y) const
        {
            return (y.base() - this->base_reference()) / stride;
        }
        
        void advance(BOOST_DEDUCED_TYPENAME super_t::difference_type n)
        {
            std::advance(this->base_reference(), stride * n);
        }
    private:
        void increment() { std::advance(this->base_reference(), stride); }
        void decrement() { std::advance(this->base_reference(), -stride); }
        
    };
    
    template <typename Iterator, typename N>
    stride_iterator<Iterator> make_stride_iterator(stride_iterator<Iterator> x, N stride)
    {
        return stride_iterator<Iterator>(x, stride);
    }
    
    template <typename Iterator, typename N>
    stride_iterator<Iterator> make_stride_iterator(Iterator x, N stride)
    {
        return stride_iterator<Iterator>(x, stride);
    }
}


// Recursive FFT based on Cormen et al, Introduction to algorithms, 3rd ed., 2009
template <typename I, typename J>
void recursive_fft(I first, I last, J y)
{
    using namespace boost::math::tools::detail;
    
    typedef BOOST_DEDUCED_TYPENAME std::iterator_traits<I>::difference_type N;
    typedef BOOST_DEDUCED_TYPENAME std::iterator_traits<J>::value_type complex_t;
    N const n = distance(first, last);
    BOOST_ASSERT((n & (n - 1)) == 0); // We require that n is a power of 2.
    if (n == 1)
    {
        *y = *first;
        return;
    }
    complex_t const i(0, 1);
    complex_t const w_n = std::exp((2.0 * M_PI * i) / static_cast<double>(n));
    complex_t w = 1;
    std::vector<complex_t> y_0, y_1;
    y_0.resize(n / 2);
    y_1.resize(n / 2);
    
    recursive_fft(make_stride_iterator(first, 2), make_stride_iterator(last, 2), y_0.begin());
    recursive_fft(make_stride_iterator(first + 1, 2), make_stride_iterator(last + 1, 2), y_1.begin());

    for (N k = 0; k < n / 2; k++)
    {
        y[k] = y_0[k] + w * y_1[k];
        y[k + n / 2] = y_0[k] - w * y_1[k];
        w *= w_n;
    }
}


// Based on code from Henry Warren Jr, Hacker's Delight, 2nd ed., 2013, p 138
template <typename N>
void increment_reversed_integer(N &x, N m)
{
    x ^= m;
    if (x < m)
    {
        do
        {
            m >>= 1;
            x ^= m;
        }
        while (x < m);
    }
}


// Iterative FFT based on Cormen et al, Introduction to algorithms, 3rd ed., 2009
template <typename I, typename J>
void iterative_fft(I first, I last, J A)
{
    typedef BOOST_DEDUCED_TYPENAME std::iterator_traits<J>::value_type complex_t;
    typedef BOOST_DEDUCED_TYPENAME std::iterator_traits<I>::difference_type N;
    // Copy and permute the input.
    // TODO: Presumably can be done more elegantly with permutation_iterator.
    N const n = std::distance(first, last);
    N i_rev = 0;
    for (; first != last; first++)
    {
        A[i_rev] = *first;
        increment_reversed_integer(i_rev, n / 2);
    }

    complex_t const i(0, 1);
    N const log2_n = log2(n); // TODO: Err... C99?
    for (N s = 1; s <= log2_n; s++)
    {
        N const m = N(1) << s;
        complex_t const w_m = std::exp((2.0 * M_PI * i) / static_cast<double>(m));
        for (N k = 0; k < n; k += m)
        {
            complex_t w = 1;
            N const m_2 = m / 2;
            for (N j = 0; j < m_2; j++)
            {
                complex_t const t = w * A[k + j + m_2];
                // TODO: Assignment into arrays probably not good for performance.
                A[k + j + m_2] = A[k + j] - t;
                A[k + j] += t;
                w *= w_m;
            }
        }
    }
}

}}}

#endif


