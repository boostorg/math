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
#include <boost/range/adaptor/strided.hpp>

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

template <typename I, typename J, typename N, typename CRU>
void inner_fft(I y_0, I y_1, J y, N n, CRU &w, CRU const &w_n)
{
    for (unsigned k = 0; k < n / 2 - 1; k++)
    {
        y[k] = y_0[k] + w * y_1[k];
        y[k + n / 2] = y_0[k] - w * y_1[k];
        w *= w_n;
    }
}

template <typename I, typename J>
void fft(I first, I last, J y)
{
    using namespace boost::math::tools::detail;
    
    typedef BOOST_DEDUCED_TYPENAME std::iterator_traits<I>::difference_type difference_type;
    typedef std::complex<double> complex_t;
    difference_type const n = distance(first, last);
    BOOST_ASSERT((n & (n - 1)) == 0); // We require that n is a power of 2.
    if (n == 1)
        return;
    complex_t const i(0, 1);
    complex_t const w_n = std::exp((2.0 * M_PI * i) / static_cast<double>(n));
    complex_t w = 1;
    std::vector<complex_t> y_0, y_1;
    y_0.resize(n / 2);
    y_1.resize(n / 2);
    
    fft(make_stride_iterator(first, 2), make_stride_iterator(last - 1, 2), y_0.begin());
    fft(make_stride_iterator(first + 1, 2), make_stride_iterator(last, 2), y_1.begin());

    for (unsigned k = 0; k < n / 2 - 1; k++)
    {
        y[k] = y_0[k] + w * y_1[k];
        y[k + n / 2] = y_0[k] - w * y_1[k];
        w *= w_n;
    }
}

}}}

#endif


