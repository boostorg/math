///////////////////////////////////////////////////////////////////
//  Copyright Eduardo Quintana 2021
//  Copyright Janek Kozicki 2021
//  Copyright Christopher Kormanyos 2021
//  Distributed under the Boost Software License,
//  Version 1.0. (See accompanying file LICENSE_1_0.txt
//  or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_FFT_HPP
  #define BOOST_MATH_FFT_HPP

  #include <iterator>
  #include <boost/math/fft/fftw_backend.hpp>
  //#include <boost/math/fft/gsl_backend.hpp>

  namespace boost { namespace math { namespace fft {

  // fftw_plan-like Fourier Transform API
  template<class T,
           template<class U> class backend>
  class dft : public backend<T>
  {
  public:
    using backend_t = backend<T>;

    using backend_t::forward;
    using backend_t::backward;

    constexpr dft(unsigned int n) : backend_t{ n } { }
  };

  // std::transform-like Fourier Transform API
  template<template<class U> class backend,
           class InputIterator,
           class OutputIterator>
  void dft_forward(InputIterator  input_begin,
                   InputIterator  input_end,
                   OutputIterator output)
  {
    using input_iterator_type  = typename std::iterator_traits<InputIterator >::value_type;
    using output_iterator_type = typename std::iterator_traits<OutputIterator>::value_type;

    static_assert(std::is_same<input_iterator_type, output_iterator_type>::value, "Input and output types mismatch");

    dft<input_iterator_type, backend> P(std::distance(input_begin,input_end));

    P.forward(input_begin, output);
  }

  // std::transform-like Fourier Transform API
  template<template<class U> class backend,
           class InputIterator,
           class OutputIterator>
  void dft_backward(InputIterator  input_begin,
                    InputIterator  input_end,
                    OutputIterator output)
  {
    using input_iterator_type  = typename std::iterator_traits<InputIterator >::value_type;
    using output_iterator_type = typename std::iterator_traits<OutputIterator>::value_type;

    static_assert(std::is_same<input_iterator_type, output_iterator_type>::value, "Input and output types mismatch");

    dft<input_iterator_type, backend> P(std::distance(input_begin,input_end));

    P.backward(input_begin, output);
  }

  } } } // namespace boost::math::fft

#endif // BOOST_MATH_FFT_HPP
