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
  //#include <boost/math/fft/fftw_backend.hpp>
  //#include <boost/math/fft/gsl_backend.hpp>
  #include <boost/math/fft/bsl_backend.hpp>

  namespace boost { namespace math { namespace fft {

  // fftw_plan-like Fourier Transform API
  
  /*
    RingType axioms:
    1. Abelian group addition (operator+)
      -> closure
      -> associativity
      -> neutral element (0)
      -> inverse (operator-)
      -> commutativity
    2. Monoid multiplication (operator*)
      -> closure
      -> associativity
      -> neutral element (1)
    3. addition and multiplication compatibility
      -> left distributivity, ie. a*(b+c) == a*b + a*c
      -> right distributivity, ie. (b+c)*a == b*a + c*a
  */
  
  template<typename RingType,
           template<class U> class BackendType = bsl_dft >
  class dft : public BackendType<RingType>
  {
  public:
    using backend_t = BackendType<RingType>;

    using backend_t::forward;
    using backend_t::backward;

    constexpr dft(unsigned int n) : backend_t{ n } { }

    template<typename InputIterator,
             typename OutputIterator>
    constexpr dft(InputIterator  input_begin,
                  InputIterator  input_end,
                  OutputIterator output) : backend_t(input_begin, input_end, output){}
  };

  // std::transform-like Fourier Transform API
  template<template<class U> class backend = bsl_dft,
           typename InputIterator,
           typename OutputIterator>
  void dft_forward(InputIterator  input_begin,
                   InputIterator  input_end,
                   OutputIterator output)
  {
    using input_value_type  = typename std::iterator_traits<InputIterator >::value_type;
    using output_value_type = typename std::iterator_traits<OutputIterator>::value_type;

    static_assert(std::is_same<input_value_type, output_value_type>::value,
      "Input and output types mismatch");

    //dft<input_value_type, backend> plan(input_begin, input_end, output);
    dft<input_value_type, backend> plan(std::distance(input_begin, input_end));

    plan.forward(input_begin, output);
  }

  // std::transform-like Fourier Transform API
  template<template<class U> class backend = bsl_dft,
           class InputIterator,
           class OutputIterator>
  void dft_backward(InputIterator  input_begin,
                    InputIterator  input_end,
                    OutputIterator output)
  {
    using input_value_type  = typename std::iterator_traits<InputIterator >::value_type;
    using output_value_type = typename std::iterator_traits<OutputIterator>::value_type;

    static_assert(std::is_same<input_value_type, output_value_type>::value, 
      "Input and output types mismatch");

    //dft<input_value_type, backend> plan(input_begin, input_end, output);
    dft<input_value_type, backend> plan(std::distance(input_begin, input_end));

    plan.backward(input_begin, output);
  }

  } } } // namespace boost::math::fft

#endif // BOOST_MATH_FFT_HPP
