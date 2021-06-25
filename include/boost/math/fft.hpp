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
           // typename Allocator = Default_Allocator 
  class dft : public BackendType<RingType>
  {
    std::vector<RingType> my_mem;
    enum class execution_type { forward, backward };
    
    template<typename InputIteratorType,
             typename OutputIteratorType>
    void execute(
      execution_type ex,
      InputIteratorType in_first, InputIteratorType in_last,
      OutputIteratorType out,
      typename std::enable_if<(   (std::is_convertible<InputIteratorType,  const RingType*>::value == true)
                               && (std::is_convertible<OutputIteratorType,       RingType*>::value == true))>::type* = nullptr)
    {
      resize(std::distance(in_first,in_last));
      
      if(ex==execution_type::backward)
        backend_t::backward(in_first,out);
      else
        backend_t::forward(in_first,out);
    }
    
    template<typename InputIteratorType,
             typename OutputIteratorType>
    void execute(
      execution_type ex,
      InputIteratorType in_first, InputIteratorType in_last,
      OutputIteratorType out,
      typename std::enable_if<(   (std::is_convertible<InputIteratorType,  const RingType*>::value == false)
                               && (std::is_convertible<OutputIteratorType,       RingType*>::value == true))>::type* = nullptr)
    {
      resize(std::distance(in_first,in_last));
      std::copy(in_first, in_last, out);
      
      if(ex==execution_type::backward)
        backend_t::backward(out,out);
      else
        backend_t::forward(out,out);
    }

    template<typename InputIteratorType,
             typename OutputIteratorType>
    void execute(
      execution_type ex,
      InputIteratorType in_first, InputIteratorType in_last,
      OutputIteratorType out,
      typename std::enable_if<(   (std::is_convertible<InputIteratorType,  const RingType*>::value == true)
                               && (std::is_convertible<OutputIteratorType,       RingType*>::value == false))>::type* = nullptr)
    {   
      resize(std::distance(in_first,in_last));
      my_mem.resize(size());
      
      if(ex==execution_type::backward)
        backend_t::backward(in_first,my_mem.data());
      else
        backend_t::forward(in_first,my_mem.data());
      
      std::copy(std::begin(my_mem), std::end(my_mem), out);
    }

    template<typename InputIteratorType,
             typename OutputIteratorType>
    void execute(
      execution_type ex,
      InputIteratorType in_first, InputIteratorType in_last,
      OutputIteratorType out,
      typename std::enable_if<(   (std::is_convertible<InputIteratorType,  const RingType*>::value == false)
                               && (std::is_convertible<OutputIteratorType,       RingType*>::value == false))>::type* = nullptr)
    {
      resize(std::distance(in_first,in_last));
      my_mem.resize(size());
      std::copy(in_first, in_last, std::begin(my_mem));
      
      if(ex==execution_type::backward)
        backend_t::backward(my_mem.data(),my_mem.data());
      else
        backend_t::forward(my_mem.data(),my_mem.data());
        
      std::copy(std::begin(my_mem),std::end(my_mem), out);
    }
    
  public:
    using backend_t = BackendType<RingType>;
    using backend_t::size;
    using backend_t::resize;

    // complex types ctor. n: the size of the dft
    constexpr dft(unsigned int n) : backend_t{ n } { }
    
    // ring types ctor. n: the size of the dft, w: an n-root of unity
    constexpr dft(unsigned int n, RingType w) : backend_t( n, w ) { }
    
    
    template<typename InputIteratorType,
             typename OutputIteratorType>
    void forward(
      InputIteratorType in_first, InputIteratorType in_last,
      OutputIteratorType out)
    {
      execute(execution_type::forward,in_first,in_last,out);
    }
    
    template<typename InputIteratorType,
             typename OutputIteratorType>
    void backward(
      InputIteratorType in_first, InputIteratorType in_last,
      OutputIteratorType out)
    {
      execute(execution_type::backward,in_first,in_last,out);
    }
  };

  // std::transform-like Fourier Transform API
  // for complex types
  template<template<class U> class backend = bsl_dft,
           typename InputIterator,
           typename OutputIterator>
  void dft_forward(InputIterator  input_begin,
                   InputIterator  input_end,
                   OutputIterator output)
  {
    using input_value_type  = typename std::iterator_traits<InputIterator >::value_type;
    dft<input_value_type, backend> plan(std::distance(input_begin, input_end));
    plan.forward(input_begin, input_end, output);
  }

  // std::transform-like Fourier Transform API
  // for complex types
  template<template<class U> class backend = bsl_dft,
           class InputIterator,
           class OutputIterator>
  void dft_backward(InputIterator  input_begin,
                    InputIterator  input_end,
                    OutputIterator output)
  {
    using input_value_type  = typename std::iterator_traits<InputIterator >::value_type;
    dft<input_value_type, backend> plan(std::distance(input_begin, input_end));
    plan.backward(input_begin, input_end, output);
  }
  
  // std::transform-like Fourier Transform API
  // for Ring types
  template<template<class U> class backend = bsl_dft,
           typename InputIterator,
           typename OutputIterator,
           typename value_type>
  void dft_forward(InputIterator  input_begin,
                   InputIterator  input_end,
                   OutputIterator output,
                   value_type w)
  {
    using input_value_type  = typename std::iterator_traits<InputIterator >::value_type;
    dft<input_value_type, backend> plan(std::distance(input_begin, input_end),w);
    plan.forward(input_begin, input_end, output);
  }

  // std::transform-like Fourier Transform API
  // for Ring types
  template<template<class U> class backend = bsl_dft,
           class InputIterator,
           class OutputIterator,
           typename value_type>
  void dft_backward(InputIterator  input_begin,
                    InputIterator  input_end,
                    OutputIterator output,
                    value_type w)
  {
    using input_value_type  = typename std::iterator_traits<InputIterator >::value_type;
    dft<input_value_type, backend> plan(std::distance(input_begin, input_end),w);
    plan.backward(input_begin, input_end, output);
  }

  } } } // namespace boost::math::fft

#endif // BOOST_MATH_FFT_HPP
