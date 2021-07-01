///////////////////////////////////////////////////////////////////
//  Copyright Eduardo Quintana 2021
//  Copyright Janek Kozicki 2021
//  Copyright Christopher Kormanyos 2021
//  Distributed under the Boost Software License,
//  Version 1.0. (See accompanying file LICENSE_1_0.txt
//  or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_FFT_HPP
  #define BOOST_MATH_FFT_HPP

  #include <algorithm>
  #include <iterator>
  #include <vector>
  #include <boost/math/fft/algorithms.hpp>

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
           template<class U> class BackendType>
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
  template<template<class U> class backend,
           typename InputIterator,
           typename OutputIterator>
  void dft_forward(InputIterator  input_begin,
                   InputIterator  input_end,
                   OutputIterator output)
  {
    using input_value_type  = typename std::iterator_traits<InputIterator >::value_type;
    dft<input_value_type, backend> plan(static_cast<unsigned int>(std::distance(input_begin, input_end)));
    plan.forward(input_begin, input_end, output);
  }

  // std::transform-like Fourier Transform API
  // for complex types
  template<template<class U> class backend,
           class InputIterator,
           class OutputIterator>
  void dft_backward(InputIterator  input_begin,
                    InputIterator  input_end,
                    OutputIterator output)
  {
    using input_value_type  = typename std::iterator_traits<InputIterator >::value_type;
    dft<input_value_type, backend> plan(static_cast<unsigned int>(std::distance(input_begin, input_end)));
    plan.backward(input_begin, input_end, output);
  }
  
  // std::transform-like Fourier Transform API
  // for Ring types
  template<template<class U> class backend,
           typename InputIterator,
           typename OutputIterator,
           typename value_type>
  void dft_forward(InputIterator  input_begin,
                   InputIterator  input_end,
                   OutputIterator output,
                   value_type w)
  {
    using input_value_type  = typename std::iterator_traits<InputIterator >::value_type;
    dft<input_value_type, backend> plan(static_cast<unsigned int>(std::distance(input_begin, input_end)),w);
    plan.forward(input_begin, input_end, output);
  }

  // std::transform-like Fourier Transform API
  // for Ring types
  template<template<class U> class backend,
           class InputIterator,
           class OutputIterator,
           typename value_type>
  void dft_backward(InputIterator  input_begin,
                    InputIterator  input_end,
                    OutputIterator output,
                    value_type w)
  {
    using input_value_type  = typename std::iterator_traits<InputIterator >::value_type;
    dft<input_value_type, backend> plan(static_cast<unsigned int>(std::distance(input_begin, input_end)),w);
    plan.backward(input_begin, input_end, output);
  }
  
  template<template<class U> class backend,
           typename InputIterator1,
           typename InputIterator2,
           typename OutputIterator>
  void convolution(InputIterator1 input1_begin,
                   InputIterator1 input1_end,
                   InputIterator2 input2_begin,
                   OutputIterator output)
  {
    using input_value_type  = typename std::iterator_traits<InputIterator1>::value_type;
    using real_value_type = typename input_value_type::value_type;
    
    const long N = std::distance(input1_begin,input1_end);
    const long N_extended = detail::is_power2(N) ? N : detail::upper_bound_power2(2*N-1);
    
    std::vector<input_value_type> In1(N_extended),In2(N_extended),Out(N_extended);
    
    std::copy(input1_begin,input1_end,In1.begin());
    
    InputIterator2 input2_end{input2_begin};
    std::advance(input2_end,N);
    std::copy(input2_begin,input2_end,In2.begin());
    
    // padding
    for(long i=N;i<N_extended;++i)
      In1[i]=In2[i]=input_value_type{0};
    
    // fake N-periodicity
    if(N!=N_extended)
    for(long i=1;i<N;++i)
      In2[N_extended-N+i] = In2[i];
    
    dft<input_value_type, backend> plan(static_cast<unsigned int>(N_extended));
    plan.forward(In1.begin(),In1.end(),In1.begin());
    plan.forward(In2.begin(),In2.end(),In2.begin());
    
    // direct convolution
    std::transform(In1.begin(),In1.end(),In2.begin(),Out.begin(),std::multiplies<input_value_type>()); 
    
    plan.backward(Out.begin(),Out.end(),Out.begin());
    
    const real_value_type inv_N = real_value_type{1}/N_extended;
    for(auto & x : Out)
        x *= inv_N;
    
    std::copy(Out.begin(),Out.begin() + N,output);
  }
  
  } } } // namespace boost::math::fft

#endif // BOOST_MATH_FFT_HPP
