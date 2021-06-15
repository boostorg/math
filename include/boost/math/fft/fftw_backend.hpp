///////////////////////////////////////////////////////////////////
//  Copyright Eduardo Quintana 2021
//  Copyright Janek Kozicki 2021
//  Copyright Christopher Kormanyos 2021
//  Distributed under the Boost Software License,
//  Version 1.0. (See accompanying file LICENSE_1_0.txt
//  or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_FFT_FFTWBACKEND_HPP
  #define BOOST_MATH_FFT_FFTWBACKEND_HPP
  
  #include <complex>
  #include <memory>

  #include <fftw3.h>

  namespace boost { namespace math {  namespace fft {

  namespace detail {
  
  template<typename T>
  struct fftw_traits_c_interface;

  template<>
  struct fftw_traits_c_interface<float>
  {
    using plan_type = fftwf_plan;

    using real_value_type = float;

    using complex_value_type = real_value_type[2U];

    static plan_type plan_construct(
      int n, complex_value_type* in, complex_value_type* out, int sign, unsigned int flags) 
    { 
      return ::fftwf_plan_dft_1d(n, in, out, sign, flags); 
    }

    static void plan_execute(
      plan_type plan, complex_value_type* in, complex_value_type* out) 
    { 
      ::fftwf_execute_dft(plan, in, out); 
    }

    static void plan_destroy(plan_type p) { ::fftwf_destroy_plan(p); }
  };

  template<>
  struct fftw_traits_c_interface<double>
  {
    using plan_type = fftw_plan;

    using real_value_type = double;

    using complex_value_type = real_value_type[2U];

    static plan_type plan_construct(
      int n, complex_value_type* in, complex_value_type* out, int sign, unsigned int flags) 
    { 
      return ::fftw_plan_dft_1d(n, in, out, sign, flags); 
    }

    static void plan_execute(
      plan_type plan, complex_value_type* in, complex_value_type* out) 
    { 
      ::fftw_execute_dft(plan, in, out); 
    }

    static void plan_destroy(plan_type p) { ::fftw_destroy_plan(p); }
  };

  template<>
  struct fftw_traits_c_interface<long double>
  {
    using plan_type = fftwl_plan;

    using real_value_type = long double;

    using complex_value_type = real_value_type[2U];

    static plan_type plan_construct(
      int n, complex_value_type* in, complex_value_type* out, int sign, unsigned int flags) 
    { 
      return ::fftwl_plan_dft_1d(n, in, out, sign, flags); 
    }

    static void plan_execute(
      plan_type plan, complex_value_type* in, complex_value_type* out) 
    { 
      ::fftwl_execute_dft(plan, in, out); 
    }

    static void plan_destroy(plan_type p) { ::fftwl_destroy_plan(p); }
  };

  } // namespace detail

  template<class NativeComplexType>
  class fftw_dft
  {
  private:
    using real_value_type    = typename NativeComplexType::value_type;
    using plan_type          = typename detail::fftw_traits_c_interface<real_value_type>::plan_type;
    using complex_value_type = std::complex<real_value_type>;

  public:
    constexpr fftw_dft(std::ptrdiff_t n)
      : my_size{ n },
        my_forward_plan{ 
          detail::fftw_traits_c_interface<real_value_type>::plan_construct
          (
            size(), 
            nullptr, 
            nullptr, 
            FFTW_FORWARD,  
            FFTW_ESTIMATE | FFTW_PRESERVE_INPUT
          ) },
        my_backward_plan { 
          detail::fftw_traits_c_interface<real_value_type>::plan_construct
          (
            size(), 
            nullptr, 
            nullptr, 
            FFTW_BACKWARD, 
            FFTW_ESTIMATE | FFTW_PRESERVE_INPUT
          ) }
    { }

    ~fftw_dft()
    {
      detail::fftw_traits_c_interface<real_value_type>::plan_destroy(my_forward_plan);
      detail::fftw_traits_c_interface<real_value_type>::plan_destroy(my_backward_plan);
    }

    constexpr std::ptrdiff_t size() const { return my_size; }
    
    void forward(const complex_value_type* in, complex_value_type* out) const
    {
      using local_complex_type = typename detail::fftw_traits_c_interface<real_value_type>::complex_value_type;
      
      if(in!=out)
        std::copy(in,in+size(),out);

      detail::fftw_traits_c_interface<real_value_type>::plan_execute
      (
        my_forward_plan,
        reinterpret_cast<local_complex_type*>(out),
        reinterpret_cast<local_complex_type*>(out)
      );
    }

    void backward(const complex_value_type* in, complex_value_type* out) const
    {
      using local_complex_type = typename detail::fftw_traits_c_interface<real_value_type>::complex_value_type;
      
      if(in!=out)
        std::copy(in,in+size(),out);
      
      detail::fftw_traits_c_interface<real_value_type>::plan_execute
      (
        my_backward_plan,
        reinterpret_cast<local_complex_type*>(out),
        reinterpret_cast<local_complex_type*>(out)
      );
    }

  private:
    const std::ptrdiff_t my_size;
          plan_type      my_forward_plan;
          plan_type      my_backward_plan;
  };

  } } } // namespace boost::math::fft

#endif // BOOST_MATH_FFT_FFTWBACKEND_HPP
