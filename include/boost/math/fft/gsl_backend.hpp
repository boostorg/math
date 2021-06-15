///////////////////////////////////////////////////////////////////
//  Copyright Eduardo Quintana 2021
//  Copyright Janek Kozicki 2021
//  Copyright Christopher Kormanyos 2021
//  Distributed under the Boost Software License,
//  Version 1.0. (See accompanying file LICENSE_1_0.txt
//  or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_FFT_GSLBACKEND_HPP
#define BOOST_MATH_FFT_GSLBACKEND_HPP

#include <gsl/gsl_fft_complex.h>
#include <complex>

namespace boost { namespace math { 
namespace fft {
    
    template<class T>
    class gsl_dft;
    
    template<>
    class gsl_dft< std::complex<double> >
    {
    private:
      using real_value_type    = double;
      using complex_value_type = std::complex<real_value_type>;
      enum plan_type { forward_plan, backward_plan };
      
      std::size_t _size; 
      gsl_fft_complex_wavetable *wtable;
      gsl_fft_complex_workspace *wspace;
        
      void execute(plan_type p, const complex_value_type* in, complex_value_type* out) const
      {
        if(in!=out)
        {
          // we avoid this extra step for in-place transforms
          // notice that if in==out, the following code has
          // undefined-behavior
          std::copy(in,in+size(),out);
        }
        
        if(p==forward_plan)
        gsl_fft_complex_forward(
          reinterpret_cast<real_value_type*>(std::addressof(*out)),
          1, _size, wtable, wspace);
        else
        gsl_fft_complex_backward(
          reinterpret_cast<real_value_type*>(std::addressof(*out)),
          1, _size, wtable, wspace);
      }
   public:
      
      gsl_dft(std::size_t n):
          _size{n}
      {
        wtable = gsl_fft_complex_wavetable_alloc(n);
        wspace = gsl_fft_complex_workspace_alloc(n);
      }
        
      ~gsl_dft()
      {
        gsl_fft_complex_wavetable_free(wtable);
        gsl_fft_complex_workspace_free(wspace);
      }
      std::size_t size() const {return _size;}
        
      void forward(const complex_value_type* in, complex_value_type* out) const
      {
        execute(forward_plan,in,out);
      }
      void backward(const complex_value_type* in, complex_value_type* out) const
      {
        execute(backward_plan,in,out);
      }
    };
    
    
} // namespace fft
} // namespace math
} // namespace boost

#endif // BOOST_MATH_FFT_GSLBACKEND_HPP


