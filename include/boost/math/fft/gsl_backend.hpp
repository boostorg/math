#ifndef BOOST_MATH_FFT_GSLBACKEND_HPP
#define BOOST_MATH_FFT_GSLBACKEND_HPP

#include <gsl/gsl_fft_complex.h>
#include <complex>

namespace boost { namespace math { 
namespace fft {
    
    template<class T>
    class gsl_dft;
    
    template<>
    class gsl_dft<std::complex<double> >
    {
        private:
        
        using value_type  = std::complex<double>;
        using cvalue_type = double;
        enum plan_type { forward_plan, backward_plan };
        
        std::size_t _size; 
        gsl_fft_complex_wavetable *wtable;
        gsl_fft_complex_workspace *wspace;
        
        template<class Iterator1, class Iterator2>
        void execute(plan_type p, Iterator1 in, Iterator2 out) const
        {
            using T1  = typename std::iterator_traits<Iterator1>::value_type;
            using T2 = typename std::iterator_traits<Iterator2>::value_type;
            
            static_assert(std::is_same<T1,T2>::value,
                "Input and output types mismatch");
                
            static_assert(std::is_same<value_type,T1>::value,
                "Plan and Input types mismatch");
                
            if(std::addressof(*in)!=std::addressof(*out))
            {
                // we avoid this extra step for in-place transforms
                // notice that if in==out, the following code has
                // undefined-behavior
                Iterator1 in_end{in};
                std::advance(in_end,size());
                std::copy(in,in_end,out);
            }
            
            if(p==forward_plan)
            gsl_fft_complex_forward(
                reinterpret_cast<cvalue_type*>(std::addressof(*out)),
                1, _size, wtable, wspace);
            else
            gsl_fft_complex_backward(
                reinterpret_cast<cvalue_type*>(std::addressof(*out)),
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
        
        template<class Iterator1, class Iterator2>
        void forward(Iterator1 in, Iterator2 out) const
        {
            execute(forward_plan,in,out);
        }
        template<class Iterator1, class Iterator2>
        void backward(Iterator1 in, Iterator2 out) const
        {
            execute(backward_plan,in,out);
        }
    };
    
    
} // namespace fft
} // namespace math
} // namespace boost

#endif // BOOST_MATH_FFT_GSLBACKEND_HPP


