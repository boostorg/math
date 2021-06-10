#ifndef BOOST_MATH_FFT_FFTWBACKEND_HPP
#define BOOST_MATH_FFT_FFTWBACKEND_HPP

#include <fftw3.h>
#include <complex>

namespace boost { namespace math { 
namespace fft {
    
    namespace detail
    {
    template<class T>
    class fftw_traits;
    
    template<>
    class fftw_traits<double>
    {
        public:
        using value_type = std::complex<double>;
        using cvalue_type = fftw_complex;
        using plan_type  = fftw_plan;
        
        static plan_type plan_ctor(
            int n, cvalue_type* in, cvalue_type* out, int sign, unsigned int flags)
        {
            return fftw_plan_dft_1d(n,in,out,sign,flags);
        }
        static void plan_exec(
            const plan_type p, cvalue_type* in, cvalue_type* out)
        {
            fftw_execute_dft(p,in,out);
        }
        static void plan_dtor(plan_type p)
        {
            fftw_destroy_plan(p);
        }
    };
    template<>
    class fftw_traits<float>
    {
        public:
        using value_type = std::complex<float>;
        using cvalue_type = fftwf_complex;
        using plan_type  = fftwf_plan;
        
        static plan_type plan_ctor(
            int n, cvalue_type* in, cvalue_type* out, int sign, unsigned int flags)
        {
            return fftwf_plan_dft_1d(n,in,out,sign,flags);
        }
        static void plan_exec(
            const plan_type p, cvalue_type* in, cvalue_type* out)
        {
            fftwf_execute_dft(p,in,out);
        }
        static void plan_dtor(plan_type p)
        {
            fftwf_destroy_plan(p);
        }
    };
    template<>
    class fftw_traits<long double>
    {
        public:
        using value_type = std::complex<long double>;
        using cvalue_type = fftwl_complex;
        using plan_type  = fftwl_plan;
        
        static plan_type plan_ctor(
            int n, cvalue_type* in, cvalue_type* out, int sign, unsigned int flags)
        {
            return fftwl_plan_dft_1d(n,in,out,sign,flags);
        }
        static void plan_exec(
            const plan_type p, cvalue_type* in, cvalue_type* out)
        {
            fftwl_execute_dft(p,in,out);
        }
        static void plan_dtor(plan_type p)
        {
            fftwl_destroy_plan(p);
        }
    };
    } // namespace detail
    
    
    template<class T>
    class fftw_dft;
    
    template<class T>
    class fftw_dft<std::complex<T> >
    {
        private:
        
        using plan_type   = typename detail::fftw_traits<T>::plan_type;
        using cvalue_type = typename detail::fftw_traits<T>::cvalue_type;
        
        unsigned int _size; 
        plan_type forward_plan;
        plan_type backward_plan;
        
        template<class Iterator1, class Iterator2>
        void execute(const plan_type p, Iterator1 in, Iterator2 out) const
        {
            Iterator1 in_end{in};
            std::advance(in_end,size());
            std::copy(in,in_end,out);
            detail::fftw_traits<T>::plan_exec(
                p,
                reinterpret_cast<cvalue_type*>(&*out),
                reinterpret_cast<cvalue_type*>(&*out));
        }
        public:
        
        fftw_dft(unsigned int n):
            _size{n},
            forward_plan{
                detail::fftw_traits<T>::plan_ctor(_size,nullptr,nullptr,FFTW_FORWARD,
                    FFTW_ESTIMATE | FFTW_PRESERVE_INPUT)},
            backward_plan{
                detail::fftw_traits<T>::plan_ctor(_size,nullptr,nullptr,FFTW_BACKWARD,
                    FFTW_ESTIMATE | FFTW_PRESERVE_INPUT)}
        {}
        
        ~fftw_dft()
        {
            detail::fftw_traits<T>::plan_dtor(forward_plan);
            detail::fftw_traits<T>::plan_dtor(backward_plan);
        }
        unsigned int size() const {return _size;}
        
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

#endif // BOOST_MATH_FFT_FFTWBACKEND_HPP


