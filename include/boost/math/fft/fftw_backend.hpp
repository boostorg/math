#ifndef BOOST_MATH_FFT_FFTWBACKEND_HPP
#define BOOST_MATH_FFT_FFTWBACKEND_HPP

#include <fftw3.h>
#include <complex>

namespace boost { namespace math { 
namespace fft {
    
    // 1-D fftw backend
    template<class T>
    class fftw_dft;
    
    
    template<>
    class fftw_dft< std::complex<double>  >
    {
        public:
        using value_type = std::complex<double>;
       
        private:
        unsigned int _size; 
        fftw_plan forward_plan;
        fftw_plan backward_plan;
        
        public:
        
        fftw_dft(unsigned int n):
            _size{n},
            forward_plan{
                fftw_plan_dft_1d(_size,nullptr,nullptr,FFTW_FORWARD,
                    FFTW_ESTIMATE | FFTW_PRESERVE_INPUT)},
            backward_plan{
                fftw_plan_dft_1d(_size,nullptr,nullptr,FFTW_BACKWARD,
                    FFTW_ESTIMATE | FFTW_PRESERVE_INPUT)}
        {}
        
        ~fftw_dft()
        {
            fftw_destroy_plan(forward_plan);
            fftw_destroy_plan(backward_plan);
        }
        unsigned int size() const {return _size;}
        
        void forward(const value_type* in, value_type* out) const
        {
            std::copy(in,in+size(),out);
            fftw_execute_dft(
                forward_plan,
                reinterpret_cast<fftw_complex*>(out),
                reinterpret_cast<fftw_complex*>(out));
        }
        void backward(const value_type* in, value_type* out) const
        {
            std::copy(in,in+size(),out);
            fftw_execute_dft(
                backward_plan,
                reinterpret_cast<fftw_complex*>(out),
                reinterpret_cast<fftw_complex*>(out));
        }
    };
    
    template<>
    class fftw_dft< std::complex<float>  >
    {
        public:
        using value_type = std::complex<float>;
       
        private:
        unsigned int _size; 
        fftwf_plan forward_plan;
        fftwf_plan backward_plan;
        
        public:
        
        fftw_dft(unsigned int n):
            _size{n},
            forward_plan{
                fftwf_plan_dft_1d(_size,nullptr,nullptr,FFTW_FORWARD,
                    FFTW_ESTIMATE | FFTW_PRESERVE_INPUT)},
            backward_plan{
                fftwf_plan_dft_1d(_size,nullptr,nullptr,FFTW_BACKWARD,
                    FFTW_ESTIMATE | FFTW_PRESERVE_INPUT)}
        {}
        
        ~fftw_dft()
        {
            fftwf_destroy_plan(forward_plan);
            fftwf_destroy_plan(backward_plan);
        }
        unsigned int size() const {return _size;}
        
        void forward(const value_type* in, value_type* out) const
        {
            std::copy(in,in+size(),out);
            fftwf_execute_dft(
                forward_plan,
                reinterpret_cast<fftwf_complex*>(out),
                reinterpret_cast<fftwf_complex*>(out));
        }
        void backward(const value_type* in, value_type* out) const
        {
            std::copy(in,in+size(),out);
            fftwf_execute_dft(
                backward_plan,
                reinterpret_cast<fftwf_complex*>(out),
                reinterpret_cast<fftwf_complex*>(out));
        }
    };
    
    // I think 'long double' does not always mean 128bit float.
    // but FFTW manual does support three types: 'float', 'double', 'long double'.
    template<>
    class fftw_dft< std::complex<long double>  >
    {
        public:
        using value_type = std::complex<long double>;
       
        private:
        unsigned int _size; 
        fftwl_plan forward_plan;
        fftwl_plan backward_plan;
        
        public:
        
        fftw_dft(unsigned int n):
            _size{n},
            forward_plan{
                fftwl_plan_dft_1d(_size,nullptr,nullptr,FFTW_FORWARD,
                    FFTW_ESTIMATE | FFTW_PRESERVE_INPUT)},
            backward_plan{
                fftwl_plan_dft_1d(_size,nullptr,nullptr,FFTW_BACKWARD,
                    FFTW_ESTIMATE | FFTW_PRESERVE_INPUT)}
        {}
        
        ~fftw_dft()
        {
            fftwl_destroy_plan(forward_plan);
            fftwl_destroy_plan(backward_plan);
        }
        unsigned int size() const {return _size;}
        
        void forward(const value_type* in, value_type* out) const
        {
            std::copy(in,in+size(),out);
            fftwl_execute_dft(
                forward_plan,
                reinterpret_cast<fftwl_complex*>(out),
                reinterpret_cast<fftwl_complex*>(out));
        }
        void backward(const value_type* in, value_type* out) const
        {
            std::copy(in,in+size(),out);
            fftwl_execute_dft(
                backward_plan,
                reinterpret_cast<fftwl_complex*>(out),
                reinterpret_cast<fftwl_complex*>(out));
        }
    };
    
} // namespace fft
} // namespace math
} // namespace boost

#endif // BOOST_MATH_FFT_FFTWBACKEND_HPP


