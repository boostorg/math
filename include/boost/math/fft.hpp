#ifndef BOOST_MATH_FFT_HPP
#define BOOST_MATH_FFT_HPP

#include <iterator>
#include <boost/math/fft/fftw_backend.hpp>

namespace boost { namespace math { 
namespace fft {
    
    // fftw_plan-like Fourier Transform API
    template<
        class T,
        template<class U> class backend = fftw_dft
        >
    class dft : public backend<T>
    {
        public:
        using backend_t = backend<T>;
        using backend_t::forward;
        using backend_t::backward;
        
        dft(unsigned int n):
            backend_t{n}
        {}
    };
    
    // std::transform-like Fourier Transform API
    template<class Iterator1, class Iterator2>
    void dft_forward(Iterator1 input_begin, Iterator1 input_end, Iterator2 output)
    {
        using T  = typename std::iterator_traits<Iterator1>::value_type;
        using T2 = typename std::iterator_traits<Iterator2>::value_type;
        
        static_assert(std::is_same<T,T2>::value,
            "Input and output types mismatch");
        
        dft<T> P(std::distance(input_begin,input_end));
        P.forward(input_begin,output);
    }
    
    // std::transform-like Fourier Transform API
    template<class Iterator1, class Iterator2>
    void dft_backward(Iterator1 input_begin, Iterator1 input_end, Iterator2 output)
    {
        using T  = typename std::iterator_traits<Iterator1>::value_type;
        using T2 = typename std::iterator_traits<Iterator2>::value_type;
        
        static_assert(std::is_same<T,T2>::value,
            "Input and output types mismatch");
        
        dft<T> P(std::distance(input_begin,input_end));
        P.backward(input_begin,output);
    }
} // namespace fft
} // namespace math
} // namespace boost
#endif // BOOST_MATH_FFT_HPP
