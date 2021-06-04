#ifndef BOOST_MATH_FFT_BFORCE_HPP
#define BOOST_MATH_FFT_BFORCE_HPP

#include <iostream>
#include <iterator>
#include <complex>

namespace boost { namespace math { 
    namespace fft
    {
        template<class T>
        struct plan
        {
           plan(int n)
           {
            std::cerr << "Fallback implementation\n";
            std::cerr << __PRETTY_FUNCTION__ << '\n';
            std::cerr << "n = " << n << "\n\n";
           }
        };
        
        template<class T>
        struct plan<std::complex<T>>
        {
            plan(int n)
            {
            std::cerr << "Special implementation complex<T>\n";
            std::cerr << __PRETTY_FUNCTION__ << '\n';
            std::cerr << "n = " << n << "\n\n";
            }
        };
        
        template<>
        struct plan<std::complex<double>>
        {
            plan(int n)
            {
            std::cerr << "Special implementation complex<double>\n";
            std::cerr << __PRETTY_FUNCTION__ << '\n';
            std::cerr << "n = " << n << "\n\n";
            }
        };
        
        template<class Iterator1, class Iterator2>
        void dft(Iterator1 beg1, Iterator1 end1, Iterator2 beg2)
        {
            using T  = typename std::iterator_traits<Iterator1>::value_type;
            using T2 = typename std::iterator_traits<Iterator2>::value_type;
            
            static_assert(std::is_same<T,T2>::value,
                "Input and output types mismatch");
            
            plan<T> P(std::distance(beg1,end1));
        }
        
    } // namespace fft
} // namespace math
} // namespace boost

#endif // BOOST_MATH_FFT_BFORCE_HPP
