/*
    boost::math::fft example 01.
    
    FFT transform-like API,
    default fft engine
*/

#include <boost/math/fft.hpp>

#include <iostream>
#include <vector>
#include <complex>

namespace fft = boost::math::fft;

template<class T>
void print(const std::vector< std::complex<T> >& V)
{
    for(auto i=0UL;i<V.size();++i)
        std::cout << "V[" << i << "] = " 
            << V[i].real() << ", " << V[i].imag() << '\n';
}
int main()
{
    std::vector< std::complex<double> > A{1.0,2.0,3.0,4.0},B(A.size());
    
    // default fft engine, forward transform, out-of-place
    boost::math::fft::dft_forward(A.cbegin(),A.cend(),B.begin());
    
    print(B);
    
    // default fft engine, backward transform, in-place
    fft::dft_backward(B.cbegin(),B.cend(),B.begin());
    
    print(B);
    
    return 0;
}
