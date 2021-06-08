#include <boost/math/fft.hpp>

#include <complex>
#include <vector>
#include <array>

template<class T,int N>
void transform()
{   
    // test same type of iterator
    std::vector<T> A(N),B(A.size());
    boost::math::fft::dft_forward(A.begin(),A.end(),B.begin());
    boost::math::fft::dft_backward(A.begin(),A.end(),B.begin());
    
    // test with raw pointers
    boost::math::fft::dft_forward(A.data(),A.data()+A.size(),B.data());
    boost::math::fft::dft_backward(A.data(),A.data()+A.size(),B.data());

    const auto & cA = A;
    // const iterator as input
    boost::math::fft::dft_forward(cA.begin(),cA.end(),B.begin());
    boost::math::fft::dft_backward(cA.begin(),cA.end(),B.begin());
    
    // const pointer as input
    boost::math::fft::dft_forward(cA.data(),cA.data()+cA.size(),B.data());
    boost::math::fft::dft_backward(cA.data(),cA.data()+cA.size(),B.data());
    
    std::array<T,N> C;
    // input as vector::iterator, output as array::iterator
    boost::math::fft::dft_forward(A.begin(),A.end(),C.begin());
    boost::math::fft::dft_backward(A.begin(),A.end(),C.begin());
    boost::math::fft::dft_forward(A.data(),A.data()+A.size(),C.data());
    boost::math::fft::dft_backward(A.data(),A.data()+A.size(),C.data());
    
    // input as array::iterator, output as vector::iterator
    boost::math::fft::dft_forward(C.begin(),C.end(),B.begin());
    boost::math::fft::dft_backward(C.begin(),C.end(),B.begin());
    boost::math::fft::dft_forward(C.data(),C.data()+C.size(),B.data());
    boost::math::fft::dft_backward(C.data(),C.data()+C.size(),B.data());
}

template<class T,int N>
void back_ends()
{
    {
        boost::math::fft::fftw_dft<T> P(N);    
        std::vector<T> A(N),B(N);
        P.forward(A.data(),B.data());
        P.backward(A.data(),B.data());
    }
    {
        boost::math::fft::dft<T> P(N);    
        std::vector<T> A(N),B(N);
        P.forward(A.data(),B.data());
        P.backward(A.data(),B.data());
    }
}

int main()
{
    transform< std::complex<float>,3 >();
    transform< std::complex<double>,3 >();
    transform< std::complex<long double>,3 >();
    
    back_ends< std::complex<double>, 3  >();
    back_ends< std::complex<float>, 3  >();
    back_ends< std::complex<long double>, 3  >();
    return 0;
}
