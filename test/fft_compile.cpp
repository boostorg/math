#include <boost/math/fft.hpp>

#include <complex>
#include <vector>
#include <array>

template<class T,int N>
void compile_test()
{   
    // test same type of iterator
    std::vector<T> A(N),B(A.size());
    boost::math::fft::dft(A.begin(),A.end(),B.begin());
    
    // test with raw pointers
    boost::math::fft::dft(A.data(),A.data()+A.size(),B.data());

    const auto & cA = A;
    // const iterator as input
    boost::math::fft::dft(cA.begin(),cA.end(),B.begin());
    
    // const pointer as input
    boost::math::fft::dft(cA.data(),cA.data()+cA.size(),B.data());
    
    std::array<T,N> C;
    // input as vector::iterator, output as array::iterator
    boost::math::fft::dft(A.begin(),A.end(),C.begin());
    boost::math::fft::dft(A.data(),A.data()+A.size(),C.data());
    
    // input as array::iterator, output as vector::iterator
    boost::math::fft::dft(C.begin(),C.end(),B.begin());
    boost::math::fft::dft(C.data(),C.data()+C.size(),B.data());
    
    // fails static_assert: D and A types are different
    // std::vector<double> D(A.size());
    // boost::math::fft::dft(A.begin(),A.end(),D.begin());
}

int main()
{
    compile_test<int,3>();
    compile_test< std::complex<int>,3 >();
    compile_test< std::complex<double>,3 >();
    return 0;
}
