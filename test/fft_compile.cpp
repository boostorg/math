#include <boost/math/fft.hpp>
#include <boost/math/fft/fftw_backend.hpp>
#include <boost/math/fft/gsl_backend.hpp>

#include <complex>
#include <vector>
#include <array>

template<class T,int N, template< class U> class Backend >
void transform_api()
{   
    // test same type of iterator
    std::vector<T> A(N),B(A.size());
    boost::math::fft::dft_forward<Backend>(A.begin(),A.end(),B.begin());
    boost::math::fft::dft_backward<Backend>(A.begin(),A.end(),B.begin());
    
    // test with raw pointers
    boost::math::fft::dft_forward<Backend>(A.data(),A.data()+A.size(),B.data());
    boost::math::fft::dft_backward<Backend>(A.data(),A.data()+A.size(),B.data());

    const auto & cA = A;
    // const iterator as input
    boost::math::fft::dft_forward<Backend>(cA.begin(),cA.end(),B.begin());
    boost::math::fft::dft_backward<Backend>(cA.begin(),cA.end(),B.begin());
    
    // const pointer as input
    boost::math::fft::dft_forward<Backend>(cA.data(),cA.data()+cA.size(),B.data());
    boost::math::fft::dft_backward<Backend>(cA.data(),cA.data()+cA.size(),B.data());
    
    std::array<T,N> C;
    // input as vector::iterator, output as array::iterator
    boost::math::fft::dft_forward<Backend>(A.begin(),A.end(),C.begin());
    boost::math::fft::dft_backward<Backend>(A.begin(),A.end(),C.begin());
    boost::math::fft::dft_forward<Backend>(A.data(),A.data()+A.size(),C.data());
    boost::math::fft::dft_backward<Backend>(A.data(),A.data()+A.size(),C.data());
    
    // input as array::iterator, output as vector::iterator
    boost::math::fft::dft_forward<Backend>(C.begin(),C.end(),B.begin());
    boost::math::fft::dft_backward<Backend>(C.begin(),C.end(),B.begin());
    boost::math::fft::dft_forward<Backend>(C.data(),C.data()+C.size(),B.data());
    boost::math::fft::dft_backward<Backend>(C.data(),C.data()+C.size(),B.data());
}

template<class T,int N, template<class U> class Backend >
void plan_api()
{
    Backend<T> P(N);    
    std::vector<T> A(N),B(N);
    P.forward(A.data(),B.data());
    P.backward(A.data(),B.data());
}

int main()
{
  
    transform_api< std::complex<float>,      3,boost::math::fft::fftw_dft >();
#if (not defined(__clang__))
    transform_api< std::complex<double>,     3,boost::math::fft::fftw_dft >();
#endif
    transform_api< std::complex<long double>,3,boost::math::fft::fftw_dft >();
    
    transform_api< std::complex<double>,3,boost::math::fft::gsl_dft >();
    
    transform_api< std::complex<float>,      4,boost::math::fft::bsl_dft >();
    transform_api< std::complex<double>,     4,boost::math::fft::bsl_dft >();
    transform_api< std::complex<long double>,4,boost::math::fft::bsl_dft >();
    
    plan_api< std::complex<double>,      3, boost::math::fft::fftw_dft >();
    plan_api< std::complex<float>,       3, boost::math::fft::fftw_dft >();
    plan_api< std::complex<long double>, 3, boost::math::fft::fftw_dft >();
    
    plan_api< std::complex<double>,3, boost::math::fft::gsl_dft >();
    
    plan_api< std::complex<double>,      4, boost::math::fft::bsl_dft >();
    plan_api< std::complex<float>,       4, boost::math::fft::bsl_dft >();
    plan_api< std::complex<long double>, 4, boost::math::fft::bsl_dft >();
    
    return 0;
}
