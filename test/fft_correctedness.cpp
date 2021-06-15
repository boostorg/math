#include "math_unit_test.hpp"
#include <boost/math/fft.hpp>
#include <boost/math/fft/fftw_backend.hpp>
#include <boost/math/fft/gsl_backend.hpp>
#include <boost/math/constants/constants.hpp>

#include <type_traits>
#include <complex>
#include <vector>
#include <limits>
#include <cmath>
#include <random>

using namespace boost::math::fft;

template<class T, template<class U> class Backend>
void test_fixed_transforms()
{
    const T tol = 4*std::numeric_limits<T>::epsilon();
    
    {
        std::vector< std::complex<T> > A{1.0},B(1);
        //dft_forward<Backend>(A.begin(),A.end(),B.begin());
        dft_forward<Backend>(A.data(),A.data()+A.size(),B.data());
        CHECK_MOLLIFIED_CLOSE(T{1.0},B[0].real(),0);
        CHECK_MOLLIFIED_CLOSE(T{0.0},B[0].imag(),0);
    }
    {
        std::vector< std::complex<T> > A{1.0,1.0},B(2);
        //dft_forward<Backend>(A.begin(),A.end(),B.begin());
        dft_forward<Backend>(A.data(),A.data()+A.size(),B.data());
        CHECK_MOLLIFIED_CLOSE(T{2.0},B[0].real(),tol);
        CHECK_MOLLIFIED_CLOSE(T{0.0},B[0].imag(),tol);
        
        CHECK_MOLLIFIED_CLOSE(T{0.0},B[1].real(),tol);
        CHECK_MOLLIFIED_CLOSE(T{0.0},B[1].imag(),tol);
    }
    {
        std::vector< std::complex<T> > A{1.0,1.0,1.0},B(3);
        //dft_forward<Backend>(A.begin(),A.end(),B.begin());
        dft_forward<Backend>(A.data(),A.data()+A.size(),B.data());
        CHECK_MOLLIFIED_CLOSE(T{3.0},B[0].real(),tol);
        CHECK_MOLLIFIED_CLOSE(T{0.0},B[0].imag(),tol);
        
        CHECK_MOLLIFIED_CLOSE(
            T{0.0},B[1].real(),tol);
        CHECK_MOLLIFIED_CLOSE(
            T{0.0},B[1].imag(),tol);
        
        CHECK_MOLLIFIED_CLOSE(
            T{0.0},B[2].real(),tol);
        CHECK_MOLLIFIED_CLOSE(
            T{0.0},B[2].imag(),tol);
    }
    {
        std::vector< std::complex<T> > A{1.0,1.0,1.0};
        //dft_forward<Backend>(A.cbegin(),A.cend(),A.begin());
        dft_forward<Backend>(A.data(),A.data()+A.size(),A.data());
        CHECK_MOLLIFIED_CLOSE(T{3.0},A[0].real(),tol);
        CHECK_MOLLIFIED_CLOSE(T{0.0},A[0].imag(),tol);
        
        CHECK_MOLLIFIED_CLOSE(
            T{0.0},A[1].real(),tol);
        CHECK_MOLLIFIED_CLOSE(
            T{0.0},A[1].imag(),tol);
        
        CHECK_MOLLIFIED_CLOSE(
            T{0.0},A[2].real(),tol);
        CHECK_MOLLIFIED_CLOSE(
            T{0.0},A[2].imag(),tol);
    }
}


template<class T, template<class U> class Backend>
void test_inverse(int N)
{
  const T tol = 4*std::numeric_limits<T>::epsilon();
  std::mt19937 rng;
  std::uniform_real_distribution<T> U(0.0,1.0);
  {
    std::vector<std::complex<T>> A(N),B(N),C(N);
    
    for(auto& x: A)
    {
        x.real( U(rng) );
        x.imag( U(rng) );
    }
    //dft_forward<Backend>(A.cbegin(),A.cend(),B.begin());
    //dft_backward<Backend>(B.cbegin(),B.cend(),C.begin());
    
    dft_forward<Backend>(A.data(),A.data()+A.size(),B.data());
    dft_backward<Backend>(B.data(),B.data()+B.size(),C.data());
    
    const T inverse_N = T{1.0}/N;
    for(auto &x : C)
      x *= inverse_N;
    
    T diff{0.0};
    
    for(size_t i=0;i<A.size();++i)
    {
        diff += std::norm(A[i]-C[i]);
    }
    diff = std::sqrt(diff)*inverse_N;
    CHECK_MOLLIFIED_CLOSE(T{0.0},diff,tol);
  }
  
  {
    std::vector<std::complex<T>> A(N),B(N),C(N);
    
    for(auto& x: A)
    {
      x.real( 1.0 );
      x.imag( 0.0 );
    }
    //dft_forward<Backend>(A.cbegin(),A.cend(),B.begin());
    //dft_backward<Backend>(B.cbegin(),B.cend(),C.begin());
    
    dft_forward<Backend>(A.data(),A.data()+A.size(),B.data());
    dft_backward<Backend>(B.data(),B.data()+B.size(),C.data());
    
    const T inverse_N = T{1.0}/N;
    for(auto &x : C)
        x *= inverse_N;
    
    T diff{0.0};
    
    for(size_t i=0;i<A.size();++i)
    {
        diff += std::norm(A[i]-C[i]);
    }
    diff = std::sqrt(diff)*inverse_N;
    CHECK_MOLLIFIED_CLOSE(T{0.0},diff,tol);
  }
}


int main()
{
  test_fixed_transforms<float,fftw_dft>();
  test_fixed_transforms<double,fftw_dft>();
  test_fixed_transforms<long double,fftw_dft>();
  
  test_fixed_transforms<double,gsl_dft>();
  
  //test_fixed_transforms<float,bsl_dft>();
  //test_fixed_transforms<double,bsl_dft>();
  //test_fixed_transforms<long double,bsl_dft>();
  
  for(int i=1;i<=(1<<12); i*=2)
  {
    test_inverse<float,fftw_dft>(i);
    test_inverse<double,fftw_dft>(i);
    test_inverse<long double,fftw_dft>(i);
    
    test_inverse<double,gsl_dft>(i);
    
    test_inverse<float,bsl_dft>(i);
    test_inverse<double,bsl_dft>(i);
    test_inverse<long double,bsl_dft>(i);
  }
  for(int i=1;i<=100'000; i*=10)
  {
    test_inverse<float,fftw_dft>(i);
    test_inverse<double,fftw_dft>(i);
    test_inverse<long double,fftw_dft>(i);
    
    test_inverse<double,gsl_dft>(i);
    
  //  test_inverse<float,bsl_dft>(i);
  //  test_inverse<double,bsl_dft>(i);
  //  test_inverse<long double,bsl_dft>(i);
  }
  for(auto i : std::vector<int>{3,5,7,11,13,17,23})
  {
    test_inverse<float,fftw_dft>(i);
    test_inverse<double,fftw_dft>(i);
    test_inverse<long double,fftw_dft>(i);
    
    test_inverse<double,gsl_dft>(i);
    
  //  test_inverse<float,bsl_dft>(i);
  //  test_inverse<double,bsl_dft>(i);
  //  test_inverse<long double,bsl_dft>(i);
  }
  // TODO: can we print a useful compilation error message for the following
  // illegal case?
  // dft<std::complex<int>> P(3);   
  return boost::math::test::report_errors();
}
