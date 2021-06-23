#include "math_unit_test.hpp"
#include <boost/math/fft.hpp>
#include <boost/math/fft/fftw_backend.hpp>
#include <boost/math/fft/gsl_backend.hpp>
#include <boost/math/fft/algorithms.hpp>
#include <boost/math/fft/abstract_ring.hpp>
#include <boost/math/constants/constants.hpp>
#ifdef BOOST_MATH_USE_FLOAT128
#include <boost/multiprecision/complex128.hpp>
#endif
#include <boost/multiprecision/cpp_complex.hpp>
#include <boost/random.hpp>

#include "fft_test_helpers.hpp"

#include <type_traits>
#include <vector>
#include <limits>

using namespace boost::math::fft;



template<class T, template<class U> class Backend>
void test_fixed_transforms()
{
  using Complex = typename detail::select_complex<T>::type;
  // TODO: increase precision of the generic dft 
  // const T tol = std::numeric_limits<T>::epsilon();
    const T tol = 8*std::numeric_limits<T>::epsilon();
    {
        std::vector< Complex > A{1.0},B(1);
        dft_forward<Backend>(A.data(),A.data()+A.size(),B.data());
        CHECK_MOLLIFIED_CLOSE(T{1.0},B[0].real(),0);
        CHECK_MOLLIFIED_CLOSE(T{0.0},B[0].imag(),0);
    }
    {
        std::vector< Complex > A{1.0,1.0},B(2);
        dft_forward<Backend>(A.data(),A.data()+A.size(),B.data());
        CHECK_MOLLIFIED_CLOSE(T{2.0},B[0].real(),tol);
        CHECK_MOLLIFIED_CLOSE(T{0.0},B[0].imag(),tol);
        
        CHECK_MOLLIFIED_CLOSE(T{0.0},B[1].real(),tol);
        CHECK_MOLLIFIED_CLOSE(T{0.0},B[1].imag(),tol);
    }
    {
        std::vector< Complex > A{1.0,1.0,1.0},B(3);
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
        std::vector< Complex > A{1.0,1.0,1.0};
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
  using Complex = typename detail::select_complex<T>::type;
  // TODO: increase precision of the generic dft 
  // const T tol = std::numeric_limits<T>::epsilon();
  const T tol = 128*std::numeric_limits<T>::epsilon();
  
  boost::random::mt19937 rng;
  boost::random::uniform_real_distribution<T> U(0.0,1.0);
  {
    std::vector<Complex> A(N),B(N),C(N);
    
    for(auto& x: A)
    {
        x.real( U(rng) );
        x.imag( U(rng) );
    }
    dft_forward<Backend>(A.data(),A.data()+A.size(),B.data());
    dft_backward<Backend>(B.data(),B.data()+B.size(),C.data());
    
    const T inverse_N = T{1.0}/N;
    for(auto &x : C)
      x *= inverse_N;
    
    T diff{0.0};
    
    for(size_t i=0;i<A.size();++i)
    {
        using std::norm;
        diff += norm(A[i]-C[i]);
    }
    using std::sqrt;
    diff = sqrt(diff)*inverse_N;
    CHECK_MOLLIFIED_CLOSE(T{0.0},diff,tol);
  }
}


int main()
{
  test_fixed_transforms<float,fftw_dft>();
  test_fixed_transforms<double,fftw_dft>();
  test_fixed_transforms<long double,fftw_dft>();
#ifdef BOOST_MATH_USE_FLOAT128
  test_fixed_transforms<boost::multiprecision::float128,fftw_dft>();
#endif
  
  test_fixed_transforms<double,gsl_dft>();
  
  test_fixed_transforms<float,bsl_dft>();
  test_fixed_transforms<double,bsl_dft>();
  test_fixed_transforms<long double,bsl_dft>();
#ifdef BOOST_MATH_USE_FLOAT128
  test_fixed_transforms<boost::multiprecision::float128,bsl_dft>();
#endif
  test_fixed_transforms<boost::multiprecision::cpp_bin_float_50,bsl_dft>();
  test_fixed_transforms<boost::multiprecision::cpp_bin_float_100,bsl_dft>();
  test_fixed_transforms<boost::multiprecision::cpp_bin_float_quad,bsl_dft>();
  // TODO:
  //test_fixed_transforms<boost::multiprecision::mpfr_float_100,bsl_dft>();
  
  for(int i=1;i<=(1<<10); i*=2)
  {
    test_inverse<float,fftw_dft>(i);
    test_inverse<double,fftw_dft>(i);
    test_inverse<long double,fftw_dft>(i);
#ifdef BOOST_MATH_USE_FLOAT128
    test_inverse<boost::multiprecision::float128,fftw_dft>(i);
#endif
    
    test_inverse<double,gsl_dft>(i);
    
    test_inverse<float,bsl_dft>(i);
    test_inverse<double,bsl_dft>(i);
    test_inverse<long double,bsl_dft>(i);
#ifdef BOOST_MATH_USE_FLOAT128
    test_inverse<boost::multiprecision::float128,bsl_dft>(i);
#endif
    test_inverse<boost::multiprecision::cpp_bin_float_50,bsl_dft>(i);
    
    test_inverse<float,test_dft_power2_dit>(i);
    test_inverse<float,test_dft_power2_dif>(i);
    
    test_inverse<double,test_dft_power2_dit>(i);
    test_inverse<double,test_dft_power2_dif>(i);
    
    test_inverse<long double,test_dft_power2_dit>(i);
    test_inverse<long double,test_dft_power2_dif>(i);

#ifdef BOOST_MATH_USE_FLOAT128
    test_inverse<boost::multiprecision::float128,test_dft_power2_dit>(i);
    test_inverse<boost::multiprecision::float128,test_dft_power2_dif>(i);
#endif

    test_inverse<boost::multiprecision::cpp_bin_float_50,test_dft_power2_dit>(i);
    test_inverse<boost::multiprecision::cpp_bin_float_50,test_dft_power2_dif>(i);
  }
  for(int i=1;i<=1000; i*=10)
  {
    test_inverse<float,fftw_dft>(i);
    test_inverse<double,fftw_dft>(i);
    test_inverse<long double,fftw_dft>(i);
    if(i <=10) {
#ifdef BOOST_MATH_USE_FLOAT128
      test_inverse<boost::multiprecision::float128,fftw_dft>(i);
#endif
    }
    
    test_inverse<double,gsl_dft>(i);
    
    test_inverse<float,bsl_dft>(i);
    test_inverse<double,bsl_dft>(i);
    test_inverse<long double,bsl_dft>(i);
    if(i <=10) {
#ifdef BOOST_MATH_USE_FLOAT128
      test_inverse<boost::multiprecision::float128,bsl_dft>(i);
#endif
      test_inverse<boost::multiprecision::cpp_bin_float_50,bsl_dft>(i);
    }
  }
  for(auto i : std::vector<int>{2,3,5,7,11,13,17,23,29,31})
  {
    test_inverse<float,fftw_dft>(i);
    test_inverse<double,fftw_dft>(i);
    test_inverse<long double,fftw_dft>(i);
#ifdef BOOST_MATH_USE_FLOAT128
    test_inverse<boost::multiprecision::float128,fftw_dft>(i);
#endif
    
    test_inverse<double,gsl_dft>(i);
    
    test_inverse<float,bsl_dft>(i);
    test_inverse<double,bsl_dft>(i);
    test_inverse<long double,bsl_dft>(i);
#ifdef BOOST_MATH_USE_FLOAT128
    test_inverse<boost::multiprecision::float128,bsl_dft>(i);
#endif
    test_inverse<boost::multiprecision::cpp_bin_float_50,bsl_dft>(i);
    
    test_inverse<float,test_dft_generic_prime_bruteForce>(i);
    test_inverse<double,test_dft_generic_prime_bruteForce>(i);
    test_inverse<long double,test_dft_generic_prime_bruteForce>(i);
#ifdef BOOST_MATH_USE_FLOAT128
    test_inverse<boost::multiprecision::float128,test_dft_generic_prime_bruteForce>(i);
#endif
    test_inverse<boost::multiprecision::cpp_bin_float_50,test_dft_generic_prime_bruteForce>(i);
    
    test_inverse<float,test_dft_complex_prime_bruteForce>(i);
    test_inverse<double,test_dft_complex_prime_bruteForce>(i);
    test_inverse<long double,test_dft_complex_prime_bruteForce>(i);
#ifdef BOOST_MATH_USE_FLOAT128
    test_inverse<boost::multiprecision::float128,test_dft_complex_prime_bruteForce>(i);
#endif
    test_inverse<boost::multiprecision::cpp_bin_float_50,test_dft_complex_prime_bruteForce>(i);
  }
  // TODO: can we print a useful compilation error message for the following
  // illegal case?
  // dft<std::complex<int>> P(3);   
  
  for(int i=1;i<=100;++i)
  {
    test_inverse<float,test_dft_generic_composite>(i);
    test_inverse<double,test_dft_generic_composite>(i);
    test_inverse<long double,test_dft_generic_composite>(i);
  }
  
  return boost::math::test::report_errors();
}
