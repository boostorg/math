///////////////////////////////////////////////////////////////////
//  Copyright Eduardo Quintana 2021
//  Copyright Janek Kozicki 2021
//  Copyright Christopher Kormanyos 2021
//  Distributed under the Boost Software License,
//  Version 1.0. (See accompanying file LICENSE_1_0.txt
//  or copy at http://www.boost.org/LICENSE_1_0.txt)
/*
    boost::math::fft example 05
    
    several engines, different complex types.
*/

#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/multiprecision/cpp_bin_float.hpp>
#include <boost/multiprecision/cpp_complex.hpp>
//#include <boost/multiprecision/mpc.hpp>
#ifdef BOOST_MATH_USE_FLOAT128
#include <boost/multiprecision/complex128.hpp>
#endif

#include <boost/math/fft.hpp>
#include <boost/math/fft/fftw_backend.hpp>
#include <boost/math/fft/gsl_backend.hpp>
#include <boost/core/demangle.hpp>
#include <iostream>
#include <vector>
#include <complex>

template<class T>
void print(const std::vector< T >& V)
{
    for(auto i=0UL;i<V.size();++i)
        std::cout << "V[" << i << "] = " << std::setprecision(std::numeric_limits<typename T::value_type>::digits10 + 4)
            << V[i].real() << ", " << V[i].imag() << '\n';
}
template<class Complex>
void test_bsl() {
    std::cout << "BSL engine with " << boost::core::demangle(typeid(Complex).name()) << "\n";
    std::cout << "Real type is    " << boost::core::demangle(typeid(typename Complex::value_type).name()) << "\n";
    std::vector< Complex > A{1.0,2.0,3.0,4.0},B(A.size());
    // forward transform, out-of-place
    boost::math::fft::dft_forward(A.cbegin(),A.cend(),B.begin());
    print(B);
    // backward transform, in-place
    boost::math::fft::dft_backward(B.cbegin(),B.cend(),B.begin());
    print(B);
}

template<class Complex>
void test_fftw() {
    std::cout << "FFTW engine with " << boost::core::demangle(typeid(Complex).name()) << "\n";
    std::vector< Complex > A{1.0,2.0,3.0,4.0},B(A.size());
    // forward transform, out-of-place
    boost::math::fft::dft_forward<boost::math::fft::fftw_dft>(A.cbegin(),A.cend(),B.begin());
    print(B);
    // backward transform, in-place
    boost::math::fft::dft_backward<boost::math::fft::fftw_dft>(B.cbegin(),B.cend(),B.begin());
    print(B);
}

template<class Complex>
void test_gsl() {
    std::cout << "GSL engine with " << boost::core::demangle(typeid(Complex).name()) << "\n";
    std::vector< Complex > A{1.0,2.0,3.0,4.0},B(A.size());
    // forward transform, out-of-place
    boost::math::fft::dft_forward<boost::math::fft::gsl_dft>(A.cbegin(),A.cend(),B.begin());
    print(B);
    // backward transform, in-place
    boost::math::fft::dft_backward<boost::math::fft::gsl_dft>(B.cbegin(),B.cend(),B.begin());
    print(B);
}

int main()
{
    test_bsl<std::complex<float>>();
    test_bsl<std::complex<double>>();
    test_bsl<std::complex<long double>>();
#ifdef BOOST_MATH_USE_FLOAT128
    test_bsl<boost::multiprecision::complex128>();
#endif
    test_bsl<boost::multiprecision::cpp_complex_50>();
    test_bsl<boost::multiprecision::cpp_complex_100>();
    test_bsl<boost::multiprecision::cpp_complex_quad>();
//  test_bsl<boost::multiprecision::number<boost::multiprecision::complex_adaptor<boost::multiprecision::cpp_dec_float<50>>, boost::multiprecision::et_off>>();
//  test_bsl<boost::multiprecision::mpc_complex_50>();
//using aa = boost::math::fft::detail::select_complex<void>::type;
//std::cout << "AA " << boost::core::demangle(typeid(aa).name()) << "\n";
    test_fftw<std::complex<float>>();
    test_fftw<std::complex<double>>();
    test_fftw<std::complex<long double>>();
#ifdef BOOST_MATH_USE_FLOAT128
    test_fftw<boost::multiprecision::complex128>();
#endif

    test_gsl<std::complex<double>>();
    return 0;
}

