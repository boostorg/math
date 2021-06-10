#include "math_unit_test.hpp"
#include <boost/math/fft.hpp>
#include <boost/math/constants/constants.hpp>


#include <complex>
#include <vector>
#include <limits>
#include <cmath>
#include <random>

using namespace boost::math::fft;

template<class T>
void test_fixed_transforms()
{
    const T pi = boost::math::constants::pi<T>(); 
    const T tol = 4*std::numeric_limits<T>::epsilon();
    
    {
        std::vector< std::complex<T> > A{1.0},B(1);
        dft_forward(A.begin(),A.end(),B.begin());
        CHECK_MOLLIFIED_CLOSE(B[0].real(),T(1.0),0);
        CHECK_MOLLIFIED_CLOSE(B[0].imag(),T(0.0),0);
    }
    {
        std::vector< std::complex<T> > A{1.0,1.0},B(2);
        dft_forward(A.begin(),A.end(),B.begin());
        CHECK_MOLLIFIED_CLOSE(B[0].real(),T(2.0),tol);
        CHECK_MOLLIFIED_CLOSE(B[0].imag(),T(0.0),tol);
        
        CHECK_MOLLIFIED_CLOSE(B[1].real(),T(0.0),tol);
        CHECK_MOLLIFIED_CLOSE(B[1].imag(),T(0.0),tol);
    }
    {
        std::vector< std::complex<T> > A{1.0,1.0,1.0},B(3);
        dft_forward(A.begin(),A.end(),B.begin());
        CHECK_MOLLIFIED_CLOSE(B[0].real(),T(3.0),tol);
        CHECK_MOLLIFIED_CLOSE(B[0].imag(),T(0.0),tol);
        
        CHECK_MOLLIFIED_CLOSE(
            B[1].real(),
            T(1.0 + std::cos(2*pi/3.0) + std::cos(4*pi/3.0)),
            tol);
        CHECK_MOLLIFIED_CLOSE(
            B[1].imag(),
            T(- std::sin(2*pi/3.0) - std::sin(4*pi/3.0)),
            tol);
        
        CHECK_MOLLIFIED_CLOSE(
            B[2].real(),
            T(1.0 + std::cos(4*pi/3.0) + std::cos(8*pi/3.0)),
            tol);
        CHECK_MOLLIFIED_CLOSE(
            B[2].imag(),
            T(- std::sin(4*pi/3.0) - std::sin(8*pi/3.0)),
            tol);
    }
}

template<class T>
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
        dft_forward(A.cbegin(),A.cend(),B.begin());
        dft_backward(B.cbegin(),B.cend(),C.begin());
        
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
        dft_forward(A.cbegin(),A.cend(),B.begin());
        dft_backward(B.cbegin(),B.cend(),C.begin());
        
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
    test_fixed_transforms<float>();
    test_fixed_transforms<double>();
    test_fixed_transforms<long double>();
   
    for(int i=1;i<=(1<<12); i*=2)
    {
        test_inverse<float>(i);
        test_inverse<double>(i);
        test_inverse<long double>(i);
    
    }
    for(int i=1;i<=100'000; i*=10)
    {
        test_inverse<float>(i);
        test_inverse<double>(i);
        test_inverse<long double>(i);
    
    }
    for(auto i : std::vector<int>{3,5,7,11,13,17,23})
    {
        test_inverse<float>(i);
        test_inverse<double>(i);
        test_inverse<long double>(i);
    
    }
    return boost::math::test::report_errors();
}
