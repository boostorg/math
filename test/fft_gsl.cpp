#include "math_unit_test.hpp"
#include <boost/math/fft.hpp>
#include <boost/math/fft/gsl_backend.hpp>
#include <boost/math/constants/constants.hpp>

#include <type_traits>
#include <complex>
#include <vector>
#include <limits>
#include <cmath>
#include <random>

using namespace boost::math::fft;

void test_gsl(int N)
{
    using T = double;
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
        //dft_forward<gsl_dft>(A.cbegin(),A.cend(),B.begin());
        //dft_backward<gsl_dft>(B.cbegin(),B.cend(),C.begin());
        
        dft_forward<gsl_dft>(A.data(),A.data()+A.size(),B.data());
        dft_backward<gsl_dft>(B.data(),B.data()+B.size(),C.data());
        
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
        //dft_forward<gsl_dft>(A.cbegin(),A.cend(),B.begin());
        //dft_backward<gsl_dft>(B.cbegin(),B.cend(),C.begin());
        
        dft_forward<gsl_dft>(A.data(),A.data()+A.size(),B.data());
        dft_backward<gsl_dft>(B.data(),B.data()+B.size(),C.data());
        
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
    for(int i=1;i<=(1<<12); i*=2)
    {
        test_gsl(i);
    }
    for(int i=1;i<=100'000; i*=10)
    {
        test_gsl(i);
    }
    for(auto i : std::vector<int>{3,5,7,11,13,17,23})
    {
        test_gsl(i);
    }
    return boost::math::test::report_errors();
}

