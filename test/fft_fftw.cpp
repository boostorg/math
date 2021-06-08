#include "math_unit_test.hpp"
#include <boost/math/fft.hpp>
#include <boost/math/constants/constants.hpp>

#include <vector>
#include <limits>
#include <cmath>

using namespace boost::math::fft;

template<class T>
void test_fixed_transforms()
{
    const T pi = boost::math::constants::pi<T>(); 
    const T tol = 4*std::numeric_limits<T>::epsilon();
    
    // { TODO: error, fftw does nothing for n=1
    //     std::vector< std::complex<T> > A{1.0},B(1);
    //     dft_forward(A.begin(),A.end(),B.begin());
    //     CHECK_MOLLIFIED_CLOSE(B[0].real(),T(1.0),0);
    //     CHECK_MOLLIFIED_CLOSE(B[0].imag(),T(0.0),0);
    // }
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

int main()
{
    test_fixed_transforms<float>();
    test_fixed_transforms<double>();
    test_fixed_transforms<long double>();
    return boost::math::test::report_errors();
}
