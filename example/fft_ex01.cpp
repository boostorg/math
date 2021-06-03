#include <boost/math/fft.hpp>

#include <complex>
#include <vector>

int main()
{
    {
        std::vector<int> A{1,2,3},B(A.size());
        boost::math::fft::dft(A.begin(),A.end(),B.begin());
        boost::math::fft::dft(A.data(),A.data()+A.size(),B.data());
    }
    {
        std::vector< std::complex<int> > A{1,2,3},B(A.size());
        boost::math::fft::dft(A.begin(),A.end(),B.begin());
        boost::math::fft::dft(A.data(),A.data()+A.size(),B.data());
    }
    {
        std::vector< std::complex<double> > A{1,2,3},B(A.size());
        boost::math::fft::dft(A.begin(),A.end(),B.begin());
        boost::math::fft::dft(A.data(),A.data()+A.size(),B.data());
    }
    return 0;
}
