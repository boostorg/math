#include <boost/math/fft.hpp>
#include <boost/math/fft/bsl_backend.hpp>

#include <complex>
#include <vector>
#include <array>

void compile_test()
{   
    // test same type of iterator
    std::vector<int> A(3);
    std::vector<double> B(A.size());
    
    // fails static_assert: D and A types are different
    boost::math::fft::dft(A.begin(),A.end(),B.begin());
}

int main()
{
    compile_test();
    return 0;
}
