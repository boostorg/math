#include "math_unit_test.hpp"
#include <boost/math/fft.hpp>
#include <boost/math/fft/bsl_backend.hpp>
#include <boost/math/fft/fftw_backend.hpp>
#include <boost/math/fft/gsl_backend.hpp>
#include <boost/math/constants/constants.hpp>

#include <algorithm>
#include <list>
#include <type_traits>
#include <complex>
#include <vector>
#include <limits>
#include <cmath>
#include <random>

using namespace boost::math::fft;

template<class T>
std::vector<std::complex<T>> random_vector(int N)
{
  using local_vector_complex_type = std::vector<std::complex<T>>;

  std::mt19937 rng;
  std::uniform_real_distribution<T> U(0.0,1.0);
  local_vector_complex_type A(N);
  for(auto& x: A)
  {
    x.real( U(rng) );
    x.imag( U(rng) );
  }
  return A;
}
template<class Container1, class Container2>
typename Container1::value_type::value_type difference(const Container1& A, const Container2& B)
{
  BOOST_MATH_ASSERT_MSG( A.size()==B.size(), "Different container sizes.");
  using ComplexType = typename Container1::value_type;
  using RealType = typename ComplexType::value_type;
  const RealType inv_N = RealType{1}/B.size();
  RealType diff =
    std::inner_product(
      A.begin(),A.end(),
      B.begin(),
      RealType{0.0},
      std::plus<RealType>(),
      [inv_N](const ComplexType & a, const ComplexType& b)
      {
        return std::norm(a-b*inv_N);
      });
  diff = std::sqrt(diff)/A.size();
  return diff;
}

template<class RealType, template<class U> class Backend>
void test_inverse(int N, int tolerance)
{
  using ComplexType = std::complex<RealType>;
  const RealType tol = tolerance*std::numeric_limits<RealType>::epsilon();
  const std::vector< ComplexType > A{random_vector<RealType>(N)};

  {
    std::vector<ComplexType> B(N), C(N);
    
    dft_forward<Backend>(std::begin(A),std::end(A),std::begin(B));
    dft_backward<Backend>(std::begin(B),std::end(B),std::begin(C));
    
    RealType diff{difference(A,C)};
    CHECK_MOLLIFIED_CLOSE(RealType{0.0},diff,tol);
  }
  {
    std::vector<ComplexType>  C(N);
    std::list<ComplexType> B;

    dft_forward<Backend>(std::begin(A),std::end(A),std::back_inserter(B));
    dft_backward<Backend>(std::begin(B),std::end(B),std::begin(C));

    RealType diff{difference(A,C)};
    CHECK_MOLLIFIED_CLOSE(RealType{0.0},diff,tol);
  }
  {
    std::list<ComplexType> C;

    dft_forward<Backend>(std::begin(A),std::end(A),std::back_inserter(C));
    dft_backward<Backend>(std::begin(C),std::end(C),std::begin(C));

    RealType diff{difference(A,C)};
    CHECK_MOLLIFIED_CLOSE(RealType{0.0},diff,tol);
  }
}

int main()
{
  for(int i=1;i<=(1<<12); i*=2)
  {
    test_inverse<float,fftw_dft>(i,1);
    test_inverse<double,fftw_dft>(i,1);
    test_inverse<long double,fftw_dft>(i,1);
    
    test_inverse<double,gsl_dft>(i,1);
    
    test_inverse<float,bsl_dft>(i,1);
    test_inverse<double,bsl_dft>(i,1);
    test_inverse<long double,bsl_dft>(i,1);
  }
  return boost::math::test::report_errors();
}
