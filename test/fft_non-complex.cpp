
///////////////////////////////////////////////////////////////////
//  Copyright Eduardo Quintana 2021
//  Copyright Janek Kozicki 2021
//  Copyright Christopher Kormanyos 2021
//  Distributed under the Boost Software License,
//  Version 1.0. (See accompanying file LICENSE_1_0.txt
//  or copy at http://www.boost.org/LICENSE_1_0.txt)
/*
  boost::math::fft test for non-complex types
  Use of DFT for Number Theoretical Transform.
*/
#include <boost/math/fft.hpp>
#include <boost/math/fft/bsl_backend.hpp>

namespace fft = boost::math::fft;

#include <iostream>
#include <vector>
#include "math_unit_test.hpp"
#include "fft_test_helpers.hpp"

class Z337
{
public:
  typedef int integer;
  static constexpr integer mod{337};
};

void test_inverse()
{
  using M_int = fft::my_modulo_lib::mint<Z337>;
  const M_int w{85};
  const M_int inv_8{M_int{8}.inverse()};

  std::vector<M_int> A{4, 3, 2, 1, 0, 0, 0, 0};
  std::vector<M_int> FT_A,FT_FT_A;

  fft::dft_forward<boost::math::fft::bsl_dft>(A.cbegin(),A.cend(),std::back_inserter(FT_A),w);

  fft::dft_backward<boost::math::fft::bsl_dft>(FT_A.cbegin(),FT_A.cend(),std::back_inserter(FT_FT_A),w);

  std::transform(FT_FT_A.begin(), FT_FT_A.end(), FT_FT_A.begin(),
                 [&inv_8](M_int x) { return x * inv_8; });

  int diff = 0;
  for (size_t i = 0; i < A.size(); ++i)
      diff += A[i] == FT_FT_A[i] ? 0 : 1;
  CHECK_EQUAL(0,diff);
}
void test_convolution()
/*
  product of two integer by means of the NTT,
  using the convolution theorem
*/
{
  typedef fft::my_modulo_lib::field_modulo<int, 337> Z337;
  using M_int = fft::my_modulo_lib::mint<Z337>;
  const M_int w{85};
  const M_int inv_8{M_int{8}.inverse()};

  // Multiplying 1234 times 5678 = 7006652
  std::vector<M_int> A{4, 3, 2, 1, 0, 0, 0, 0};
  std::vector<M_int> B{8, 7, 6, 5, 0, 0, 0, 0};

  // forward FFT
  fft::dft_forward<boost::math::fft::bsl_dft>(A.cbegin(),A.cend(),A.begin(), w);
  fft::dft_forward<boost::math::fft::bsl_dft>(B.cbegin(),B.cend(),B.begin(), w);

  // convolution in Fourier space
  std::vector<M_int> AB;
  std::transform(A.begin(), A.end(), B.begin(),
                 std::back_inserter(AB),
                 [](M_int x, M_int y) { return x * y; });

  // backwards FFT
  fft::dft_backward<boost::math::fft::bsl_dft>(AB.cbegin(),AB.cend(),AB.begin(),w);
  std::transform(AB.begin(), AB.end(), AB.begin(),
                 [&inv_8](M_int x) { return x * inv_8; });

  // carry the remainders in base 10
  std::vector<int> C;
  M_int r{0};
  for (auto x : AB)
  {
    auto y = x + r;
    C.emplace_back(int(y) % 10);
    r = M_int(int(y) / 10);
  }
  // yields 7006652
  CHECK_EQUAL(8,static_cast<int>(C.size()));
  CHECK_EQUAL(2,C[0]);
  CHECK_EQUAL(5,C[1]);
  CHECK_EQUAL(6,C[2]);
  CHECK_EQUAL(6,C[3]);
  CHECK_EQUAL(0,C[4]);
  CHECK_EQUAL(0,C[5]);
  CHECK_EQUAL(7,C[6]);
  CHECK_EQUAL(0,C[7]);
}
int main()
{
  test_inverse();
  test_convolution();
  return boost::math::test::report_errors();
}
