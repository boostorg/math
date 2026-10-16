//  Copyright Nick Thompson 2026.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// Basic sanity check that header
// #includes all the files that it needs to.
//
#include <boost/math/tools/fft.hpp>
//
// Note this header includes no other headers, this is
// important if this test is to be meaningful:
//
#include "test_compile_result.hpp"

void compile_and_link_test()
{
   boost::math::tools::fft<std::complex<double>> fd(12);
   std::vector<std::complex<double>> v(12);
   fd.transform(v);
   fd.inverse_transform(v);
   check_result<double>(v[0].real());

   boost::math::tools::real_fft<std::complex<float>> rf(8);
   std::vector<float> x(8);
   std::vector<std::complex<float>> X(5);
   rf.transform(x, X);
   check_result<float>(X[0].real());
}
