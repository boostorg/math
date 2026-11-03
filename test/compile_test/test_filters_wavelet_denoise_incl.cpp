//  Copyright Nick Thompson 2026.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// Basic sanity check that header <boost/math/filters/wavelet_denoise.hpp>
// #includes all the files that it needs to.
//
#include <boost/math/filters/wavelet_denoise.hpp>
#include <vector>

void compile_and_link_test()
{
   using namespace boost::math::filters;
   const std::vector<double> y(64, 1.0);
   const orthogonal_wavelet<double> w = wavelet::symlet<double>(4);
   std::vector<double> x = wavelet_denoise(y, w);
   x = wavelet_denoise(y, wavelet::daubechies<double>(4), threshold_rule::sure, 0.5);
   x = wavelet_denoise(y);
   x = denoise(y);
   x = inverse_wavelet_transform(wavelet_transform(y, w), w);
   x = wavelet_transform(y, w, 2);
   const std::vector<float> yf(64, 1.0f);
   std::vector<float> xf = wavelet_denoise(yf, wavelet::daubechies<float>(2), threshold_rule::wiener);
   xf = denoise(yf);
}
