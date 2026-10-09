//  (C) Copyright Nick Thompson 2026.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// Writes 2D ulps plots of ligamma(a, x) in double over a log-log grid.
// The reference is log(tgamma(a, x)) in 50 digits, which can't underflow.
// Requires lodepng: link lodepng.cpp and put lodepng.h on the include path.

#include <iostream>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/tools/ulps_heatmap.hpp>
#include <boost/multiprecision/cpp_bin_float.hpp>

using boost::math::tools::ulps_heatmap;
using PreciseReal = boost::multiprecision::cpp_bin_float_50;
using CoarseReal = double;

int main()
{
   auto precise = [](PreciseReal const & a, PreciseReal const & x) { return log(boost::math::tgamma(a, x)); };
   auto ligamma = [](CoarseReal a, CoarseReal x) { return boost::math::ligamma(a, x); };

   auto h = ulps_heatmap<decltype(precise), PreciseReal, CoarseReal>(precise, 1e-3, 1e4, 1e-3, 1e4, 600, 400, true, true);
   h.add_fn(ligamma);
   std::cout << h.summary();
   h.x_label("a").y_label("x");

   h.title("ligamma(a, x), double: error in ulps");
   h.log_scale(false).write_ulps("ligamma_ulps.svg");
   h.log_scale(true).write_ulps("ligamma_ulps_log.svg");

   h.title("ligamma(a, x), double: error relative to the conditioning envelope");
   h.log_scale(false).write_ulps_envelope("ligamma_ulps_envelope.svg");
   h.log_scale(true).write_ulps_envelope("ligamma_ulps_envelope_log.svg");
}
