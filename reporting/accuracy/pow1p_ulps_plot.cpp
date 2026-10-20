//  (C) Copyright Nick Thompson 2026.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// Writes the ULPs plots of pow1p used in doc/sf/pow1p.qbk.
// pow1p is drawn in steelblue and the naive pow(1 + x, y) in orange.

#include <cmath>
#include <string>
#include <boost/math/special_functions/pow1p.hpp>
#include <boost/math/tools/ulps_plot.hpp>
#include <boost/multiprecision/cpp_bin_float.hpp>

using boost::math::tools::ulps_plot;
using PreciseReal = boost::multiprecision::cpp_bin_float_100;
using CoarseReal = double;

void plot(CoarseReal y, CoarseReal a, CoarseReal b, std::string const & title, std::string const & filename)
{
   auto precise = [y](PreciseReal const & x) { return exp(PreciseReal(y) * log1p(x)); };
   auto pow1p = [y](CoarseReal x) { return boost::math::pow1p(x, y); };
   auto naive = [y](CoarseReal x) { return std::pow(1 + x, y); };

   auto p = ulps_plot<decltype(precise), PreciseReal, CoarseReal>(precise, a, b, 2000, false, 42);
   p.clip(100).width(700).vertical_lines(6);
   p.background_color("white").font_color("black");
   p.title(title);
   p.add_fn(naive, "orange");
   p.add_fn(pow1p, "steelblue");
   p.write(filename);
}

int main()
{
   plot(1000, -0.5, 0.5, "pow1p(x, 1000) vs pow(1 + x, 1000), double", "pow1p_ulps_y1000.svg");
   plot(1048576, -1e-4, 1e-4, "pow1p(x, 2^20) vs pow(1 + x, 2^20), double", "pow1p_ulps_y2e20.svg");
}
