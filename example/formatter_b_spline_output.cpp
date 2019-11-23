//  (C) Copyright John Maddock 2019.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <boost/math/tools/formatting.hpp>
#include <iostream>
#include <vector>
#include <boost/math/interpolators/cardinal_cubic_b_spline.hpp>

int main() {
   std::vector<double> v(10);
   for (size_t i = 0; i < v.size(); ++i)
   {
      v[i] = 5 + 6 * i;
   }

   double step = 0.1;
   double a = 5;
   boost::math::interpolators::cardinal_cubic_b_spline<double> spline(v.data(), v.size(), a, step);

   boost::math::tools::text_printer printer(std::cout);
   printer << spline << std::endl;

   boost::math::tools::latex_printer lprinter(std::cout);
   lprinter << spline << std::endl;

   boost::math::tools::html_printer hprinter(std::cout);
   hprinter << spline << std::endl;

   return 0;
}
