//  (C) Copyright John Maddock 2019.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <boost/math/tools/formatting.hpp>
#include <iostream>
#include <vector>
#include <boost/math/interpolators/cardinal_cubic_b_spline.hpp>
#include <boost/math/interpolators/cardinal_quadratic_b_spline.hpp>
#include <boost/math/interpolators/cardinal_quintic_b_spline.hpp>

int main() {
   std::vector<double> v(10);
   for (size_t i = 0; i < v.size(); ++i)
   {
      v[i] = static_cast<double>(5 + 6 * i);
   }

   double step = 0.1;
   double a = 5;
   {
      std::cout << "Quadratic B Spline:\n\n";
      boost::math::interpolators::cardinal_quadratic_b_spline<double> spline(v.data(), v.size(), a, step);

      boost::math::tools::text_printer printer(std::cout);
      printer << spline << std::endl << std::endl;

      boost::math::tools::latex_printer lprinter(std::cout);
      lprinter << spline << std::endl << std::endl;

      boost::math::tools::html_printer hprinter(std::cout);
      hprinter << spline << std::endl << std::endl;

      boost::math::tools::docbook_printer dprinter(std::cout);
      dprinter << spline << std::endl << std::endl;
   }
   {
      std::cout << "Cubic B Spline:\n\n";
      boost::math::interpolators::cardinal_cubic_b_spline<double> spline(v.data(), v.size(), a, step);

      boost::math::tools::text_printer printer(std::cout);
      printer << spline << std::endl << std::endl;

      boost::math::tools::latex_printer lprinter(std::cout);
      lprinter << spline << std::endl << std::endl;

      boost::math::tools::html_printer hprinter(std::cout);
      hprinter << spline << std::endl << std::endl;

      boost::math::tools::docbook_printer dprinter(std::cout);
      dprinter << spline << std::endl << std::endl;
   }
   {
      std::cout << "Quintic B Spline:\n\n";
      boost::math::interpolators::cardinal_quintic_b_spline<double> spline(v.data(), v.size(), a, step);

      boost::math::tools::text_printer printer(std::cout);
      printer << spline << std::endl << std::endl;

      boost::math::tools::latex_printer lprinter(std::cout);
      lprinter << spline << std::endl << std::endl;

      boost::math::tools::html_printer hprinter(std::cout);
      hprinter << spline << std::endl << std::endl;

      boost::math::tools::docbook_printer dprinter(std::cout);
      dprinter << spline << std::endl << std::endl;
   }
   return 0;
}
