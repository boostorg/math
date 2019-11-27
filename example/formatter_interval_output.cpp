//  (C) Copyright John Maddock 2019.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <boost/math/tools/formatting.hpp>
#include <iostream>
#include <boost/multiprecision/mpfi.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/math/tools/polynomial.hpp>

int main() 
{
   boost::multiprecision::mpfi_float_50 f = boost::math::constants::pi<boost::multiprecision::mpfi_float_50>();

   std::cout << std::setprecision(std::numeric_limits<boost::multiprecision::mpfi_float_50>::digits10);

   boost::math::tools::text_printer t_print(std::cout);
   t_print << f << std::endl;

   boost::math::tools::latex_printer l_print(std::cout);
   l_print << f << std::endl;

   boost::math::tools::docbook_printer d_print(std::cout);
   d_print << f << std::endl;

   boost::math::tools::polynomial<boost::multiprecision::mpfi_float_50> poly1 = { 2, 3, 4, -5, -4, 12 };
   boost::math::tools::polynomial<boost::multiprecision::mpfi_float_50> poly2 = { -3, 2, 4, 5 };
   poly1 *= poly2;
   poly1 /= 3;

   t_print << poly1 << std::endl;

   l_print << poly1 << std::endl;

   d_print << poly1 << std::endl;

   return 0;
}
