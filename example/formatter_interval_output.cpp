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
   typedef boost::multiprecision::number<boost::multiprecision::mpfi_float_backend<20> > mp_type;
   mp_type f = boost::math::constants::pi<mp_type>();

   std::cout << std::setprecision(std::numeric_limits<mp_type>::digits10);

   boost::math::tools::text_printer t_print(std::cout);
   t_print << f << std::endl;

   boost::math::tools::latex_printer l_print(std::cout);
   l_print << f << std::endl;

   boost::math::tools::docbook_printer d_print(std::cout);
   d_print << f << std::endl;

   boost::math::tools::polynomial<mp_type> poly1 = { 2, 3, 4, -5, -4, 12 };
   boost::math::tools::polynomial<mp_type> poly2 = { -3, 2, 4, 5 };
   poly1 *= poly2;
   poly1 /= 3;

   t_print << poly1 << std::endl;

   l_print << poly1 << std::endl;

   d_print << poly1 << std::endl;

   return 0;
}
