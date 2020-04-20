//  (C) Copyright John Maddock 2019.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <boost/math/tools/formatting.hpp>
#include <complex>
#include <fstream>
#include <boost/math/tools/polynomial.hpp>


int main(int argc, const char* argv[])
{

   //[formatting_eg1
   //` By way of example, lets suppose we have a polynomial with complex coefficients:
   boost::math::tools::polynomial<std::complex<double> > poly;
   //` And then do something to populate those coefficients...
   //<-
   poly = boost::math::tools::polynomial<std::complex<double> >{ {2.25, -3.5 },  {-12.5, -4.5 }, { 23.34 }, { 0, 34.5 } };
   //->
   //` For debugging purposes, we may wish to output the polynomial to std::cout.
   //` Using `text_printer` we can do just that:

   boost::math::tools::text_printer printer(std::cout);
   printer << poly << std::endl;
   
   //` By default, `text_printer` outputs using Unicode, so even console output is properly formatted:
   //`  
   //` [$../images/console1.png]
   //`     
   //` Note however, that not all platforms support Unicode console output, so if required
   //` we can use the `ascii_text_output` format manipulator to switch to plain ASCII output:

   printer << boost::math::tools::ascii_text_output << poly << std::endl;

   //` The result is rather less pretty, but the best we can do under the circumstances:
   //` 
   //` [$../images/console2.png]
   //`
   //]
   {
      //[formatting_eg2
      //` 
      //`Later, lets suppose we've now finalised our polynomial and wish to format it for print,
      //`in that case we can simply change the printer to `latex_printer` and stream out as before:
      //`
      boost::math::tools::latex_printer printer(std::cout);
      printer << poly << std::endl;
      //` 
      //` After running the output through your favorite Tex processor the result is a nicely formatted
      //` equation:
      //`  
      //` [$../images/polynomial_latex_1.svg]
      //`  
      //]
   }
   {
      //[formatting_eg3
      //` 
      //` The output though may not quite be as we would like, let's suppose
      //` instead we wish to write a complete LaTex document to file and change the
      //` formatting so that the values are printed to 2 decimal places, and the
      //` imaginary unit i is italisized as is common in many historical texts.
      //` The code now looks like:
      //`  
      std::ofstream fs("polynomial.tex");
      boost::math::tools::latex_printer printer(fs);
      // Document preamble:
      printer << "\\documentclass{article}" << std::endl;
      printer << "\\begin{document}" << std::endl << std::endl;
      //
      // Stream out some manipulators to set the formatting:
      //
      printer << std::setprecision(2) << std::fixed << boost::math::tools::slanted_i;
      printer << poly << std::endl;
      printer << "\\end{document}\n\n";

      //]
   }
   return 0;
}
