//  (C) Copyright John Maddock 2019.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <boost/math/tools/formatting.hpp>
#include <complex>
#include <fstream>

void print(std::ostream& os)
{
   boost::math::tools::text_printer printer(os);

   printer << "Integers:\n\n";

   int ival = 0;
   printer << std::setw(20) << "Value" << std::setw(20) << "Base" << std::setw(20) << "Result" << std::endl;
   printer.stream() << std::setw(20) << ival << std::setw(20) << "default" << std::setw(20);
   printer << ival << std::endl;
   ival = -23;
   printer.stream() << std::setw(20) << ival << std::setw(20) << "default" << std::setw(20);
   printer << ival << std::endl;
   ival = 23;
   printer.stream() << std::setw(20) << ival << std::setw(20) << "hex" << std::setw(20);
   printer << std::hex << ival << std::endl;
   printer.stream() << std::setw(20) << ival << std::setw(20) << "oct" << std::setw(20);
   printer << std::oct << ival << std::endl;

   printer << std::endl << std::endl;

   printer << "Basic floating point values:\n\n";
   double fval = 3;
   printer << std::setw(20) << "Value" << std::setw(20) << "Precision" << std::setw(20) << "Format" << std::setw(20) << "Result" << std::endl;
   printer.stream() << std::setw(20) << fval << std::setw(20) << "default" << std::setw(20) << "default" << std::setw(20);
   printer << fval << std::endl;
   fval = 3.14;
   printer.stream() << std::setw(20) << fval << std::setw(20) << "default" << std::setw(20) << "default" << std::setw(20);
   printer << fval << std::endl;
   fval = -1.2345678765e-24;
   printer.stream() << std::setw(20) << fval << std::setw(20) << "default" << std::setw(20) << "default" << std::setw(20);
   printer << fval << std::endl;
   printer.stream() << std::setw(20) << fval << std::setw(20) << "3" << std::setw(20) << "scientific" << std::setw(20);
   printer << std::scientific << std::setprecision(3) << fval << std::endl;

   printer << std::endl << std::endl << std::defaultfloat;

   printer << "Complex Values:\n\n";
   std::complex<double> cval(3.25, 4.67);
   printer << std::setw(20) << "Value" << std::setw(20) << "Precision" << std::setw(20) << "Format" << std::setw(30) << "Result" << std::endl;
   printer.stream() << std::setw(20) << cval << std::setw(20) << "default" << std::setw(20) << "default" << std::setw(30);
   printer << cval << std::endl;
   cval = 3.14;
   printer.stream() << std::setw(20) << cval << std::setw(20) << "default" << std::setw(20) << "default" << std::setw(30);
   printer << cval << std::endl;
   cval = std::complex<double>(1.23, -1.2345678765e-24);
   printer.stream() << std::setw(20) << cval << std::setw(20) << "default" << std::setw(20) << "default" << std::setw(30);
   printer << cval << std::endl;
   printer.stream() << std::setw(20) << cval << std::setw(20) << "3" << std::setw(20) << "scientific" << std::setw(30);
   printer << std::scientific << std::setprecision(3) << cval << std::endl;
}

int main(int argc, const char* argv[])
{
   if (argc > 1)
   {
      std::ofstream ofs(argv[1]);
      print(ofs);
   }
   else
      print(std::cout);
   return 0;
}
