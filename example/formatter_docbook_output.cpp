//  (C) Copyright John Maddock 2019.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <boost/math/tools/formatting.hpp>
#include <complex>
#include <fstream>

void print(std::ostream& os)
{
   boost::math::tools::docbook_printer printer(os);

   printer << "[template integer_formatting_examples[]\n";

   int ival = 0;
   printer << "[table:integer_fmt_examples Basic Integer Values\n[[Value][Base][Result]]\n";
   printer.stream() << "[[" << ival << "][default]['''";
   printer << ival << "''']]\n";
   ival = -23;
   printer.stream() << "[[" << ival << "][default]['''";
   printer << ival << "''']]\n";
   ival = 23;
   printer.stream() << std::dec << "[[" << ival << "][hex]['''";
   printer << std::hex << ival << "''']]\n";
   printer.stream() << std::dec << "[[" << ival << "][oct]['''";
   printer << std::oct << ival << "''']]\n]\n]\n";

   printer << "[template float_formatting_examples[]\n";
   double fval = 3;
   printer << "[table:float_fmt_examples Basic Float Values\n[[Value][Precision][Format][Result]]";
   printer.stream() << "[[" << fval << "][default][default]['''";
   printer << fval << "''']]\n";
   fval = 3.14;
   printer.stream() << "[[" << fval << "][default][default]['''";
   printer << fval << "''']]\n";
   fval = -1.2345678765e-24;
   printer.stream() << "[[" << fval << "][default][default]['''";
   printer << fval << "''']]\n";
   printer.stream() << "[[" << fval << "][3][scientific]['''";
   printer << std::setprecision(3) << std::scientific << fval << "''']]\n";
   printer.stream() << "[[" << fval << "][3][scientific + multiply_x]['''";
   printer << std::setprecision(3) << std::scientific << boost::math::tools::multiply_x << fval << "''']]\n";
   printer.stream() << "[[" << fval << "][3][scientific + multiply_dot]['''";
   printer << std::setprecision(3) << std::scientific << boost::math::tools::multiply_dot << fval << "''']]\n" << boost::math::tools::multiply_times;

   printer << "]\n]\n" << std::defaultfloat;

   printer << "[template complex_formatting_examples[]\n";
   std::complex<double> cval(3.25, 4.67);
   printer << "[table:complex_fmt_examples Basic Complex Values\n[[Value][Precision][Format][Result]]";
   printer.stream() << "[[" << cval << "][default][default]['''";
   printer << cval << "''']]\n";
   cval = 3.14;
   printer.stream() << "[[" << cval << "][default][default]['''";
   printer << cval << "''']]\n";
   cval = std::complex<double>(1.23, -1.2345678765e-24);
   printer.stream() << "[[" << cval << "][default][default]['''";
   printer << cval << "''']]\n";
   printer.stream() << "[[" << cval << "][3][scientific]['''";
   printer << std::scientific << std::setprecision(3) << cval << "''']]\n";
   printer.stream() << "[[" << cval << "][12][default + slanted_i]['''";
   printer << boost::math::tools::slanted_i << std::defaultfloat << std::setprecision(12) << cval << "''']]\n";
   printer.stream() << "[[" << cval << "][12][default + doublestruck_i]['''";
   printer << boost::math::tools::doublestruck_i << std::defaultfloat << std::setprecision(12) << cval << "''']]\n";
   printer.stream() << "[[" << cval << "][12][default + doublestruck_i + multiply_x]['''";
   printer << boost::math::tools::doublestruck_i << std::defaultfloat << std::setprecision(12) << boost::math::tools::multiply_x << cval << "''']]\n";
   printer.stream() << "[[" << cval << "][12][default + doublestruck_i + multiply_dot]['''";
   printer << boost::math::tools::doublestruck_i << std::defaultfloat << std::setprecision(12) << boost::math::tools::multiply_dot << cval << "''']]\n" << boost::math::tools::multiply_x;
   printer << "\n]\n]\n";
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
