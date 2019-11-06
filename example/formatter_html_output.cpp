//  (C) Copyright John Maddock 2019.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <boost/math/tools/formatting.hpp>
#include <complex>
#include <fstream>

void print(std::ostream& os)
{
   boost::math::tools::html_printer printer(os);

   printer << "<!DOCTYPE html>\n<html>\n<body>" << std::endl;
   printer << "<style>  span.number { font-weight:bold; } span.integer { color:red; } span.float { color:blue; } span.complex { color:green; } td {  text-align:right; } th { text-align:right; } </style>" << std::endl;

   printer << "<h2>Basic Integers Values:</h2>\n" << boost::math::tools::full_styling;

   int ival = 0;
   printer << "<table style=\"width:90%\"><tr><th>Value</th><th>Base</th><th>Result</th></tr>";
   printer.stream() << "<tr><td>" << ival << "</td><td>default</td><td>";
   printer << ival << "</td></tr>\n";
   ival = -23;
   printer.stream() << "<tr><td>" << ival << "</td><td>default</td><td>";
   printer << ival << "</td></tr>\n";
   ival = 23;
   printer.stream() << std::dec << "<tr><td>" << ival << "</td><td>hex</td><td>";
   printer << std::hex << ival << "</td></tr>\n";
   printer.stream() << std::dec << "<tr><td>" << ival << "</td><td>oct</td><td>";
   printer << std::oct << ival << "</td></tr>\n</table>\n\n";

   printer << "<h2>Basic Floating Point Values:</h2>\n\n";
   double fval = 3;
   printer << "<table style=\"width:90%\"><tr><th>Value</th><th>Precision</th><th>Format</th><th>Result</th></tr>";
   printer.stream() << "<tr><td>" << fval << "</td><td>default</td><td>default</td><td>";
   printer << fval << "</td></tr>\n";
   fval = 3.14;
   printer.stream() << "<tr><td>" << fval << "</td><td>default</td><td>default</td><td>";
   printer << fval << "</td></tr>\n";
   fval = -1.2345678765e-24;
   printer.stream() << "<tr><td>" << fval << "</td><td>default</td><td>default</td><td>";
   printer << fval << "</td></tr>\n";
   printer.stream() << "<tr><td>" << fval << "</td><td>3</td><td>scientific</td><td>";
   printer << std::setprecision(3) << std::scientific << fval << "</td></tr>\n";
   printer.stream() << "<tr><td>" << fval << "</td><td>default</td><td>multiply_x</td><td>";
   printer << boost::math::tools::multiply_x << fval << "</td></tr>\n";
   printer.stream() << "<tr><td>" << fval << "</td><td>3</td><td>scientific + multiply_x</td><td>";
   printer << std::setprecision(3) << std::scientific << fval << "</td></tr>\n";
   printer.stream() << "<tr><td>" << fval << "</td><td>default</td><td>multiply_dot</td><td>";
   printer << boost::math::tools::multiply_dot << fval << "</td></tr>\n";
   printer.stream() << "<tr><td>" << fval << "</td><td>3</td><td>scientific + multiply_dot</td><td>";
   printer << std::setprecision(3) << std::scientific << fval << "</td></tr>\n" << boost::math::tools::multiply_times;

   printer << "</table>\n\n" << std::defaultfloat;

   printer << "<h2>Complex Values:</h2>\n\n";
   std::complex<double> cval(3.25, 4.67);
   printer << "<table style=\"width:90%\"><tr><th>Value</th><th>Precision</th><th>Format</th><th>Result</th></tr>";
   printer.stream() << "<tr><td>" << cval << "</td><td>default</td><td>default</td><td>";
   printer << cval << "</td></tr>\n";
   cval = 3.14;
   printer.stream() << "<tr><td>" << cval << "</td><td>default</td><td>default</td><td>";
   printer << cval << "</td></tr>\n";
   cval = std::complex<double>(1.23, -1.2345678765e-24);
   printer.stream() << "<tr><td>" << cval << "</td><td>default</td><td>default</td><td>";
   printer << cval << "</td></tr>\n";
   printer.stream() << "<tr><td>" << cval << "</td><td>3</td><td>scientific</td><td>";
   printer << std::scientific << std::setprecision(3) << cval << "</td></tr>\n";
   printer.stream() << "<tr><td>" << cval << "</td><td>12</td><td>default + slanted_i</td><td>";
   printer << boost::math::tools::slanted_i << std::defaultfloat << std::setprecision(12) << cval << "</td></tr>\n";
   printer.stream() << "<tr><td>" << cval << "</td><td>12</td><td>default + doublestruck_i</td><td>";
   printer << boost::math::tools::doublestruck_i << std::defaultfloat << std::setprecision(12) << cval << "</td></tr>\n";

   printer.stream() << "<tr><td>" << cval << "</td><td>12</td><td>default + doublestruck_i + multiply_x</td><td>";
   printer << boost::math::tools::doublestruck_i << boost::math::tools::multiply_x << std::defaultfloat << std::setprecision(12) << cval << "</td></tr>\n";
   printer.stream() << "<tr><td>" << cval << "</td><td>12</td><td>default + doublestruck_i + multiply_dot</td><td>";
   printer << boost::math::tools::doublestruck_i << std::defaultfloat << std::setprecision(12) << boost::math::tools::multiply_dot << cval << "</td></tr>\n" << boost::math::tools::multiply_times;
   printer << "\n</table>\n\n";

   printer << "</body></html>\n\n";
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
