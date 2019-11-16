//  (C) Copyright John Maddock 2019.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <boost/math/tools/formatting.hpp>
#include <complex>
#include <fstream>
#include <boost/math/tools/polynomial.hpp>
#include <boost/rational.hpp>
#include <boost/multiprecision/cpp_int.hpp>

void print(std::ostream& os)
{
   boost::math::tools::html_printer printer(os);

   printer << "<!DOCTYPE html>\n<html>\n<body>" << std::endl;
   printer << "<style>  span.number { font-weight:bold; } span.integer { color:red; } span.float { color:blue; } span.complex { color:green; } span.polynomial { color:Olive; } span.rational {color:gray; } td {  text-align:right; } th { text-align:right; } </style>" << std::endl;

   printer << "<h2>Synopsis</h2>\n<para>Sample HTML output for various number types, these are then styled in bold, and various (somewhat untasteful!) colors just because we can, and to check that our markup generation is working correctly.</para>\n";

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
   printer << std::oct << ival << "</td></tr>\n</table>\n\n" << std::dec;

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
   printer << std::setprecision(3) << std::scientific << fval << "</td></tr>\n" << boost::math::tools::multiply_times << std::defaultfloat;

   // Infinities:
   fval = std::numeric_limits<double>::infinity();
   printer.stream() << "<tr><td>" << fval << "</td><td>default</td><td>default</td><td>";
   printer << fval << "</td></tr>\n";
   fval = -std::numeric_limits<double>::infinity();
   printer.stream() << "<tr><td>" << fval << "</td><td>default</td><td>default</td><td>";
   printer << fval << "</td></tr>\n";
   // NaN's:
   fval = std::numeric_limits<double>::quiet_NaN();
   printer.stream() << "<tr><td>" << fval << "</td><td>default</td><td>default</td><td>";
   printer << fval << "</td></tr>\n";
   // Zeros:
   fval = 0;
   printer.stream() << "<tr><td>" << fval << "</td><td>default</td><td>default</td><td>";
   printer << fval << "</td></tr>\n";
   fval = -fval;
   printer.stream() << "<tr><td>" << fval << "</td><td>default</td><td>default</td><td>";
   printer << fval << "</td></tr>\n";

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
   printer << boost::math::tools::doublestruck_i << std::defaultfloat << std::setprecision(12) << boost::math::tools::multiply_dot << cval << "</td></tr>\n" << boost::math::tools::multiply_times << boost::math::tools::upright_i;
   printer << "\n</table>\n\n";

   printer << "<h2>Complex Special Values:</h2>\n\n";
   printer << "<table style=\"width:90%\"><tr><th>Value</th><th>Precision</th><th>Format</th><th>Result</th></tr>";
   cval = std::complex<double>(0, 0);
   printer.stream() << "<tr><td>" << cval << "</td><td>default</td><td>default</td><td>";
   printer << cval << "</td></tr>\n";
   printer.stream() << "<tr><td>" << cval << "</td><td>default</td><td>show_zero_components</td><td>";
   printer << boost::math::tools::show_zero_components << cval << boost::math::tools::hide_zero_components << "</td></tr>\n";
   cval = std::complex<double>(3.14, 0);
   printer.stream() << "<tr><td>" << cval << "</td><td>default</td><td>default</td><td>";
   printer << cval << "</td></tr>\n";
   printer.stream() << "<tr><td>" << cval << "</td><td>default</td><td>show_zero_components</td><td>";
   printer << boost::math::tools::show_zero_components << cval << boost::math::tools::hide_zero_components << "</td></tr>\n";
   cval = std::complex<double>(-3.14, 0);
   printer.stream() << "<tr><td>" << cval << "</td><td>default</td><td>default</td><td>";
   printer << cval << "</td></tr>\n";
   printer.stream() << "<tr><td>" << cval << "</td><td>default</td><td>show_zero_components</td><td>";
   printer << boost::math::tools::show_zero_components << cval << boost::math::tools::hide_zero_components << "</td></tr>\n";
   cval = std::complex<double>(0, 25.5);
   printer.stream() << "<tr><td>" << cval << "</td><td>default</td><td>default</td><td>";
   printer << cval << "</td></tr>\n";
   printer.stream() << "<tr><td>" << cval << "</td><td>default</td><td>show_zero_components</td><td>";
   printer << boost::math::tools::show_zero_components << cval << boost::math::tools::hide_zero_components << "</td></tr>\n";
   cval = std::complex<double>(0, -25.5);
   printer.stream() << "<tr><td>" << cval << "</td><td>default</td><td>default</td><td>";
   printer << cval << "</td></tr>\n";
   printer.stream() << "<tr><td>" << cval << "</td><td>default</td><td>show_zero_components</td><td>";
   printer << boost::math::tools::show_zero_components << cval << boost::math::tools::hide_zero_components << "</td></tr>\n";
   cval = std::complex<double>(std::numeric_limits<double>::infinity(), 0);
   printer.stream() << "<tr><td>" << cval << "</td><td>default</td><td>default</td><td>";
   printer << cval << "</td></tr>\n";
   cval = std::complex<double>(-std::numeric_limits<double>::infinity(), 0);
   printer.stream() << "<tr><td>" << cval << "</td><td>default</td><td>default</td><td>";
   printer << cval << "</td></tr>\n";
   cval = std::complex<double>(25.5, std::numeric_limits<double>::infinity());
   printer.stream() << "<tr><td>" << cval << "</td><td>default</td><td>default</td><td>";
   printer << cval << "</td></tr>\n";
   cval = std::complex<double>(25.5, -std::numeric_limits<double>::infinity());
   printer.stream() << "<tr><td>" << cval << "</td><td>default</td><td>default</td><td>";
   printer << cval << "</td></tr>\n";

   cval = std::complex<double>(std::numeric_limits<double>::quiet_NaN(), 0);
   printer.stream() << "<tr><td>" << cval << "</td><td>default</td><td>default</td><td>";
   printer << cval << "</td></tr>\n";
   cval = std::complex<double>(-std::numeric_limits<double>::quiet_NaN(), 0);
   printer.stream() << "<tr><td>" << cval << "</td><td>default</td><td>default</td><td>";
   printer << cval << "</td></tr>\n";
   cval = std::complex<double>(25.5, std::numeric_limits<double>::quiet_NaN());
   printer.stream() << "<tr><td>" << cval << "</td><td>default</td><td>default</td><td>";
   printer << cval << "</td></tr>\n";
   cval = std::complex<double>(25.5, -std::numeric_limits<double>::quiet_NaN());
   printer.stream() << "<tr><td>" << cval << "</td><td>default</td><td>default</td><td>";
   printer << cval << "</td></tr>\n";

   cval = std::complex<double>(std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::infinity());
   printer.stream() << "<tr><td>" << cval << "</td><td>default</td><td>default</td><td>";
   printer << cval << "</td></tr>\n";
   cval = std::complex<double>(-std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::infinity());
   printer.stream() << "<tr><td>" << cval << "</td><td>default</td><td>default</td><td>";
   printer << cval << "</td></tr>\n";
   cval = std::complex<double>(std::numeric_limits<double>::infinity(), std::numeric_limits<double>::quiet_NaN());
   printer.stream() << "<tr><td>" << cval << "</td><td>default</td><td>default</td><td>";
   printer << cval << "</td></tr>\n";
   cval = std::complex<double>(std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::quiet_NaN());
   printer.stream() << "<tr><td>" << cval << "</td><td>default</td><td>default</td><td>";
   printer << cval << "</td></tr>\n";

   printer << "\n</table>\n\n";

   printer << "<h2>Rationals:</h2>\n\n";

   printer << "<table style=\"width:90%\"><tr><th>Value</th><th>Result</th></tr>";
   boost::rational<int> rat(1, 3);
   printer.stream() << "<tr><td>" << rat << "</td><td>";
   printer << rat << "</td></tr>\n";
   rat = -rat;
   printer.stream() << "<tr><td>" << rat;
   printer << "</td><td>" << rat << "</td></tr>\n";
   rat *= 345634;
   rat /= 565;
   printer.stream() << "<tr><td>" << rat;
   printer << "</td><td>" << rat << "</td></tr>\n";
   rat = 0;
   printer.stream() << "<tr><td>" << rat;
   printer << "</td><td>" << rat << "</td></tr>\n";
   rat = -23;
   printer.stream() << "<tr><td>" << rat;
   printer << "</td><td>" << rat << "</td></tr>\n";
   boost::multiprecision::cpp_rational rat2(1);
   for (unsigned i = 1; i < 20; i += 2)
   {
      rat2 *= i;
      rat2 /= i + 1;
   }
   printer.stream() << "<tr><td>" << rat2;
   printer << "</td><td>" << rat2 << "</td></tr>\n";
   printer << "\n</table>\n\n";

   printer << "<h2>Polynomials:</h2>\n\n";

   printer << "<table style=\"width:90%\"><tr><th>Type</th><th>Result</th></tr>";
   boost::math::tools::polynomial<int> poly1 = { 2, -3, 4, 5 };
   printer << "<tr><td>Integer</td><td>" << poly1 << "</td></tr>\n";
   boost::math::tools::polynomial<double> poly2 = { 2.4, -34.25, 4.2e-6, -5.34e-67 };
   printer << "<tr><td>Float</td><td>" << poly2 << "</td></tr>\n";
   boost::math::tools::polynomial<std::complex<double> > poly3 = { { 2.4, 3.25 }, {-34.25 }, { 0, 4.2e-6 }, { -5.34e-67, 4.65e-20 } };
   printer << "<tr><td>Complex</td><td>" << poly3 << "</td></tr>\n";
   boost::math::tools::polynomial<boost::rational<int>> poly4 = { {2, 3}, {-3, 23}, {4, 56}, {5, 32} };
   printer << "<tr><td>Polynomial</td><td>" << poly4 << "</td></tr>\n";
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
