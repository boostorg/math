//  (C) Copyright John Maddock 2019.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <boost/math/tools/formatting.hpp>
#include <boost/math/tools/polynomial.hpp>
#include <boost/rational.hpp>
#include <boost/multiprecision/cpp_int.hpp>
#include <complex>
#include <fstream>

template <class charT>
void print(std::basic_ostream<charT>& os)
{
   boost::math::tools::basic_numeric_printer<boost::math::tools::text_format, charT> printer(os);

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
   printer << std::oct << ival << std::endl << std::dec;

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
   printer.stream() << std::scientific << std::setw(20) << fval << std::setw(20) << "3" << std::setw(20) << "scientific" << std::setw(20);
   printer << std::setprecision(3) << fval << std::endl << std::defaultfloat;
   fval = 0;
   printer.stream() << std::setw(20) << fval << std::setw(20) << "default" << std::setw(20) << "default" << std::setw(20);
   printer << fval << std::endl;
   fval = -fval;
   printer.stream() << std::setw(20) << fval << std::setw(20) << "default" << std::setw(20) << "default" << std::setw(20);
   printer << fval << std::endl;
   fval = 0;
   printer.stream() << std::scientific << std::setw(20) << fval << std::setw(20) << "default" << std::setw(20) << "scientific" << std::setw(20);
   printer << fval << std::endl;
   fval = -fval;
   printer.stream() << std::setw(20) << fval << std::setw(20) << "default" << std::setw(20) << "scientific" << std::setw(20);
   printer << fval << std::endl << std::defaultfloat;
   fval = std::numeric_limits<double>::infinity();
   printer.stream() << std::setw(20) << fval << std::setw(20) << "default" << std::setw(20) << "default" << std::setw(20);
   printer << fval << std::endl;
   fval = -std::numeric_limits<double>::infinity();
   printer.stream() << std::setw(20) << fval << std::setw(20) << "default" << std::setw(20) << "default" << std::setw(20);
   printer << fval << std::endl;
   fval = std::numeric_limits<double>::quiet_NaN();
   printer.stream() << std::setw(20) << fval << std::setw(20) << "default" << std::setw(20) << "default" << std::setw(20);
   printer << fval << std::endl;

   printer << std::endl << std::endl;

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
   printer << std::scientific << std::setprecision(3) << cval << std::endl << std::endl << std::defaultfloat;
   printer.stream() << std::setw(20) << cval << std::setw(20) << "slanted_i" << std::setw(20) << "default" << std::setw(30);
   printer << boost::math::tools::slanted_i << cval << std::endl;
   printer.stream() << std::setw(20) << cval << std::setw(20) << "doublestruck_i" << std::setw(20) << "default" << std::setw(30);
   printer << boost::math::tools::doublestruck_i << cval << std::endl << boost::math::tools::upright_i;

   std::cout << std::endl << std::endl;

   printer << "Complex Zeros:\n\n";
   cval = std::complex<double>(0);
   printer << std::setw(20) << "Value" << std::setw(20) << "Precision" << std::setw(23) << "Format" << std::setw(30) << "Result" << std::endl;
   printer.stream() << std::setw(20) << cval << std::setw(20) << "default" << std::setw(23) << "default" << std::setw(30);
   printer << cval << std::endl;
   printer.stream() << std::setw(20) << cval << std::setw(20) << "default" << std::setw(23) << "show_zero_components" << std::setw(30);
   printer << boost::math::tools::show_zero_components << cval << std::endl << boost::math::tools::hide_zero_components;
   cval = 3.14;
   printer.stream() << std::setw(20) << cval << std::setw(20) << "default" << std::setw(23) << "default" << std::setw(30);
   printer << cval << std::endl;
   printer.stream() << std::setw(20) << cval << std::setw(20) << "default" << std::setw(23) << "show_zero_components" << std::setw(30);
   printer << boost::math::tools::show_zero_components << cval << std::endl << boost::math::tools::hide_zero_components;
   cval = -3.14;
   printer.stream() << std::setw(20) << cval << std::setw(20) << "default" << std::setw(23) << "default" << std::setw(30);
   printer << cval << std::endl;
   printer.stream() << std::setw(20) << cval << std::setw(20) << "default" << std::setw(23) << "show_zero_components" << std::setw(30);
   printer << boost::math::tools::show_zero_components << cval << std::endl << boost::math::tools::hide_zero_components;
   cval = std::complex<double>(0, 25.25);
   printer.stream() << std::setw(20) << cval << std::setw(20) << "default" << std::setw(23) << "default" << std::setw(30);
   printer << cval << std::endl;
   printer.stream() << std::setw(20) << cval << std::setw(20) << "default" << std::setw(23) << "show_zero_components" << std::setw(30);
   printer << boost::math::tools::show_zero_components << cval << std::endl << boost::math::tools::hide_zero_components;
   cval = std::complex<double>(0, -25.25);
   printer.stream() << std::setw(20) << cval << std::setw(20) << "default" << std::setw(23) << "default" << std::setw(30);
   printer << cval << std::endl;
   printer.stream() << std::setw(20) << cval << std::setw(20) << "default" << std::setw(23) << "show_zero_components" << std::setw(30);
   printer << boost::math::tools::show_zero_components << cval << std::endl << boost::math::tools::hide_zero_components;
   // Infinities:
   cval = std::complex<double>(std::numeric_limits<double>::infinity(), -25.25);
   printer.stream() << std::setw(20) << cval << std::setw(20) << "default" << std::setw(23) << "default" << std::setw(30);
   printer << cval << std::endl;
   cval = std::complex<double>(2, std::numeric_limits<double>::infinity());
   printer.stream() << std::setw(20) << cval << std::setw(20) << "default" << std::setw(23) << "default" << std::setw(30);
   printer << cval << std::endl;
   cval = std::complex<double>(-std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity());
   printer.stream() << std::setw(20) << cval << std::setw(20) << "default" << std::setw(23) << "default" << std::setw(30);
   printer << cval << std::endl;
   // NaN's:
   cval = std::complex<double>(std::numeric_limits<double>::quiet_NaN(), -25.25);
   printer.stream() << std::setw(20) << cval << std::setw(20) << "default" << std::setw(23) << "default" << std::setw(30);
   printer << cval << std::endl;
   cval = std::complex<double>(3, std::numeric_limits<double>::quiet_NaN());
   printer.stream() << std::setw(20) << cval << std::setw(20) << "default" << std::setw(23) << "default" << std::setw(30);
   printer << cval << std::endl;
   cval = std::complex<double>(std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::infinity());
   printer.stream() << std::setw(20) << cval << std::setw(20) << "default" << std::setw(23) << "default" << std::setw(30);
   printer << cval << std::endl;
   cval = std::complex<double>(std::numeric_limits<double>::infinity(), std::numeric_limits<double>::quiet_NaN());
   printer.stream() << std::setw(20) << cval << std::setw(20) << "default" << std::setw(23) << "default" << std::setw(30);
   printer << cval << std::endl;


   printer << "\nRationals:\n\n";
   boost::rational<int> rat(1, 3);
   printer << std::setw(20) << "Value" << std::setw(20) << "Result" << std::endl;
   printer.stream() << std::setw(20) << rat;
   printer << std::setw(20) << rat << std::endl;
   rat = -rat;
   printer.stream() << std::setw(20) << rat;
   printer << std::setw(20) << rat << std::endl;
   rat *= 345634;
   rat /= 565;
   printer.stream() << std::setw(20) << rat;
   printer << std::setw(20) << rat << std::endl;
   boost::multiprecision::cpp_rational rat2(1);
   for (unsigned i = 1; i < 20; i += 2)
   {
      rat2 *= i;
      rat2 /= i + 1;
   }
   printer.stream() << std::setw(20) << rat2;
   printer << std::setw(20) << rat2 << std::endl;

   printer << "\nPolynomials:\n\n";
   boost::math::tools::polynomial<int> poly1 = { 2, -3, 4, 5 };
   printer << "Integer: " << poly1 << std::endl;
   boost::math::tools::polynomial<double> poly2 = { 2.4, -34.25, 4.2e-6, -5.34e-67 };
   printer << "Float: " << poly2 << std::endl;
   boost::math::tools::polynomial<std::complex<double> > poly3 = { { 2.4, 3.25 }, {-34.25 }, { 0, 4.2e-6 }, { -5.34e-67, 4.65e-20 } };
   printer << "Complex: " << poly3 << std::endl;
   boost::math::tools::polynomial<boost::rational<int>> poly4 = { {2, 3}, {-3, 23}, {4, 56}, {5, 32} };
   printer << "Rational: " << poly4 << std::endl;
}

#ifdef _WIN32
#include <Windows.h>
#endif

int main(int argc, const char* argv[])
{
   if (argc > 1)
   {
      std::ofstream ofs(argv[1]);
      print(ofs);
   }
   else
   {
#ifdef _WIN32
      SetConsoleOutputCP(CP_UTF8);
#endif
      print(std::cout);
   }
   return 0;
}
