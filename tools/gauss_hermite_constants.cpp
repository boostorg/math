//  Copyright Jacob Hass 2026.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// Prints the tables of abscissas and weights in boost/math/quadrature/gauss_hermite.hpp, or with the argument
// laguerre, in boost/math/quadrature/gauss_laguerre.hpp.
// Category 999 forces the values to be computed on demand, at the precision of mp_type.

#include <iomanip>
#include <iostream>
#include <string>
#include <cstring>
#include <boost/math/quadrature/gauss_hermite.hpp>
#include <boost/math/quadrature/gauss_laguerre.hpp>
#include <boost/multiprecision/cpp_bin_float.hpp>

template <class T>
void print_gauss_hermite_constants(const char* suffix, int prec, int tag, const char* family)
{
   const std::string suffix_append = std::string(suffix) + ")";
   const char* type = prec > 40 ? "T" : "storage_type";
   const auto& ab = T::abscissa();
   const auto& w = T::weights();
   std::cout << std::setprecision(prec) << std::scientific;
   // The Hermite tables hold the nonnegative abscissas only.
   std::size_t order = std::strcmp(family, "hermite") ? ab.size() : (ab[0] == 0) ? (ab.size() * 2) - 1 : ab.size() * 2;
   std::cout <<
      "template <class T>\n"
      "class gauss_" << family << "_detail<T, " << order << ", " << tag << ">\n"
      "{\n"
      "   using storage_type = typename quadrature_constant_category<T>::storage_type;\n"
      "   public:\n"
      "      static std::array<" << type << ", " << ab.size() << "> const & abscissa()\n"
      "      {\n"
      "         static std::array<" << type << ", " << ab.size() << "> data = {\n";
   for (unsigned i = 0; i < ab.size(); ++i)
      std::cout << "            " << (prec > 40 ? "BOOST_MATH_HUGE_CONSTANT(T, 0, " : "static_cast<storage_type>(") << ab[i] << (prec > 40 ? ")" : suffix_append) << ",\n";
   std::cout <<
      "         };\n"
      "         return data;\n"
      "      }\n"
      "      static std::array<" << type << ", " << w.size() << "> const & weights()\n"
      "      {\n"
      "         static std::array<" << type << ", " << w.size() << "> data = {\n";
   for (unsigned i = 0; i < w.size(); ++i)
      std::cout << "            " << (prec > 40 ? "BOOST_MATH_HUGE_CONSTANT(T, 0, " : "static_cast<storage_type>(") << w[i] << (prec > 40 ? ")" : suffix_append) << ",\n";
   std::cout <<
      "         };\n"
      "         return data;\n"
      "      }\n"
      "   };\n\n";
}

using mp_type = boost::multiprecision::number<boost::multiprecision::cpp_bin_float<150> >;

template <class rule>
void print_all(const char* family)
{
   std::cout << "#ifndef BOOST_HAS_FLOAT128\n";
   print_gauss_hermite_constants<rule>("L", 35, 0, family);
   std::cout << "#else\n";
   print_gauss_hermite_constants<rule>("Q", 35, 0, family);
   std::cout << "#endif\n";
   print_gauss_hermite_constants<rule>("", 115, 4, family);
}

int main(int argc, char* argv[])
{
   using namespace boost::math::quadrature::detail;
   if (argc > 1 && !std::strcmp(argv[1], "laguerre"))
   {
      print_all<gauss_laguerre_detail<mp_type, 7, 999> >("laguerre");
      print_all<gauss_laguerre_detail<mp_type, 10, 999> >("laguerre");
      print_all<gauss_laguerre_detail<mp_type, 15, 999> >("laguerre");
      return 0;
   }
   print_all<gauss_hermite_detail<mp_type, 7, 999> >("hermite");
   print_all<gauss_hermite_detail<mp_type, 10, 999> >("hermite");
   print_all<gauss_hermite_detail<mp_type, 15, 999> >("hermite");
   return 0;
}
