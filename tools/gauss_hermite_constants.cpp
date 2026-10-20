//  Copyright Jacob Hass 2026.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// Prints the tables of abscissas and weights in boost/math/quadrature/gauss_hermite.hpp.
// Category 999 forces the values to be computed on demand, at the precision of mp_type.

#include <iomanip>
#include <iostream>
#include <string>
#include <boost/math/quadrature/gauss_hermite.hpp>
#include <boost/multiprecision/cpp_bin_float.hpp>

template <class T>
void print_gauss_hermite_constants(const char* suffix, int prec, int tag)
{
   const std::string suffix_append = std::string(suffix) + ")";
   const char* type = prec > 40 ? "T" : "storage_type";
   const auto& ab = T::abscissa();
   const auto& w = T::weights();
   std::cout << std::setprecision(prec) << std::scientific;
   std::size_t order = (ab[0] == 0) ? (ab.size() * 2) - 1 : ab.size() * 2;
   std::cout <<
      "template <class T>\n"
      "class gauss_hermite_detail<T, " << order << ", " << tag << ">\n"
      "{\n"
      "   using storage_type = typename gauss_hermite_constant_category<T>::storage_type;\n"
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

template <unsigned N>
void print_all(void)
{
   using mp_type = boost::multiprecision::number<boost::multiprecision::cpp_bin_float<150> >;
   using rule = boost::math::quadrature::detail::gauss_hermite_detail<mp_type, N, 999>;
   std::cout << "#ifndef BOOST_HAS_FLOAT128\n";
   print_gauss_hermite_constants<rule>("L", 35, 0);
   std::cout << "#else\n";
   print_gauss_hermite_constants<rule>("Q", 35, 0);
   std::cout << "#endif\n";
   print_gauss_hermite_constants<rule>("", 115, 4);
}

int main()
{
   print_all<7>();
   print_all<10>();
   print_all<15>();
   return 0;
}
