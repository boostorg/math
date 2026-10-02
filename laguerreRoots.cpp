#include <vector>
#include <cmath>
#include <iostream>
#include <boost/math/special_functions/laguerre.hpp>
#include <boost/multiprecision/cpp_bin_float.hpp>
int main()
{
   using cpp_bin_float_100 = boost::multiprecision::cpp_bin_float_100;
   unsigned n = 5; // Degree of the Laguerre polynomial
   std::vector<cpp_bin_float_100> roots = boost::math::laguerre_zeros<cpp_bin_float_100>(n);

   std::cout << "Roots of the Laguerre polynomial of degree " << n << ":\n";
   for (const auto& root : roots)
   {
      std::cout << std::setprecision(100) << root << "\n";
   }

   return 0;
}