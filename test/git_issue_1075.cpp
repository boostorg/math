//  Copyright John Maddock 2017.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <iostream>
#include <boost/math/quadrature/gauss_kronrod.hpp>
#include <boost/multiprecision/cpp_bin_float.hpp>

template <typename T>
auto f1(T t) {
   return exp(-t * t / 2);
};

template <class T>
void test()
{
   T error = 0;

   auto r = boost::math::quadrature::gauss_kronrod<T, 15>::integrate(&f1<T>, 0, 1, 0, 0, &error);
   std::cout << r << std::endl;
   r = boost::math::quadrature::gauss_kronrod<T, 21>::integrate(&f1<T>, 0, 1, 0, 0, &error);
   std::cout << r << std::endl;
   r = boost::math::quadrature::gauss_kronrod<T, 31>::integrate(&f1<T>, 0, 1, 0, 0, &error);
   std::cout << r << std::endl;
   r = boost::math::quadrature::gauss_kronrod<T, 41>::integrate(&f1<T>, 0, 1, 0, 0, &error);
   std::cout << r << std::endl;
   r = boost::math::quadrature::gauss_kronrod<T, 51>::integrate(&f1<T>, 0, 1, 0, 0, &error);
   std::cout << r << std::endl;
   r = boost::math::quadrature::gauss_kronrod<T, 61>::integrate(&f1<T>, 0, 1, 0, 0, &error);
   std::cout << r << std::endl;
}


int main()
{
   test<boost::multiprecision::cpp_bin_float_double>();
   test<boost::multiprecision::cpp_bin_float_double_extended>();
   test<boost::multiprecision::cpp_bin_float_quad>();

   return 0;
}
