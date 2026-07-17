//  (C) Copyright Matt Borland 2026.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// Module build check for lambert_w. The lookup tables in
// detail/lambert_w_lookup_table.ipp are namespace-scope arrays that are indexed
// at run time (an odr-use), so if they were left with internal linkage a module
// consumer could not resolve them. Exercising lambert_w0 and lambert_wm1 here
// forces that resolution when the test is built against the boost.math module.

#ifndef BOOST_MATH_BUILD_MODULE
#include <boost/math/special_functions/lambert_w.hpp>
#else
import boost.math;
#endif

#include <cmath>
#include "math_unit_test.hpp"

template <class T>
void test_spots(T, const char* t)
{
   std::cout << "Testing lambert_w for type " << t << std::endl;

   using std::exp;
   const T tolerance = 16;

   // W0 at 1 is the omega constant.
   CHECK_ULP_CLOSE(::boost::math::lambert_w0(static_cast<T>(1)),
      static_cast<T>(0.56714329040978387299996866221035554975381578718651L), tolerance);

   // W0 identity: w * e^w == z. Runtime indexing here reaches the w0 tables.
   for (const T z : {static_cast<T>(0.5), static_cast<T>(2), static_cast<T>(10),
                     static_cast<T>(100), static_cast<T>(1000)})
   {
      const T w = ::boost::math::lambert_w0(z);
      CHECK_ULP_CLOSE(w * exp(w), z, tolerance);
   }

   // W-1 identity on its branch (z in [-1/e, 0)); reaches the wm1 tables.
   for (const T z : {static_cast<T>(-0.3), static_cast<T>(-0.1),
                     static_cast<T>(-0.01), static_cast<T>(-0.001)})
   {
      const T w = ::boost::math::lambert_wm1(z);
      CHECK_ULP_CLOSE(w * exp(w), z, tolerance);
   }
}

int main()
{
   test_spots(0.0F, "float");
   test_spots(0.0, "double");

   return boost::math::test::report_errors();
}
