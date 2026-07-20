//  (C) Copyright Matt Borland 2026.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// Module build check for ccmath::nextafter. Its detail dispatch defaults a
// template parameter to boost::math::tools::detail::has_backend_type_v<T>, so
// building this against the boost.math module exercises the variable templates
// in tools/traits.hpp from a module consumer.

#ifndef BOOST_MATH_BUILD_MODULE
#include <boost/math/ccmath/next.hpp>
#else
import boost.math;
#endif

#ifndef BOOST_MATH_BUILD_MODULE
#include <cmath>
#endif
#include "math_unit_test.hpp"

template <class T>
void test_spots(T, const char* type_name)
{
   std::cout << "Testing ccmath next for type " << type_name << std::endl;

   CHECK_EQUAL(boost::math::ccmath::nextafter(static_cast<T>(1), static_cast<T>(2)), std::nextafter(static_cast<T>(1), static_cast<T>(2)));
   CHECK_EQUAL(boost::math::ccmath::nextafter(static_cast<T>(1), static_cast<T>(0)), std::nextafter(static_cast<T>(1), static_cast<T>(0)));
}

int main()
{
   test_spots(0.0F, "float");
   test_spots(0.0, "double");

   return boost::math::test::report_errors();
}
