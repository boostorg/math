//  (C) Copyright Matt Borland 2026.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// Module build check for students_t. find_degrees_of_freedom(t, p) reads the
// namespace-scope constant detail::df_hint_fallback; building this against the
// boost.math module confirms that constant is usable from a module consumer.

#ifndef BOOST_MATH_BUILD_MODULE
#include <boost/math/distributions/students_t.hpp>
#else
import boost.math;
#endif

#include "math_unit_test.hpp"

template <class T>
void test_spots(T, const char* type_name)
{
   std::cout << "Testing students_t find_degrees_of_freedom for type " << type_name << std::endl;

   using dist_t = boost::math::students_t_distribution<T>;

   // The (t, p) overload runs the Edgeworth warm-start that seeds the root
   // finder from detail::df_hint_fallback. Recovering df and feeding it back
   // through the cdf is self-checking.
   const T df {dist_t::find_degrees_of_freedom(static_cast<T>(2), static_cast<T>(0.95))};
   CHECK_GE(df, static_cast<T>(0));

   const dist_t d {df};
   CHECK_ABSOLUTE_ERROR(boost::math::cdf(d, static_cast<T>(2)), static_cast<T>(0.95), static_cast<T>(1e-4));
}

int main()
{
   test_spots(0.0, "double");

   return boost::math::test::report_errors();
}
