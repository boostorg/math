//  (C) Copyright John Maddock 2006.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_STANDALONE
#include <boost/math/concepts/real_concept.hpp>
#endif
#include "math_unit_test.hpp"
#include <boost/math/tools/rational.hpp>
#include <iostream>

template <class T, class U>
void do_test_spots(T, U);

template <class T>
void test_spots(T t, const char* n)
{
   std::cout << "Testing basic sanity checks for type " << n << std::endl;
   do_test_spots(t, int(0));
   do_test_spots(t, unsigned(0));
#ifdef BOOST_HAS_LONG_LONG
   do_test_spots(t, (boost::ulong_long_type)(0));
#endif
   do_test_spots(t, float(0));
   do_test_spots(t, T(0));
}

int main()
{
   test_spots(0.0F, "float");
   test_spots(0.0, "double");
#ifndef BOOST_MATH_NO_LONG_DOUBLE_MATH_FUNCTIONS
   test_spots(0.0L, "long double");
#if !defined(BOOST_MATH_NO_REAL_CONCEPT_TESTS) && !defined(BOOST_MATH_STANDALONE)
   test_spots(boost::math::concepts::real_concept(0.1), "real_concept");
#endif
#else
   std::cout << "<note>The long double tests have been disabled on this platform "
      "either because the long double overloads of the usual math functions are "
      "not available at all, or because they are too inaccurate for these tests "
      "to pass.</note>" << std::endl;
#endif
   return boost::math::test::report_errors();
}


