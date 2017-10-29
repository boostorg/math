//  (C) Copyright John Maddock 2006.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <pch_light.hpp>
#include "test_1F0.hpp"

#include <boost/multiprecision/cpp_bin_float.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>

void expected_results()
{
   //
   // Define the max and mean errors expected for
   // various compilers and platforms.
   //
   const char* largest_type;
#ifndef BOOST_MATH_NO_LONG_DOUBLE_MATH_FUNCTIONS
   if(boost::math::policies::digits<double, boost::math::policies::policy<> >() == boost::math::policies::digits<long double, boost::math::policies::policy<> >())
   {
      largest_type = "(long\\s+)?double";
   }
   else
   {
      largest_type = "long double";
   }
#else
   largest_type = "(long\\s+)?double";
#endif

   //
   // Finish off by printing out the compiler/stdlib/platform names,
   // we do this to make it easier to mark up expected error rates.
   //
   std::cout << "Tests run with " << BOOST_COMPILER << ", " 
      << BOOST_STDLIB << ", " << BOOST_PLATFORM << std::endl;
}

BOOST_AUTO_TEST_CASE( test_main )
{
   expected_results();
   BOOST_MATH_CONTROL_FP;

#ifndef BOOST_MATH_BUGGY_LARGE_FLOAT_CONSTANTS
   test_spots(0.0F);
#endif
   test_spots(0.0);
#ifndef BOOST_MATH_NO_LONG_DOUBLE_MATH_FUNCTIONS
   test_spots(0.0L);
#ifndef BOOST_MATH_NO_REAL_CONCEPT_TESTS
   test_spots(boost::math::concepts::real_concept(0.1));
#endif
#endif
   test_spots(boost::multiprecision::cpp_bin_float_quad());
   test_spots(boost::multiprecision::cpp_dec_float_50());
}



