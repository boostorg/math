//  (C) Copyright John Maddock 2006.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include "test_1F1_log.hpp"

#include <boost/multiprecision/cpp_bin_float.hpp>

// The peak relative errors, in epsilon, of each type for the large random values;
// float is evaluated in double, so it is exact:
std::uintmax_t max_error(const std::string& type_name, const std::string&)
{
   return type_name == "float" ? 1 : type_name == "cpp_bin_float_quad" ? 500000 : type_name == "dec_40" ? 20000000 : 400000;
}

int main()
{
   std::cout << "Tests run with " << BOOST_COMPILER << ", " << BOOST_STDLIB << ", " << BOOST_PLATFORM << std::endl;

#ifndef BOOST_MATH_BUGGY_LARGE_FLOAT_CONSTANTS
#if !defined(TEST) || (TEST == 2)
   test_spots(0.0F, "float");
#endif
#endif
#if !defined(TEST) || (TEST == 3)
   test_spots(0.0, "double");
#endif
#ifndef BOOST_MATH_NO_LONG_DOUBLE_MATH_FUNCTIONS
#if (!defined(TEST) || (TEST == 4)) && (DBL_MAX_EXP != LDBL_MAX_EXP)
   test_spots(0.0L, "long double");
#endif
#ifndef BOOST_MATH_NO_REAL_CONCEPT_TESTS
#if !defined(TEST) || (TEST == 5)
   test_spots(boost::math::concepts::real_concept(0.1), "real_concept");
#endif
#endif
#endif
#if (!defined(TEST) || (TEST == 6)) && !defined(__clang__)
   #if defined(BOOST_MATH_RUN_MP_TESTS)
   test_spots(boost::multiprecision::cpp_bin_float_quad(), "cpp_bin_float_quad");
   #endif
#endif
#if (!defined(TEST) || (TEST == 7)) && !defined(__clang__)
   #if defined(BOOST_MATH_RUN_MP_TESTS)
   typedef boost::multiprecision::number<boost::multiprecision::cpp_bin_float<40> > dec_40;
   test_spots(dec_40(), "dec_40");
   #endif
#endif
   return boost::math::test::report_errors();
}
