//  (C) Copyright John Maddock 2006.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifdef _MSC_VER
#  pragma warning(disable: 4127) // conditional expression is constant.
#  pragma warning(disable: 4512) // assignment operator could not be generated.
#  pragma warning(disable: 4996) // 'std::char_traits<char>::copy' was declared deprecated.
#endif


#include <boost/math/concepts/real_concept.hpp>
#include <boost/test/included/test_exec_monitor.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/math/special_functions/factorials.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/tools/stats.hpp>
#include <boost/math/tools/test.hpp>

template <class T>
void test_spots(T)
{
   //
   // Basic sanity checks, tolerance is 10 decimal places expressed as a percentage,
   // One check per domain of the implementation:
   //
   T tolerance = boost::math::tools::epsilon<T>() * 100 * 2;  // 2 eps as a percent.
   BOOST_CHECK_CLOSE(
      ::boost::math::factorial<T>(0),
      static_cast<T>(1), tolerance);
   BOOST_CHECK_CLOSE(
      ::boost::math::factorial<T>(1),
      static_cast<T>(1), tolerance);
   BOOST_CHECK_CLOSE(
      ::boost::math::factorial<T>(10),
      static_cast<T>(3628800L), tolerance);
   BOOST_CHECK_CLOSE(
      ::boost::math::unchecked_factorial<T>(0),
      static_cast<T>(1), tolerance);
   BOOST_CHECK_CLOSE(
      ::boost::math::unchecked_factorial<T>(1),
      static_cast<T>(1), tolerance);
   BOOST_CHECK_CLOSE(
      ::boost::math::unchecked_factorial<T>(10),
      static_cast<T>(3628800L), tolerance);

   tolerance = boost::math::tools::epsilon<T>() * 100 * 20;  // 20 eps as a percent.
   unsigned i = boost::math::max_factorial<T>::value;
   BOOST_CHECK_CLOSE(
      ::boost::math::unchecked_factorial<T>(i),
      boost::math::tgamma(static_cast<T>(i+1)), tolerance);

   i += 10;
   while(boost::math::lgamma(static_cast<T>(i+1)) < boost::math::tools::log_max_value<T>())
   {
      BOOST_CHECK_CLOSE(
         ::boost::math::factorial<T>(i),
         boost::math::tgamma(static_cast<T>(i+1)), tolerance);
      i += 10;
   }

}

int test_main(int, char* [])
{
   test_spots(0.0F);
   test_spots(0.0);
   test_spots(0.0L);
   test_spots(boost::math::concepts::real_concept(0.1));
   return 0;
}

/*

Output is:

------ Build started: Project: test_factorials, Configuration: Debug Win32 ------
Compiling...
test_factorials.cpp
Linking...
Autorun "i:\boost-06-05-03-1300\libs\math\test\Math_test\debug\test_factorials.exe"
Running 1 test case...
*** No errors detected
Build Time 0:05
Build log was saved at "file://i:\boost-06-05-03-1300\libs\math\test\Math_test\test_factorials\Debug\BuildLog.htm"
test_factorials - 0 error(s), 0 warning(s)
========== Build: 1 succeeded, 0 failed, 0 up-to-date, 0 skipped ==========

*/



